import time
import os
import sys

import tensorflow as tf
import numpy as np
from matplotlib import pyplot as plt

from pinnschwarz.pde import PDE_1D_Steady_AdvecDiff
from pinnschwarz.pinn import PINN_Architecture, FD_1D_Steady, PINN_Schwarz_Steady

# Set data type
DTYPE = "float64"
tf.keras.backend.set_floatx(DTYPE)


class Logger:
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.log = open(logfile, "w")

    def write(self, output):
        self.terminal.write(output)
        self.log.write(output)
        self.log.flush()


def trainer(params, outdir, make_figs=False):
    # pull parameters for brevity, cast to proper type
    n_colloc = int(params["n_colloc"])
    domain_bounds = [float(param) for param in params["domain_bounds"]]
    percent_overlap = float(params["percent_overlap"])
    n_subdomains = int(params["n_subdomains"])
    Pe = float(params["peclet"])

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    BC_label = "_WDBC"
    BC_type = params["BC_type"]
    # Choose type of BC enforcement for system BCs (strong (1) or weak (0))
    if params["sys_bcs"] == "weak":
        sysBC = 0
    elif params["sys_bcs"] == "strong":
        sysBC = 1
    else:
        raise ValueError("'sys_bcs' must be either 'weak' or 'strong'")

    # Choose type of BC enforcement for Schwarz BCs (strong (1) or weak (0))
    if params["sch_bcs"] == "weak":
        schBC = 0
    elif params["sch_bcs"] == "strong":
        schBC = 1
    else:
        raise ValueError("'sch_bcs' must be either 'weak' or 'strong'")

    # assertion to catch SChwarz boundary enforcement when only domain is used
    if schBC:
        assert (
            n_subdomains > 1
        ), "Schwarz boundaries cannot be strongly enforced because you are using a single-domain model which has no Schwarz boundaries."

    # Store BC enforcement booleans together
    SDBC = [sysBC, schBC]
    if all(SDBC):
        BC_label = "_SDBC_both"
    elif SDBC[0]:
        BC_label = "_SDBC_sys"
    elif SDBC[1]:
        BC_label = "_SDBC_schwarz"

    # check to ensure mixed DBC is only associated with Dirichirlet-Dirichlet BC
    if (BC_label == "_SBC_sys" or BC_label == "_SDBC_schwarz") and not BC_type == "DD":
        raise Exception(
            "Mixed stong/weak BCs only compatible with Dirichlet-Dirichlet BC enformement. "
            "Change BC type to DD if running mixed DBCs"
        )
    if BC_label == "_SBC_both" and BC_type == "RR" and not n_subdomains == 2:
        raise Exception("Strong RBC only funcitonal for 2 subdomains")
    if BC_label == "_SBC_both" and BC_type == "DN":
        Warning("DNBC strong enforcement is not consistent")
    if percent_overlap == 0 and BC_type == "DD":
        raise Exception(
            "Dirchlet-Dirichlet must be used with overlapping domain, change percent overlap to be non-zero"
        )

    NN_label = "PINN"
    # number of snapshots per subdomain used to aid training, if using snapshots
    if bool(params["snap_loss"]):
        NN_label = "NN"
        snap = 2**6
    else:
        snap = 0

    # Set subdomains for FOM modeling, if any
    FOM_label = ""
    sub_FOM = np.zeros((n_subdomains,))
    if "fom_loc" in params:
        if params["fom_loc"] == "left":
            sub_FOM[0] = 1
            FOM_label = "_FOM_left"
        elif params["fom_loc"] == "right":
            sub_FOM[-1] = 1
            FOM_label = "_FOM_right"

    # Calculate the size of subdomains and overlap to achieve the parameters defined above
    domain_size = 1 / (n_subdomains * (1 - percent_overlap) + percent_overlap)
    overlap = percent_overlap * domain_size

    # Construct subdomain set
    sub = ()
    for i in range(n_subdomains):
        step = i * (domain_size - overlap)
        sub += ([step, step + domain_size],)

    # Truncate subdomain boundaries to 8 decimal precision
    sub = np.round(sub, 8)

    # Generate boundary points for each subdomain boundary
    X_b_om = [[tf.constant(np.tile([i], (1, 1)), dtype=DTYPE) for i in sub[j]] for j in range(len(sub))]
    lambda_xb = [tf.constant(np.tile(0, (1, 1)), dtype=DTYPE) for _ in range(len(sub))]

    # Set random seed for reproducible results
    tf.keras.utils.set_random_seed(0)

    # Declare nu based on Peclet number
    nu = 1 / Pe
    # Declare an instance of the PDE class
    pde1 = PDE_1D_Steady_AdvecDiff(nu=nu, order=int(params["order"]))

    # FOM value "true solution" to use as a reference
    x_true = tf.constant(
        np.linspace(domain_bounds[0], domain_bounds[1], num=n_colloc), shape=(n_colloc, 1), dtype=DTYPE
    )
    u_true = pde1.f(x_true)

    # Set number of FD points and initialize the step size
    n_FD = int(n_colloc / 2)

    # Initialize tuple to store internal points and models in each subdomain
    X_r_om = ()
    model_om = ()

    # Build and store neural networks for applicable subdomains; Store FD class for FOM sub-domains
    for i in range(sub.shape[0]):
        # If subdomain is to be modeled by NN, add random uniform points to subdomain, include sub-domain bounds
        if not sub_FOM[i]:
            xl = sub[i][0]
            xr = sub[i][1]
            temp = np.random.uniform(low=xl, high=xr, size=(n_colloc,))
            temp = np.append(temp, sub[i])
            X_r_om += (tf.constant(temp, shape=(temp.shape[0], 1), dtype=DTYPE),)

            model_om += (
                PINN_Architecture(
                    xl=xl,
                    xr=xr,
                    num_hidden_layers=int(params["n_layers"]),
                    num_neurons_per_layer=int(params["n_neurons"]),
                ),
            )
            model_om[i].build(input_shape=(None, X_r_om[i].shape[1]))

        else:
            # If subdomain is to be modeled by FD, construct uniform line space for FD on subdomain
            xl = sub[i][0]
            xr = sub[i][1]
            x_FD = np.linspace(xl, xr, num=n_FD)

            # Add line space for subdomain to internal point storage
            X_r_om += (tf.constant(x_FD, shape=(x_FD.shape[0], 1), dtype=DTYPE),)

            model = FD_1D_Steady(x_FD, domain_bounds, pde1)

            model_om += (model,)

    # Initialize schwarz loop operators
    schwarz_conv = 1
    ref_err = 1
    iterCount = 0
    np.random.seed(0)
    x_schwarz = [tf.constant(np.linspace(s[0], s[1], num=n_FD), shape=(n_FD, 1), dtype=DTYPE) for s in sub]
    u_i_minus1 = [tf.constant(np.random.rand(n_FD, 1), shape=(n_FD, 1), dtype=DTYPE) for _ in x_schwarz]
    u_i = [tf.constant(np.zeros((n_FD, 1)), shape=(n_FD, 1), dtype=DTYPE) for _ in x_schwarz]
    du_i = [tf.constant(np.zeros((n_FD, 1)), shape=(n_FD, 1), dtype=DTYPE) for _ in x_schwarz]

    # Initialize variables for plotting Schwarz results
    if make_figs:
        fig = plt.figure(layout="tight")
        fig.set_size_inches(6, 3)
        plt.xlabel("x", fontsize=14)
        plt.ylabel("u", fontsize=14)
        (ref,) = plt.plot(x_true, u_true, "k--")
        subdomain_plots = [plt.plot([], []) for _ in range(n_subdomains)]

    figdir = os.path.join(
        outdir,
        NN_label
        + BC_label
        + FOM_label
        + "_Pe_{:d}_nSub_{:d}_over{:f}".format(int(params["beta"] / nu), n_subdomains, percent_overlap),
    )
    if not os.path.isdir(figdir):
        os.mkdir(figdir)

    logger = Logger(os.path.join(figdir, "log.txt"))

    framecount = 1

    # Begin recording time
    time_start = time.time()

    # Record u_hat at each boundary to stack on each iteration
    u_gamma = np.zeros(sub.shape)
    du_gamma = np.zeros(sub.shape)

    theta = 1

    # Main Schwarz loop
    while (schwarz_conv > float(params["schwarz_tol"])) or (ref_err > float(params["err_tol"])):
        # initialize error check variables to 0 as they are to be added to during sub-domain iteration
        schwarz_conv = 0
        ref_err = 0

        # Add to Schwarz iteration count
        iterCount += 1
        L_count = -1

        # Update title for Schwarz iter
        if make_figs:
            plt.title("Schwarz iteration {:d}; Pe = {:.2f}".format(iterCount, Pe), fontsize=14)

        # loop over each model for training
        for s in range(len(model_om)):
            # Current model domain points
            X_r = X_r_om[s]
            # Current model boundary points
            X_b = X_b_om[s]

            # Current model
            model_r = model_om[s]
            # Adjacent models for interface conditions
            model_i = [model_om[s - 1 : s], model_om[s + 1 : s + 2]]

            # update the current models BCs according to adjacent models and save to model for SDBCs
            if any(SDBC):
                if BC_type == "DD":
                    # overlap
                    if model_i[0]:
                        u_gamma[s][0] = np.interp(sub[s][0], x_schwarz[s - 1][:, 0], u_i[s - 1][:, 0])
                        du_gamma[s][0] = np.interp(sub[s][0], x_schwarz[s - 1][:, 0], du_i[s - 1][:, 0])
                    if model_i[1]:
                        u_gamma[s][1] = np.interp(sub[s][1], x_schwarz[s + 1][:, 0], u_i[s + 1][:, 0])
                        du_gamma[s][1] = np.interp(sub[s][1], x_schwarz[s + 1][:, 0], du_i[s + 1][:, 0])
                else:
                    # non overlap
                    if model_i[0]:
                        u_gamma[s][0] = u_i[s - 1][-1, 0]
                        du_gamma[s][0] = du_i[s - 1][-1, 0]
                    if model_i[1]:
                        u_gamma[s][1] = u_i[s + 1][0, 0]
                        du_gamma[s][1] = du_i[s + 1][0, 0]
            model_r.u_gamma = u_gamma[s]
            model_r.du_gamma = du_gamma[s]

            if BC_type == "DN":
                if iterCount % 2 == 0 and model_i[0]:
                    L_count += 1
                    u_xb = model_i[0][0](X_b[0])
                    lambda_xb[L_count] = theta * u_xb + (1 - theta) * lambda_xb[L_count]
                elif not iterCount % 2 == 0 and model_i[1]:
                    L_count += 1
                    u_xb = model_i[1][0](X_b[1])
                    lambda_xb[L_count] = theta * u_xb + (1 - theta) * lambda_xb[L_count]

            # Initialize solver
            p = PINN_Schwarz_Steady(
                pde1,
                model_r,
                model_i,
                SDBC,
                X_r,
                X_b,
                float(params["alpha"]),
                snap,
                lambda_xb[L_count],
                BC_type,
                iterCount,
            )

            # Solve model for current sub-domain
            p.solve(tf.keras.optimizers.Adam(learning_rate=float(params["learn_rate"])), int(params["n_epochs"]))

            # If current model is FD, output FD Error, update current iteration solution and convergence metrics
            if isinstance(model_r, FD_1D_Steady):
                print("Model {:d}: ".format(s + 1))
                print("\t" + "Finite Difference error = {:10.8e}".format(p.err))

                u_i[s] = model_r(x_schwarz[s])
                schwarz_conv += tf.math.reduce_euclidean_norm(u_i[s] - u_i_minus1[s]) / tf.math.reduce_euclidean_norm(
                    u_i[s]
                )
                ref_err += tf.math.reduce_euclidean_norm(u_i[s] - pde1.f(x_schwarz[s])) / tf.math.reduce_euclidean_norm(
                    pde1.f(x_schwarz[s])
                )

            # If current model is NN output the model loss, update current iteration solution and convergence metrics
            else:
                print("Model {:d}: ".format(s + 1))
                print("\t" + "Residual loss = {:10.8e}".format(p.phi_r))
                if not schBC:
                    print("\t" + "Interface loss = {:10.8e}".format(p.phi_i))
                if not sysBC:
                    print("\t" + "Boundary loss = {:10.8e}".format(p.phi_b))
                if snap:
                    print("\t" + "Snapshot loss = {:10.8e}".format(p.phi_s))
                print("\t" + "Total loss = {:10.8e}".format(p.loss))

                if any(SDBC):
                    u_i[s] = p.get_u_hat(x_schwarz[s], model_r)
                else:
                    u_i[s] = model_r(x_schwarz[s])

                du_i[s], _ = p.get_gradients(model_r, x_schwarz[s])

                schwarz_conv += tf.math.reduce_euclidean_norm(u_i[s] - u_i_minus1[s]) / tf.math.reduce_euclidean_norm(
                    u_i[s]
                )
                ref_err += tf.math.reduce_euclidean_norm(u_i[s] - pde1.f(x_schwarz[s])) / tf.math.reduce_euclidean_norm(
                    pde1.f(x_schwarz[s])
                )

            # update frame for animation
            subdomain_plots[s][0].set_data(x_schwarz[s][:, 0], u_i[s])

            # update frame for animation
            if make_figs:
                subdomain_plots[s][0].set_data(x_schwarz[s][:, 0], u_i[s])

        # Calculate the normalized difference between u for the current iteration and u for the previous iteration
        schwarz_conv = schwarz_conv / len(u_i)
        ref_err = ref_err / len(u_i)

        # Update the value of u stored for the previous iteration
        u_i_minus1 = u_i.copy()

        # Save current frame as image if needed for beamer animations
        if make_figs:
            figfile = os.path.join(figdir, f"fig_{framecount}.png")
            fig.savefig(figfile, dpi=100)
            framecount += 1

        # Output current Schwarz error
        print("")
        logger.write(
            "Schwarz iteration {:d}: Convergence error = {:10.8e}, Reference Error = {:10.8e}\n".format(
                iterCount, schwarz_conv, ref_err
            )
        )

        # Cut off simulat at a particular number of Schwarz iterations
        if iterCount == 100:
            break

    # End time recording
    cpu_time = time.time() - time_start

    if make_figs:
        plt.close(fig)

    return cpu_time, iterCount, ref_err
