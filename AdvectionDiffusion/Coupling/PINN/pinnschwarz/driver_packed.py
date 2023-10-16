# This file is meant to be the cleaned driver.
######################IMPORT LIBRARIES#########################################
import time
import argparse
import os
import sys
import tensorflow as tf
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import yaml

from pinnschwarz.pde import PDE_1D_Steady_AdvecDiff
from pinnschwarz.pinn import PINN_Architecture, FD_1D_Steady, PINN_Schwarz_Steady
######################END OF IMPORT############################################


######################SET UP DATA AND SYSTEM###################################
# Set the system
os.environ["KMP_AFFINITY"] = "noverbose"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# Set data type
DTYPE = 'float64'
tf.keras.backend.set_floatx(DTYPE)
######################END OF DATA AND SYSTEM###################################

# This class is used to log information into the terminal
class Logger():
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.log = open(logfile, "w")

    def write(self, output):
        self.terminal.write(output)
        self.log.write(output)
        self.log.flush()

# This class drives the major computation
class Driver():
    
    #parameter_file is the csv file
    parameter_file = ''
    #output directory, specified when running
    outdir = ''
    #the yaml file
    hyper = {}
    #whether we want to plot figures
    make_fig = False
    #record CPU time
    start = time.time()
    end = time.time()
    
    #This specifies the mdoel type
    mod_type = ''
    
    #This contains the training result
    result = {}
    
    #This function initializes with the above parameters
    def __init__(self, parameter_file, outdir, hyper_path, make_fig, mod_type="PINN"):
        self.parameter_file = parameter_file
        self.outdir = outdir
        with open(hyper_path, 'r') as file:
            self.hyper = yaml.safe_load(file)
        self.make_fig = make_fig
        self.mod_type = mod_type
        
    # This function replace the original parameter sweep block.
    # It returns a dictionary (label_dic) containing the following:
    # - BC_label: used in ploting figures
    # - NN_label: used in ploting figures
    # - SDBC: dictionary used in main loop
    # - snap: used in main loop
    # - FOM_label: used in ploting figures
    # - n_subdomains: used in plot fig and domain_size
    # - percent_overlap: used in domain_size
    def param(self):
        
        # Initialize dataframe for parameters and results storage
        ParameterSweep = pd.read_csv(self.parameter_file)
        
        # The dictionary storing label values
        label_dic = {}
        
        #For now, just one row
        z = 0
        
        # Set percentage overlap
        percent_overlap = ParameterSweep.loc[z, 'Percent Overlap']

        # Set the number of subdomains and the desired percentage overlap
        n_subdomains = ParameterSweep.loc[z, 'Sub-domains']

        BC_label = '_WDBC'
        # Choose type of BC enforcement for system BCs (strong (1) or weak (0))
        if ParameterSweep.loc[z, 'System BCs'] == 'weak':
            sysBC = 0
        elif ParameterSweep.loc[z, 'System BCs'] == 'strong':
            sysBC = 1

        # Choose type of BC enforcement for Schwarz BCs (strong (1) or weak (0))
        if ParameterSweep.loc[z, 'Schwarz BCs'] == 'weak':
            schBC = 0
        elif ParameterSweep.loc[z, 'Schwarz BCs'] == 'strong':
            schBC = 1

        # assertion to catch SChwarz boundary enforcement when only domain is used
        if schBC:
            assert n_subdomains > 1, "Schwarz boundaries cannot be strongly enforced because you are using a single-domain model which has no Schwarz boundaries."

        # Store BC enforcement booleans together
        SDBC = [sysBC, schBC]
        if all(SDBC):
            BC_label = '_SDBC_both'
        elif SDBC[0]:
            BC_label = '_SDBC_sys'
        elif SDBC[1]:
            BC_label = '_SDBC_schwarz'

        NN_label = 'PINN'
        # number of snapshots per subdomain used to aid training, if using snapshots
        if ParameterSweep.loc[z, 'Snapshots']:
            NN_label = 'NN'
            snap = 2**6
        else:
            snap=0

        # Set subdomains for FOM modeling, if any
        FOM_label = ''
        sub_FOM = np.zeros((n_subdomains,))
        if ParameterSweep.loc[z, 'FOM'] == 'left':
            sub_FOM[0] = 1
            FOM_label = '_FOM_left'
        elif ParameterSweep.loc[z, 'FOM'] == 'right':
            sub_FOM[-1] = 1
            FOM_label = '_FOM_right'
        
        Pe = ParameterSweep.loc[z, 'Peclet Number']
        
        label_dic = {
            'BC_label': BC_label,
            'NN_label': NN_label,
            'FOM_label': FOM_label,
            "SDBC": SDBC,
            "snap": snap,
            "n_subdomains": n_subdomains,
            "percent_overlap": percent_overlap,
            "Pe": Pe
            }
        
        return label_dic, sub_FOM
    
    #This function initializes desired pde and returns the pde and X_b_om.
    def initialize_pde(self, label_dic):
        
        # Calculate the size of subdomains and overlap to achieve the parameters defined above
        domain_size =  1 / ( label_dic["n_subdomains"]*(1 - label_dic["percent_overlap"]) 
                            + label_dic["percent_overlap"] )
        overlap = label_dic["percent_overlap"]*label_dic["n_subdomains"]

        # Construct subdomain set
        sub = ()
        for i in range(label_dic["n_subdomains"]):
            step = i*(domain_size - overlap)
            sub += ([step, step+domain_size],)

        # Truncate subdomain boundaries to 8 decimal precision
        sub = np.round(sub, 8)

        # Generate boundary points for each subdomain boundary
        X_b_om = [ [ tf.constant(np.tile([i], (2**int(self.hyper['PDE']['b_power']), 1)), dtype=DTYPE) for i in sub[j] ] for j in range(len(sub))]

        # Set random seed for reproducible results
        tf.random.set_seed(0)
        np.random.seed(0)

        # Declare nu based on Peclet number
        nu = 1/label_dic["Pe"]
        # Declare an instance of the PDE class
        pde1 = PDE_1D_Steady_AdvecDiff(nu=nu, order=self.hyper['PDE']['order'])
        
        return pde1, X_b_om, sub
    
    #This function initializes NN models in each subdomain
    def initialize_nn(self, pde1, sub, sub_FOM):
        
        domain = self.hyper['PDE']['domain']
        N = 2**int(self.hyper['PDE']['power'])
        
        # FOM value "true solution" to use as a reference
        x_true = tf.constant(np.linspace(domain[0], domain[1], num=N), shape=(N, 1), dtype=DTYPE)
        #u_true = pde1.f(x_true)

        # Set number of FD points and initialize the step size
        n_FD = int(N/2)
        
        # Initialize tuple to store internal points and models in each subdomain
        X_r_om = ()
        model_om = ()
        
        hl = self.hyper['NN_para']['hl']
        nl = self.hyper['NN_para']['nl']

        # Build and store neural networks for applicable subdomains; Store FD class for FOM sub-domains
        for i in range(sub.shape[0]):

            # If subdomain is to be modeled by NN, add random uniform points to subdomain, include sub-domain bounds
            if not sub_FOM[i]:
                xl = sub[i][0]
                xr = sub[i][1]
                temp = np.random.uniform(low=xl, high=xr, size=(N,))
                temp = np.append(temp, sub[i])
                X_r_om += ( tf.constant(temp, shape=(temp.shape[0],1), dtype=DTYPE), )

                model_om += ( PINN_Architecture(xl=xl, xr=xr, num_hidden_layers=hl, num_neurons_per_layer=nl), )
                model_om[i].build(input_shape=(None, X_r_om[i].shape[1]))

            else:

                # If subdomain is to be modeled by FD, construct uniform line space for FD on subdomain
                xl = sub[i][0]
                xr = sub[i][1]
                x_FD = np.linspace(xl, xr, num=n_FD)

                # Add line space for subdomain to internal point storage
                X_r_om += ( tf.constant(x_FD, shape=(x_FD.shape[0],1), dtype=DTYPE), )

                model = FD_1D_Steady(x_FD, domain, pde1)

                model_om += ( model, )
                
        return X_r_om, model_om, x_true
        
        
    def train(self, x_true, u_true, label_dic, model_om, X_r_om, X_b_om, sub, pde1):
        
        # ----- START INPUTS FROM YAML -----

        # Initialize list of points which lie on the system boundaries
        beta = self.hyper['PDE']['beta']
    
        # Declare constant hyperparameters
        alpha = self.hyper['hyper']['alpha']
        numEpochs = 2**int(self.hyper['hyper']['epochs_power'])
        learn_rate = self.hyper['hyper']['learn_rate']
        schwarz_tol = self.hyper['hyper']['schwarz_tol']
        err_tol = self.hyper['hyper']['err_tol']
    
        # ----- END FIXED INPUTS -----
    
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
                
        n_FD = int(2**int(self.hyper['PDE']['power']) / 2)
        

        # Initialize schwarz loop operators
        schwarz_conv = 1
        ref_err = 1
        iterCount = 0
        np.random.seed(0)
        x_schwarz = [tf.constant(np.linspace(s[0], s[1], num=n_FD), shape=(n_FD, 1), dtype=DTYPE) for s in sub]
        u_i_minus1 = [tf.constant(np.random.rand(n_FD,1), shape=(n_FD, 1), dtype=DTYPE) for _ in x_schwarz]
        u_i = [tf.constant(np.zeros((n_FD,1)), shape=(n_FD, 1), dtype=DTYPE) for _ in x_schwarz]

        # Initialize variables for plotting Schwarz results
        if self.make_fig:
            fig = plt.figure(layout='tight')
            fig.set_size_inches(6, 3)
            plt.xlabel('x', fontsize=14)
            plt.ylabel('u', fontsize=14)
            ref, = plt.plot(x_true, u_true, 'k--')
            subdomain_plots = [plt.plot([], []) for _ in range(label_dic["n_subdomains"])]
    
        figdir = os.path.join(
            self.outdir,
            label_dic["NN_label"]+label_dic["BC_label"]+label_dic["FOM_label"]
            +"_Pe_{:d}_nSub_{:d}_over{:f}".format(int(beta/(1/label_dic["Pe"])), label_dic["n_subdomains"], label_dic["percent_overlap"])
        )
        
        if not os.path.isdir(figdir):
            os.mkdir(figdir)
            
        logger = Logger(os.path.join(figdir, "log.txt"))
        
        framecount = 1

        # Begin recording time
        self.start = time.time()

        # Record u_hat at each boundary to stack on each iteration
        u_gamma = np.zeros(sub.shape)

        # Main Schwarz loop
        while (schwarz_conv > schwarz_tol or ref_err > err_tol):

            # initialize error check variables to 0 as they are to be added to during sub-domain iteration
            schwarz_conv = 0
            ref_err = 0

            # Add to Schwarz iteration count
            iterCount += 1

            # Update title for Schwarz iter
            plt.title('Schwarz iteration {:d}; Pe = {:d}'.format(iterCount, label_dic["Pe"]), fontsize=14)

            # loop over each model for training
            for s in range(len(model_om)):

                # Current model domain points
                X_r = X_r_om[s]
                # Current model boundary points
                X_b = X_b_om[s]

                # Current model
                model_r = model_om[s]
                # Adjacent models for interface conditions
                model_i = [model_om[s-1:s], model_om[s+1:s+2]]

                # update the current models BCs according to adjacent models and save to model for SDBCs
                if any(label_dic["SDBC"]):
                    if model_i[0]:
                        u_gamma[s][0] = np.interp(sub[s][0], x_schwarz[s-1][:,0], u_i[s-1][:,0])
                    if model_i[1]:
                        u_gamma[s][1] = np.interp(sub[s][1], x_schwarz[s+1][:,0], u_i[s+1][:,0])
                model_r.u_gamma = u_gamma[s]

                # Initialize solver
                p = PINN_Schwarz_Steady(pde1, model_r, model_i, label_dic["SDBC"], X_r, X_b, alpha, label_dic["snap"])

                # Solve model for current sub-domain
                p.solve(tf.keras.optimizers.Adam(learning_rate=learn_rate), numEpochs)

                MSE_tag = "MSE Sub-domain {:d}".format(s+1)

                # If current model is FD, output FD Error, update current iteration solution and convergence metrics
                if isinstance(model_r, FD_1D_Steady):
                    print('Model {:d}: '.format(s+1))
                    print('\t'+'Finite Difference error = {:10.8e}'.format(p.err))

                    self.result[MSE_tag] = np.float64(p.err)

                    u_i[s] = model_r(x_schwarz[s])
                    schwarz_conv += tf.math.reduce_euclidean_norm(u_i[s] - u_i_minus1[s])/tf.math.reduce_euclidean_norm(u_i[s])
                    ref_err += tf.math.reduce_euclidean_norm(u_i[s] - pde1.f(x_schwarz[s]))/tf.math.reduce_euclidean_norm(pde1.f(x_schwarz[s]))

                # If current model is NN output the model loss, update current iteration solution and convergence metrics
                else:
                    print('Model {:d}: '.format(s+1))
                    print('\t'+'Residual loss = {:10.8e}'.format(p.phi_r))
                    if not label_dic["SDBC"][1]:
                        print('\t'+'Interface loss = {:10.8e}'.format(p.phi_i))
                    if not label_dic["SDBC"][0]:
                        print('\t'+'Boundary loss = {:10.8e}'.format(p.phi_b))
                    if label_dic["snap"]:
                        print('\t'+'Snapshot loss = {:10.8e}'.format(p.phi_s))
                    print('\t'+'Total loss = {:10.8e}'.format(p.loss))

                    self.result[MSE_tag] = np.float64(p.loss)

                    if any(label_dic["SDBC"]):
                        u_i[s] = p.get_u_hat(x_schwarz[s], model_r)
                    else:
                        u_i[s] = model_r(x_schwarz[s])

                    schwarz_conv += tf.math.reduce_euclidean_norm(u_i[s] - u_i_minus1[s])/tf.math.reduce_euclidean_norm(u_i[s])
                    ref_err += tf.math.reduce_euclidean_norm(u_i[s] - pde1.f(x_schwarz[s]))/tf.math.reduce_euclidean_norm(pde1.f(x_schwarz[s]))


                # update frame for animation
                if self.make_fig:
                    subdomain_plots[s][0].set_data(x_schwarz[s][:,0], u_i[s])


            # Calculate the normalized difference between u for the current iteration and u for the previous iteration
            schwarz_conv = schwarz_conv/len(u_i)
            ref_err = ref_err/len(u_i)

            # Update the value of u stored for the previous iteration
            u_i_minus1 = u_i.copy()

            # Save current frame as image if needed for beamer animations
            if self.make_fig:
                figfile = os.path.join(figdir, f"fig_{framecount}.png")
                fig.savefig(figfile, dpi=100)
                framecount += 1

            # Output current Schwarz error
            print("")
            logger.write('Schwarz iteration {:d}: Convergence error = {:10.8e}, Reference Error = {:10.8e}\n'.format(iterCount, schwarz_conv, ref_err))

            # Cut off simulat at a particular number of Schwarz iterations
            if (iterCount == 100):
                break

        # End time recording
        self.end = time.time()

        # Record results for coupled model
        self.result['CPU Time (s)'] = self.end - self.start
        self.result['Schwarz Iterations'] = iterCount
        self.result['Avg L2 Error'] = np.float64(ref_err)

        # Export final results to CSV
        #ParameterSweep.to_csv(os.path.join(self.outdir, "ParamResults.csv"))
        
        if self.make_fig:
            plt.close(fig)
            
        return self.end - self.start

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='PINN-Schwarz',
        description='Schwarz-based training of PINN-PINN coupling')

    parser.add_argument('parameter_file')
    parser.add_argument('outdir')
    args = parser.parse_args()

    train_mod = Driver(args.parameter_file, args.outdir, 'hyper.yaml', False)
    label_dic, sub_FOM = train_mod.param()
    pde1, X_b_om, sub = train_mod.initialize_pde(label_dic)
    X_r_om, model_om, x_true = train_mod.initialize_nn(pde1, sub, sub_FOM)
    
    time = train_mod.train(x_true, pde1.f(x_true), label_dic, model_om, X_r_om, X_b_om, sub, pde1)
    result = train_mod.result
    print(result)

