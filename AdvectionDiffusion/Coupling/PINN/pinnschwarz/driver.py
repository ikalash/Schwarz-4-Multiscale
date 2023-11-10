import argparse
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import yaml
from ray import tune

from pinnschwarz.trainer import trainer


def parse_param(param_dict, param_str, can_opt=False, default=None):

    if param_str in param_dict:
        param = param_dict[param_str]

        # define search space for parameter
        if isinstance(param, dict):
            assert can_opt, f"Must set optimizer params if setting search space for parameter: {param_str}"

            try:
                space_type = str(param["space_type"])
                space_defs = list(param["space_defs"])
            except KeyError:
                print(f"Must define a 'space_type' and 'space_defs' for optimized parameter {param_str}")
                raise

            # expand these as necessary
            if space_type == "uniform":
                param_out = tune.uniform(space_defs[0], space_defs[1])
            elif space_type == "quniform":
                param_out = tune.quniform(space_defs[0], space_defs[1], space_defs[2])
            elif space_type == "randn":
                param_out = tune.randn(space_defs[0], space_defs[1])
            elif space_type == "randint":
                param_out = tune.randint(space_defs[0], space_defs[1])
            elif space_type == "choice":
                param_out = tune.choice(space_defs)
            elif space_type == "grid":
                param_out = tune.grid_search(space_defs)
            else:
                raise ValueError(f"Invalid 'space_type' for parameter {param_str}: {space_type}")

        # single-value input
        else:
            param_out = param
            if default is not None:
                assert isinstance(param_out, type(default))

    # fall back to default, if provided
    elif default is not None:
        param_out = default

    else:
        raise ValueError(f"No input or default value for parameter: {param_str}")

    return param_out


def append_param_space(param_space, param_dict, param_str, can_opt, default=None):

    assert param_str not in param_space, f"Invalid double insertion of parameter {param_str} into parameter space"

    param_val = parse_param(param_dict, param_str, can_opt=can_opt, default=default)
    param_space[param_str] = param_val

    return param_space


def parse_input(param_file):

    with open(param_file, "r") as file:
        param_dict = yaml.safe_load(file)
    param_space = {}

    # optimization parameters
    if "optimizer" in param_dict.keys():
        opt_dict = param_dict["optimizer"]
    else:
        opt_dict = None

    # PDE parameters
    pde_dict = param_dict["pde"]
    param_space = append_param_space(param_space, pde_dict, "order", False, default=2)
    param_space = append_param_space(param_space, pde_dict, "n_colloc", False, default=1000)
    param_space = append_param_space(param_space, pde_dict, "domain_bounds", False, default=[0.0, 1.0])
    param_space = append_param_space(param_space, pde_dict, "beta", False, default=1.0)
    param_space = append_param_space(param_space, pde_dict, "peclet", False)

    # Schwarz parameters
    schwarz_dict = param_dict["schwarz"]
    param_space = append_param_space(param_space, schwarz_dict, "n_subdomains", True)
    param_space = append_param_space(param_space, schwarz_dict, "percent_overlap", True)
    param_space = append_param_space(param_space, schwarz_dict, "sys_bcs", True)
    param_space = append_param_space(param_space, schwarz_dict, "sch_bcs", True)
    param_space = append_param_space(param_space, schwarz_dict, "snap_loss", True)
    param_space = append_param_space(param_space, schwarz_dict, "schwarz_tol", False, 0.001)
    param_space = append_param_space(param_space, schwarz_dict, "err_tol", False, 0.005)

    # neural network parameters
    nn_dict = param_dict["nn"]
    param_space = append_param_space(param_space, nn_dict, "alpha", True, 0.2)
    param_space = append_param_space(param_space, nn_dict, "n_epochs", True, 1000)
    param_space = append_param_space(param_space, nn_dict, "learn_rate", True, 0.001)
    param_space = append_param_space(param_space, nn_dict, "n_layers", True, 2)
    param_space = append_param_space(param_space, nn_dict, "n_neurons", True, 20)

    return opt_dict, param_space


def launch_training(param_file, outdir):

    # check inputs
    assert os.path.isfile(param_file), f"Input file not found at {param_file}"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    opt_dict, param_space = parse_input(param_file)

    # single run
    if opt_dict is None:
        trainer(param_space, outdir)

    # optimizer run
    else:
        raise ValueError("Optimizer not finished")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="PINN-Schwarz", description="Schwarz-based training of PINN-PINN coupling")

    parser.add_argument("parameter_file")
    parser.add_argument("outdir")
    args = parser.parse_args()

    launch_training(args.parameter_file, args.outdir)
