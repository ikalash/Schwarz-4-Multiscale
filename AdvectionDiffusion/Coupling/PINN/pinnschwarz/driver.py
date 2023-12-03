import argparse
import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import yaml
import ray
from ray import tune
from ray.tune.search.basic_variant import BasicVariantGenerator
from ray.tune.search.hyperopt import HyperOptSearch
from ray.tune.search.bayesopt import BayesOptSearch

from pinnschwarz.trainer import trainer
from pinnschwarz.utils import get_resources


def parse_param(param_dict, param_str, can_opt=False, has_opt=False, default=None):
    if param_str in param_dict:
        param = param_dict[param_str]

        # define search space for parameter
        if isinstance(param, dict):
            assert can_opt and has_opt, f"Must set optimizer params if setting search space for parameter: {param_str}"

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


def append_param_space(param_space, param_dict, param_str, can_opt, has_opt, default=None):
    assert param_str not in param_space, f"Invalid double insertion of parameter {param_str} into parameter space"

    param_val = parse_param(param_dict, param_str, can_opt=can_opt, has_opt=has_opt, default=default)
    param_space[param_str] = param_val

    return param_space


def parse_input(param_file):
    with open(param_file, "r") as file:
        param_dict = yaml.safe_load(file)
    param_space = {}

    # optimization parameters
    if "optimizer" in param_dict.keys():
        opt = True
        opt_dict = param_dict["optimizer"]
        # some error checking, update as necessary
        try:
            assert opt_dict["algo"] in ["random", "grid", "hyperopt", "bayes"]
        except KeyError:
            opt_dict["algo"] = "random"
        except AssertionError:
            print("Valid optimizer 'algo's: random, grid, hyperopt, bayes")
            raise

        if "n_samples" not in opt_dict:
            opt_dict["n_samples"] = 100

    else:
        opt = False
        opt_dict = None

    # PDE parameters
    pde_dict = param_dict["pde"]
    param_space = append_param_space(param_space, pde_dict, "order", False, opt, default=2)
    param_space = append_param_space(param_space, pde_dict, "n_colloc", False, opt, default=1000)
    param_space = append_param_space(param_space, pde_dict, "domain_bounds", False, opt, default=[0.0, 1.0])
    param_space = append_param_space(param_space, pde_dict, "beta", False, opt, default=1.0)
    param_space = append_param_space(param_space, pde_dict, "peclet", opt, False)

    # Schwarz parameters
    schwarz_dict = param_dict["schwarz"]
    param_space = append_param_space(param_space, schwarz_dict, "n_subdomains", True, opt)
    param_space = append_param_space(param_space, schwarz_dict, "percent_overlap", True, opt)
    param_space = append_param_space(param_space, schwarz_dict, "sys_bcs", True, opt)
    param_space = append_param_space(param_space, schwarz_dict, "sch_bcs", True, opt)
    param_space = append_param_space(param_space, schwarz_dict, "BC_type", False, opt, "DD")
    param_space = append_param_space(param_space, schwarz_dict, "snap_loss", True, opt)
    param_space = append_param_space(param_space, schwarz_dict, "schwarz_tol", False, opt, 0.001)
    param_space = append_param_space(param_space, schwarz_dict, "err_tol", False, opt, 0.005)

    # neural network parameters
    nn_dict = param_dict["nn"]
    param_space = append_param_space(param_space, nn_dict, "alpha", True, opt, 0.2)
    param_space = append_param_space(param_space, nn_dict, "n_epochs", True, opt, 1000)
    param_space = append_param_space(param_space, nn_dict, "learn_rate", True, opt, 0.001)
    param_space = append_param_space(param_space, nn_dict, "n_layers", True, opt, 2)
    param_space = append_param_space(param_space, nn_dict, "n_neurons", True, opt, 20)

    return opt_dict, param_space


def launch_training(param_file, outdir):
    # check inputs
    assert os.path.isfile(param_file), f"Input file not found at {param_file}"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    opt_dict, param_space = parse_input(param_file)

    # single run
    if opt_dict is None:
        time, iters, err = trainer(param_space, outdir, make_figs=True)

    # optimizer run
    else:

        def objective(config):
            time, iters, err = trainer(config, outdir, make_figs=False)
            return {"time": time}

        # establish available CPU/memory resources
        avail_cpu, avail_mem = get_resources()
        mem_per_cpu = avail_mem / avail_cpu
        ray.init(
            include_dashboard=False,
            num_cpus=avail_cpu,
            object_store_memory=0.3 * avail_mem,
        )

        # set up sampling algorithm
        algo = opt_dict["algo"]
        if algo in ["random", "grid"]:
            if algo == "grid":
                grid_flag = True
            else:
                grid_flag = False
            algo_obj = BasicVariantGenerator(constant_grid_search=grid_flag, random_state=0)

        elif algo == "hyperopt":
            if "n_initial_points" in opt_dict:
                n_initial_points = int(opt_dict["n_initial_points"])
            else:
                n_initial_points = 20
            algo_obj = HyperOptSearch(n_initial_points=n_initial_points, random_state_seed=0)

        elif algo == "bayes":
            for _, param in param_space.items():
                if issubclass(type(param), ray.tune.search.sample.Domain):
                    assert isinstance(param.sampler, ray.tune.search.sample.Float._Uniform), "Bayes only permits `uniform` search spaces"
            algo_obj = BayesOptSearch(random_state=0)

        # fit and report best parameter
        tuner = tune.Tuner(
            tune.with_resources(objective, {"cpu": 1, "memory": mem_per_cpu}),
            tune_config=tune.TuneConfig(
                search_alg=algo_obj,
                num_samples=opt_dict["n_samples"],
                metric="time",
                mode="min",
                reuse_actors=False,
            ),
            param_space=param_space,
        )
        results = tuner.fit()
        print(results.get_best_result().config)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="PINN-Schwarz", description="Schwarz-based training of PINN-PINN coupling")

    parser.add_argument("parameter_file")
    parser.add_argument("outdir")
    args = parser.parse_args()

    launch_training(args.parameter_file, args.outdir)
