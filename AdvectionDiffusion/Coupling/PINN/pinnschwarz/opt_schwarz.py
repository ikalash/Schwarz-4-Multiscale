from ray import tune
import driver_schwarz
from ray.tune.search.bayesopt import BayesOptSearch
from ray.tune.search.hyperopt import HyperOptSearch
from hyperopt import hp


def objective(config):
    cpu_time = driver_schwarz.Driver(
        "/Users/maguo/Desktop/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/cases/example/input.csv",
        "/Users/maguo/Desktop/test_optimizer",
        "/Users/maguo/Desktop/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/hyper.yaml",
        make_fig=False,
        n_subdomains=config["n_subdomains"],
        percent_overlap=config["percentage_overlap"],
        nl=config["nl"],
        hl=config["hl"],
    ).train()
    return {"cpu": cpu_time}  # TODO: need to look at loss


def grid_search():
    search_space = {
        # enter param space
        "n_subdomains": tune.grid_search([2, 3, 4, 5]),
        "percentage_overlap": tune.grid_search([0.1, 0.2, 0.24, 0.25, 0.3, 0.4]),
    }

    tuner = tune.Tuner(objective, param_space=search_space)
    results = tuner.fit()

    print(results.get_best_result(metric="cpu", mode="min").config)


# Uncomment the below lines to get bayes method


def bayes():
    search = {"percentage_overlap": (0.1, 0.5)}
    bayesopt = BayesOptSearch(search, metric="cpu", mode="min")
    newtunner = tune.Tuner(
        objective,
        tune_config=tune.TuneConfig(
            search_alg=bayesopt,
        ),
    )

    result = newtunner.fit()
    print(result.get_best_result(metric="cpu", mode="min").config)


# Try hyperopt
search = {
    "percentage_overlap": hp.uniform("percentage_overlap", 0.1, 0.4),
    "n_subdomains": hp.quniform("n_subdomains", 2, 6, 1),
    "nl": hp.quniform("nl", 20, 100, 1),
    "hl": hp.quniform("hl", 2, 20, 1),
}

current_best = [
    {
        "percentage_overlap": 0.2,
        "n_subdomains": 2,
        "nl": 20,
        "hl": 2,
    }
]


hot = HyperOptSearch(search, metric="cpu", mode="min", points_to_evaluate=current_best)
newtunner = tune.Tuner(
    objective,
    tune_config=tune.TuneConfig(
        search_alg=hot,
    ),
)

result = newtunner.fit()
print(result.get_best_result(metric="cpu", mode="min").config)
