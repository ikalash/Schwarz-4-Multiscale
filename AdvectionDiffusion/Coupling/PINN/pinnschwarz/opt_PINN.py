from ray import tune
import driver_schwarz

search_space = {
    # enter param space
    "hl": tune.grid_search(range(2, 11)),
    "nl": tune.grid_search(range(10, 85, 10))
}

def objective(config):
    cpu_time, n_iter = driver_schwarz.Driver(
        "/home/users/jinnyc/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/cases/example/input.csv",
        "/scratch/users/jinnyc/test/optimize",
        "/home/users/jinnyc/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/hyper.yaml",
        make_fig=False,
        percent_overlap=0.10,
        n_subdomains=3,
        hl=config['hl'],
        nl=config['nl']).train()
    return {"cpu": cpu_time, "iterations": n_iter}

tuner = tune.Tuner(
    tune.with_resources(objective, {"cpu": 10, "memory": 160000000000}),
    param_space=search_space
    )
results = tuner.fit()

print(results.get_best_result(metric="cpu", mode="min").config)

