from ray import tune
import driver

search_space = {
    # enter param space
    "hl": tune.grid_search([2, 3]),
    "nl": tune.grid_search([10, 20])
}

def objective(config):
    cpu_time = driver.Driver("/Users/jinnychung/Documents/Stanford/Sandia/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/cases/example/input.csv",
                             "/Users/jinnychung/Documents/Stanford/Sandia/test",
                             "/Users/jinnychung/Documents/Stanford/Sandia/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/hyper.yaml",
                             make_fig=False,
                             hl=config['hl'],
                             nl=config['nl']).train()
    return {"cpu": cpu_time} # TODO: need to look at loss

tuner = tune.Tuner(objective, param_space=search_space)
results = tuner.fit()

print(results.get_best_result(metric="cpu", mode="min").config)

