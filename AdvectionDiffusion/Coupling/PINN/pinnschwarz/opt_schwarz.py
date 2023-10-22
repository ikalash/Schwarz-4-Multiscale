from ray import tune
import driver_schwartz
from ray.tune.search.bayesopt import BayesOptSearch

search_space = {
    # enter param space
    "n_subdomains": tune.grid_search([2, 3, 4, 5]),
    "percentage_overlap": tune.grid_search([0.1, 0.2, 0.24, 0.25, 0.3, 0.4])
}

def objective(config):
    cpu_time = driver_schwartz.Driver("/Users/maguo/Desktop/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/cases/example/input.csv",
                             "/Users/maguo/Desktop/test",
                             "/Users/maguo/Desktop/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/hyper.yaml",
                             make_fig=False,
                             n_subdomains=config['n_subdomains'],
                             percent_overlap=config['percentage_overlap']).train()
    return {"cpu": cpu_time} # TODO: need to look at loss


tuner = tune.Tuner(objective, param_space=search_space)
results = tuner.fit()

print(results.get_best_result(metric="cpu", mode="min").config)


#Uncomment the below lines to get bayes method
"""
search = {
    "percentage_overlap": (0.1, 0.5)
}
bayesopt = BayesOptSearch(search, metric="cpu", mode="min")
newtunner = tune.Tuner(
    objective, 
    tune_config=tune.TuneConfig(
        search_alg=bayesopt,
    )
)
    
result = newtunner.fit()
print(result.get_best_result(metric="cpu", mode="min").config)
"""
