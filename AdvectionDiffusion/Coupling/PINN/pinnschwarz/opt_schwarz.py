from ray import tune
import driver_schwarz
from ray.tune.search.bayesopt import BayesOptSearch
from ray.tune.search.hyperopt import HyperOptSearch
from hyperopt import hp
from ray.train import RunConfig


class Optimizer():
    
    def __init__(self, opt_file):
        self.opt_file = opt_file
        self.keys = self.opt_file.keys()

    def objective(config):
        cpu_time = driver_schwarz.Driver("/Users/maguo/Desktop/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/cases/example/input.csv",
                                 "/Users/maguo/Desktop/test_optimizer",
                                 "/Users/maguo/Desktop/Schwarz-4-Multiscale/AdvectionDiffusion/Coupling/PINN/pinnschwarz/hyper.yaml",
                                 make_fig=False,
                                 n_subdomains=config['n_subdomains'],
                                 percent_overlap=config['percentage_overlap'],
                                 nl=config['nl'],
                                 hl=config['hl']).train()
        return {"cpu": cpu_time} # TODO: need to look at loss
    
    def grid_search(self, keys, objective):
        search_space = {}
        for key in keys:
            search_space[key] = tune.grid_search(list(self.opt_file[key]))
    
        tuner = tune.Tuner(objective, param_space=search_space)
        results = tuner.fit()
        
        print(results.get_best_result(metric="cpu", mode="min").config)
    
    
    #Uncomment the below lines to get bayes method
    
    def bayes():
    
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
    
    def hyperopt():
        #Try hyperopt
        search = {
            "percentage_overlap": hp.uniform('percentage_overlap', 0.1,0.4),
            "n_subdomains": hp.quniform("n_subdomains", 2, 6, 1),
            "nl": hp.quniform("nl", 20, 100, 1),
            "hl": hp.quniform("hl", 2, 20, 1)
        }
        
        current_best = [{
            "percentage_overlap": 0.2,
            "n_subdomains": 4,
            "nl": 23,
            "hl": 12,
            }]
        
        
        hot = HyperOptSearch(search, metric="cpu", mode="min", points_to_evaluate=current_best)
        newtunner = tune.Tuner(
            objective, 
            tune_config=tune.TuneConfig(
                search_alg=hot,
            ),
            run_config=RunConfig(name = "loggers"))
            
        result = newtunner.fit()
        print(result.get_best_result(metric="cpu", mode="min").config)
    
