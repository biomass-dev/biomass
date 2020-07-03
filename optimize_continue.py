import os
import sys
import re
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

from biomass.models.Nakakuki_Cell_2010 import SearchParam, objective
from biomass.models.Nakakuki_Cell_2010 import __path__ as MODEL_PATH
from biomass.ga import GeneticAlgorithmContinue

if __name__ == '__main__':
    ga_continue = GeneticAlgorithmContinue(
        model_path=MODEL_PATH[0],
        sp=SearchParam(),
        obj_func=objective
    )
    args = sys.argv
    if len(args) == 2:
        ga_continue.run(int(args[1]))
    elif len(args) == 3:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_continue.run, range(int(args[1]), int(args[2]) + 1))
        p.close()