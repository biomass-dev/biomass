import os
import sys
import shutil
import re
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

from biomass.models.Nakakuki_Cell_2010 import (get_search_region,
                                                decode_gene2val, objective)
from biomass.ga import GeneticAlgorithmInit

if __name__ == '__main__':
    ga_init = GeneticAlgorithmInit(
        get_search_region, decode_gene2val, objective
    )
    args = sys.argv
    if len(args) == 2:
        ga_init.run(int(args[1]))
    elif len(args) == 3:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_init.run, range(int(args[1]), int(args[2]) + 1))
        p.close()