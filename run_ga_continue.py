import os
import sys
import re
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

from biomass import optimize,optimize_continue

def run_ga_continue(nth_paramset):
    if not os.path.isdir('./out/%d'%(nth_paramset)):
        os.mkdir('./out/%d'%(nth_paramset))
        optimize(nth_paramset)
    else:
        optimize_continue(nth_paramset)
        
if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        run_ga_continue(int(args[1]))
    elif len(args) == 3:
        processes = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes)
        p.map(run_ga_continue, range(int(args[1]),int(args[2])+1))
        p.close()