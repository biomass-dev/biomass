import os
import sys
import shutil
import re
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

from biomass import optimize


def run_ga(nth_paramset):
    if not os.path.isdir('./out'):
        os.mkdir('./out')
    try:
        files = os.listdir('./out/%d'%(nth_paramset))
        for file in files:
            if any(map(file.__contains__,('.npy','log'))):
                os.remove('./out/%d/%s'%(nth_paramset,file))
    except:
        os.mkdir('./out/%d'%(nth_paramset))

    optimize(nth_paramset)
    
if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        run_ga(int(args[1]))
    elif len(args) == 3:
        processes = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes)
        p.map(run_ga, range(int(args[1]),int(args[2])+1))
        p.close()