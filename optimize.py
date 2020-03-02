import os
import sys
import shutil
import re
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

from biomass.param_estim import optimize


def run_ga(nth_paramset):
    if not os.path.isdir('./out'):
        os.mkdir('./out')
    try:
        files = os.listdir(
            './out/{:d}'.format(
                nth_paramset
            )
        )
        for file in files:
            if any(map(file.__contains__, ('.npy', '.log'))):
                os.remove(
                    './out/{:d}/{}'.format(
                        nth_paramset, file
                    )
                )
    except FileNotFoundError:
        os.mkdir(
            './out/{:d}'.format(
                nth_paramset
            )
        )

    optimize(nth_paramset)

if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        run_ga(int(args[1]))
    elif len(args) == 3:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(run_ga, range(int(args[1]), int(args[2]) + 1))
        p.close()