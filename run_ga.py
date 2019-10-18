import os
import sys
import shutil
import re
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
    if args[1].isdigit():
        run_ga(int(args[1]))