import os
import sys
import re
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
    if args[1].isdigit():
        run_ga_continue(int(args[1]))