import os
import sys
import re
import warnings
warnings.filterwarnings('ignore')

if not os.path.isdir('./figure'):
    os.mkdir('./figure')

from biomass.param_estim.dynamics import simulate_all

n_file = []
if os.path.isdir('./out'):
    fit_param_files = os.listdir('./out')
    for file in fit_param_files:
        if re.match(r'\d', file):
            n_file.append(int(file))
                
if __name__ == '__main__':
    args = sys.argv
    if len(args) == 1 or 4 < len(args):
        raise TypeError(
            '\n$ python run_sim.py [viz_type] [show_all] [stdev]\n'
        )
    else:
        if str(args[1]) not in ['best', 'average', 'original', 'experiment']:
            try:
                int(args[1])
            except ValueError:
                print(
                    "Avairable viz_type are: 'best','average','original','experiment',n(=1,2,...)"
                )
            if len(n_file) < int(args[1]):
                raise ValueError(
                    'n ({:d}) must be smaller than n_fit_param ({:d})'.format(
                        int(args[1]), len(n_file)
                    )
                )
        if len(args) == 2:
            show_all, stdev = (False, False)
        elif len(args) == 3:
            if str(args[2]) == 'show_all':
                show_all, stdev = (True, False)
            elif str(args[2]) == 'stdev':
                show_all, stdev = (False, True)
            else:
                raise ValueError(
                    '[show_all] or [stdev]\n'
                )
        elif len(args) == 4:
            if ((str(args[2]) == 'show_all' and str(args[3]) == 'stdev')
                or (str(args[2]) == 'stdev' and str(args[3]) == 'show_all')):
                show_all, stdev = (True, True)
            else:
                raise ValueError(
                    '[show_all] or [stdev]\n'
                )

    simulate_all(viz_type=str(args[1]), show_all=show_all, stdev=stdev)