import os
import warnings
import sys
warnings.filterwarnings('ignore')
from biomass.analysis import reaction, nonzero_init

if not os.path.isdir('./figure'):
    os.mkdir('./figure')

if __name__ == '__main__':
    args = sys.argv
    if len(args) > 2:
        print(
                "Too many arguments! Try one of these: 'integral', 'amplitude', 'duration'"
            )
        sys.exit()
    else:
        if str(args[1]) not in ['integral', 'amplitude', 'duration']:
            raise ValueError(
                "Available arguments are: 'integral', 'amplitude', 'duration'"
            )
        else:
            reaction.sensitivity_barplot(metric=str(args[1]))
#           reaction.sensitivity_heatmap(metric=str(args[1]))
#           nonzero_init.sensitivity_barplot(metric=str(args[1]))
#           nonzero_init.sensitivity_heatmap(metric=str(args[1])) 