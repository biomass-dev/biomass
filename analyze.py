import os
import warnings
import sys
warnings.filterwarnings('ignore')
from biomass.analysis import reaction, nonzero_init

if not os.path.isdir('./figure'):
    os.mkdir('./figure')

if __name__ == '__main__':
    args = sys.argv
    if len(args) > 3:
        raise TypeError(
            "Too many arguments!"
        )
    else:
        if str(args[1]) not in ['integral', 'amplitude', 'duration']:
            raise ValueError(
                "Available arguments are: 'integral', 'amplitude', 'duration'"
            )
        if str(args[2]) not in ['barplot', 'heatmap']:
            raise ValueError(
                "Available arguments are: 'barplot', 'heatmap'"
            )
        reaction.analyze(metric=str(args[1]), style=str(args[2]))
#       nonzero_init.analyze(metric=str(args[1]), style=str(args[2]))
