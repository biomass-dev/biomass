import os
import warnings
warnings.filterwarnings('ignore')

if not os.path.isdir('./figure'):
    os.mkdir('./figure')

from  biomass.param_estim.dynamics import simulate_all

if __name__ == '__main__':
    simulate_all(viz_type='average',show_all=False,stdev=True)