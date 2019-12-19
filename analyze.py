import os
import warnings
warnings.filterwarnings('ignore')
from biomass.analysis import reaction, nonzero_init

if not os.path.isdir('./figure'):
    os.mkdir('./figure')
    
if __name__ == '__main__':
    reaction.sensitivity_barplot(metric='integral')
#   reaction.sensitivity_heatmap(metric='integral')
#   nonzero_init.sensitivity_barplot(metric='integral')
#   nonzero_init.sensitivity_heatmap(metric='integral')