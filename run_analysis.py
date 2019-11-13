import os
import warnings
warnings.filterwarnings('ignore')
from biomass.analysis import reaction

if not os.path.isdir('./figure'):
    os.mkdir('./figure')
    
if __name__ == '__main__':
    reaction.sensitivity_barplot()
    reaction.sensitivity_heatmap()