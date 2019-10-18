import os
import warnings
warnings.filterwarnings('ignore')
from biomass.analysis import visualize_sensivity

if not os.path.isdir('./figure'):
    os.mkdir('./figure')
    
if __name__ == '__main__':
    visualize_sensivity()