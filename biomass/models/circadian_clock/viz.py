from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from .observable import observables


class Visualization(object):
    """ Plotting parameters for customizing figure properties

    Attributes
    ----------
    timecourse_options : list of dict
        Plotting options for figure/simulation/<viz_type>/<each_observable>.
    
    multiplot_options : dict
        Plotting options for figure/simulation/<viz_type>/multiplot_observables.

    sensitivity_options : dict
        Plotting options for figure/sensitivity.

    """
    def __init__(self):
        self.cm = plt.cm.get_cmap('tab20')

        self.timecourse_options = [
            {
                'divided_by' : 1,  # to convert time unit. (e.g. sec -> min)
                'xlim' : (),
                'xticks' : None,
                'xlabel': 'Time',
                'ylim' : (),
                'yticks' : None,
                'ylabel': observables[i].replace('__', '\n').replace('_', ' '),
                'exp_data' : True,  # if False, experimental data will not be shown
                'cmap' : [self.cm.colors[j] for j in range(20)],
                'shape' : Line2D.filled_markers,
                'dont_show' : [],  # conditions you don't want to plot
            } for i, _ in enumerate(observables)]
        
        self.multiplot_options = {
            'fig_name' : 'multiplot_observables',
            'observables' : [],
            'condition' : None,
            'xlim' : (),
            'xticks' : None,
            'xlabel': 'Time',
            'ylim' : (),
            'yticks' : None,
            'ylabel': '',
            'cmap' : [self.cm.colors[j] for j in range(20)],
            'shape' : Line2D.filled_markers,
        }
    
        self.sensitivity_options = {
            'figsize' : (12, 5),
            'width' : 0.3,
            'legend_loc' : 'upper left',
            'cmap' : [self.cm.colors[j] for j in range(20)],
        }

    def get_timecourse_options(self):
        for i, _ in enumerate(observables):
            self.timecourse_options[i]['xlim'] = (0, 72)
            self.timecourse_options[i]['xticks'] = [0, 12, 24, 36, 48, 60, 72]
            self.timecourse_options[i]['xlabel'] = 'Time (h)'
        
        return self.timecourse_options

    def multiplot_observables(self):
        self.multiplot_options['fig_name'] = 'circadian_oscillations_in_continuous_darkness'
        self.multiplot_options['condition'] = 'DD'
        self.multiplot_options['observables'] = [
            'Per_mRNA',
            'Cry_mRNA',
            'Bmal1_mRNA',
        ]
        self.multiplot_options['xlim'] = (0, 72)
        self.multiplot_options['xticks'] = [0, 12, 24, 36, 48, 60, 72]
        self.multiplot_options['xlabel'] = 'Time (h)'
        self.multiplot_options['ylim'] = (0, 10)
        
        return self.multiplot_options
    
    @staticmethod
    def set_timecourse_rcParams():
        """ figure/simulation
        """
        plt.rcParams['font.size'] = 18
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.major.width'] = 1.5
        plt.rcParams['ytick.major.width'] = 1.5
        plt.rcParams['lines.linewidth'] = 1.8
        plt.rcParams['lines.markersize'] = 12
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_param_range_rcParams():
        """ figure/param_range
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def set_sensitivity_rcParams():
        """ figure/sensitivity
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """ figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        '''
        if name == 'MP':
            return 'Per mRNA'
        elif name == 'MC':
            return 'Cry mRNA'
        elif name == 'MB':
            return 'Bmal1 mRNA'
        '''
        return name