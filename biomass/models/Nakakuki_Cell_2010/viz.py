from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from .observable import observables


class Visualization(object):
    """
    Plotting parameters for customizing figure properties

    Attributes
    ----------
    cm : matplotlib.colors.ListedColormap (default: plt.cm.get_cmap('tab10'))
        Choosing colormaps for 'cmap'.

    timecourse_options : list of dict
        Plotting options for figure/simulation/<viz_type>/<each_observable>.

            Keys
            ----
            * 'divided_by' : int or float (default: 1)
                Convert time unit. (e.g. sec -> min).

            * 'xlim' : tuple
                Set the x limits of the current axes.

            * 'xticks' : list (default: None)
                Set the current tick locations of the x-axis.

            * 'xlabel' : str (default: 'Time')
                Set the label for the x-axis.

            * 'ylim' : tuple
                Set the y limits of the current axes.

            * 'yticks' : list (default: None)
                Set the current tick locations of the y-axis.

            * 'ylabel' : str (default: observables[i].replace('__', '\n').replace('_', ' '))
                Set the label for the y-axis.

            * 'exp_data' : bool (default: True)
                if False, experimental data will not be shown.

            * 'legend_loc' : Location String (default: None)
                Set the location of the legend. If 'legend_loc' is None,
                pyplot.legend will not be shown.

            * 'cmap' : list or tuple
                Set colormap.

            * 'shape' : list or tuple of strings (default: Line2D.filled_markers)
                Set markers.

            * 'dont_show' : list of strings
                Set conditions you don't want to plot.

    multiplot_options : dict
        Plotting options for figure/simulation/<viz_type>/multiplot_observables.

    sensitivity_options : dict
        Plotting options for figure/sensitivity.

    """

    def __init__(self):
        self.cm = plt.cm.get_cmap("tab10")

        self.timecourse_options = [
            {
                "divided_by": 1,
                "xlim": (),
                "xticks": None,
                "xlabel": "Time",
                "ylim": (),
                "yticks": None,
                "ylabel": observables[i].replace("__", "\n").replace("_", " "),
                "exp_data": True,
                "legend_loc": None,
                "cmap": [self.cm.colors[j] for j in range(10)],
                "shape": Line2D.filled_markers,
                "dont_show": [],
            }
            for i, _ in enumerate(observables)
        ]

        self.multiplot_options = {
            "fig_name": "multiplot_observables",
            "observables": [],
            "condition": None,
            "xlim": (),
            "xticks": None,
            "xlabel": "Time",
            "ylim": (),
            "yticks": None,
            "ylabel": "",
            "cmap": [self.cm.colors[j] for j in range(10)],
            "shape": Line2D.filled_markers,
        }

        self.sensitivity_options = {
            "figsize": (12, 5),
            "width": 0.3,
            "legend_loc": "upper left",
            "cmap": [self.cm.colors[j] for j in range(10)],
        }

    def get_timecourse_options(self):
        """
        for i, _ in enumerate(observables):
            self.timecourse_options[i]['divided_by'] = 60  # sec. -> min.
            self.timecourse_options[i]['xlim'] = (-5, 95)
            self.timecourse_options[i]['xticks'] = [0, 30, 60, 90]
            self.timecourse_options[i]['xlabel'] = 'Time (min)'
            self.timecourse_options[i]['ylim'] = (-0.1, 1.3)
            self.timecourse_options[i]['yticks'] = [0.0, 0.3, 0.6, 0.9, 1.2]
            self.timecourse_options[i]['cmap'] = ['mediumblue', 'red']
            self.timecourse_options[i]['shape'] = ['D', 's']
            self.timecourse_options[i]['dont_show'] = []

        self.timecourse_options[
            observables.index('Phosphorylated_MEKc')
        ]['ylabel'] = 'Phosphorylated MEK\n(cytoplasm)'

        self.timecourse_options[
            observables.index('Phosphorylated_ERKc')
        ]['ylabel'] = 'Phosphorylated ERK\n(cytoplasm)'

        self.timecourse_options[
            observables.index('Phosphorylated_RSKw')
        ]['ylabel'] = 'Phosphorylated RSK\n(whole cell)'

        self.timecourse_options[
            observables.index('Phosphorylated_CREBw')
        ]['ylabel'] = 'Phosphorylated CREB\n(whole cell)'

        self.timecourse_options[
            observables.index('dusp_mRNA')
        ]['ylabel'] = r'$\it{dusp}$'+' mRNA\nexpression'

        self.timecourse_options[
            observables.index('cfos_mRNA')
        ]['ylabel'] = r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA\nexpression'

        self.timecourse_options[
            observables.index('cFos_Protein')
        ]['ylabel'] = 'c-Fos Protein\nexpression'

        self.timecourse_options[
            observables.index('Phosphorylated_cFos')
        ]['ylabel'] = 'Phosphorylated c-Fos\nProtein expression'
        """
        return self.timecourse_options

    def multiplot_observables(self):
        """
        self.multiplot_options['observables'] = [
            'Phosphorylated_ERKc',
            'Phosphorylated_CREBw',
            'Phosphorylated_cFos',
            'dusp_mRNA'
        ]
        self.multiplot_options['condition'] = 'EGF'
        self.multiplot_options['xlim'] = (-5, 95)
        self.multiplot_options['xticks'] = [0, 30, 60, 90]
        self.multiplot_options['xlabel'] = 'Time (min)'
        self.multiplot_options['ylim'] = (-0.1, 1.3)
        self.multiplot_options['yticks'] = [0.0, 0.3, 0.6, 0.9, 1.2]
        self.multiplot_options['ylabel'] = 'Intensity (a.u.)'
        """
        return self.multiplot_options

    @staticmethod
    def set_timecourse_rcParams():
        """figure/simulation"""
        plt.rcParams["font.size"] = 20
        plt.rcParams["axes.linewidth"] = 1.5
        plt.rcParams["xtick.major.width"] = 1.5
        plt.rcParams["ytick.major.width"] = 1.5
        plt.rcParams["lines.linewidth"] = 1.8
        plt.rcParams["lines.markersize"] = 12
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_param_range_rcParams():
        """figure/param_range"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def set_sensitivity_rcParams():
        """figure/sensitivity"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        """
        if name == 'ERKc':
            return 'ERK (cytoplasm)'
        elif name == 'RSKc':
            return 'RSK (cytoplasm)'
        elif name == 'CREBn':
            return 'CREB (nucleus)'
        elif name == 'Elk1n':
            return 'Elk1 (nucleus)'
        """
        return name
