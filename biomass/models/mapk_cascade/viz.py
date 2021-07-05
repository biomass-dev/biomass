from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from .observable import Observable


class Visualization(Observable):
    """
    Plotting parameters for customizing figure properties

    Attributes
    ----------
    cm : matplotlib.colors.ListedColormap (default: plt.cm.get_cmap('tab10'))
        Choosing colormaps for 'cmap'.

    timecourse_options : list of dict
        Plotting options for figure/simulation/<viz_type>/<each_observable>.

        * 'divided_by' : int or float (default: 1)
            Convert time unit. (e.g. sec -> min).

        * 'figsize' : tuple (default: (4, 3))
                Width, height in inches.

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

        * 'ylabel' : str (default: self.obs_names[i].replace('__', '\n').replace('_', ' '))
            Set the label for the y-axis.

        * 'exp_data' : bool (default: True)
            if False, experimental data will not be shown.

        * 'legend_kws' : dict, optional
            Keyword arguments to pass to matplotlib.pyplot.legend().

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
        super().__init__()

        self.cm = plt.cm.get_cmap("tab10")

        self.timecourse_options = [
            {
                "divided_by": 1,
                "figsize": (4, 3),
                "xlim": (),
                "xticks": None,
                "xlabel": "Time",
                "ylim": (),
                "yticks": None,
                "ylabel": self.obs_names[i].replace("__", "\n").replace("_", " "),
                "exp_data": True,
                "legend_kws": None,
                "cmap": [self.cm.colors[j] for j in range(10)],
                "shape": Line2D.filled_markers,
                "dont_show": [],
            }
            for i, _ in enumerate(self.obs_names)
        ]

        self.multiplot_options = {
            "fname": "multiplot_observables",
            "figsize": (4, 3),
            "observables": [],
            "condition": None,
            "xlim": (),
            "xticks": None,
            "xlabel": "Time",
            "ylim": (),
            "yticks": None,
            "ylabel": "",
            "legend_kws": dict(
                bbox_to_anchor=(1.05, 1),
                loc="upper left",
                borderaxespad=0,
                labelspacing=1.25,
                frameon=False,
                fontsize=12,
            ),
            "cmap": [self.cm.colors[j] for j in range(10)],
            "shape": Line2D.filled_markers,
        }

        self.sensitivity_options = {
            "figsize": (12, 5),
            "width": 0.3,
            "legend_kws": dict(loc="upper left", frameon=False),
            "cmap": [self.cm.colors[j] for j in range(10)],
        }

    def get_timecourse_options(self):
        for i, _ in enumerate(self.obs_names):
            self.timecourse_options[i]["divided_by"] = 60
            self.timecourse_options[i]["xlim"] = (0, 150)
            self.timecourse_options[i]["xticks"] = [20 * i for i in range(8)]
            self.timecourse_options[i]["xlabel"] = "TIME (min)"
            self.timecourse_options[i]["ylim"] = (0, 310)
            self.timecourse_options[i]["yticks"] = [50 * i for i in range(7)]
            self.timecourse_options[i]["exp_data"] = False

            self.timecourse_options[self.obs_names.index("biphosphorylated_MAPK")][
                "ylabel"
            ] = "ERK-PP"
            self.timecourse_options[self.obs_names.index("unphosphorylated_MAPK")][
                "ylabel"
            ] = "ERK"

        return self.timecourse_options

    def multiplot_observables(self):
        self.multiplot_options["fname"] = "oscillations_in_the_MAPK_concentrations"
        self.multiplot_options["observables"] = [
            "biphosphorylated_MAPK",
            "unphosphorylated_MAPK",
        ]
        self.multiplot_options["condition"] = "control"
        self.multiplot_options["xlim"] = (0, 150)
        self.multiplot_options["xticks"] = [20 * i for i in range(8)]
        self.multiplot_options["xlabel"] = "TIME (min)"
        self.multiplot_options["ylim"] = (0, 310)
        self.multiplot_options["yticks"] = [50 * i for i in range(7)]
        self.multiplot_options["ylabel"] = "MAPK concentrations (nM)"

        return self.multiplot_options

    @staticmethod
    def set_timecourse_rcParams():
        """figure/simulation"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.5
        plt.rcParams["xtick.major.width"] = 1.5
        plt.rcParams["ytick.major.width"] = 1.5
        plt.rcParams["lines.linewidth"] = 1.8
        plt.rcParams["lines.markersize"] = 12
        plt.rcParams["savefig.bbox"] = "tight"
        # plt.rcParams["savefig.format"] = "pdf"
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_sensitivity_rcParams():
        """figure/sensitivity"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        plt.rcParams["savefig.bbox"] = "tight"
        plt.rcParams["savefig.format"] = "pdf"
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        """
        if name == 'MKKK':
            return 'Raf'
        elif name == 'MKK':
            return 'MEK'
        elif name == 'MAPK':
            return 'ERK'
        """
        return name
