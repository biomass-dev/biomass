import matplotlib as mpl
import matplotlib.pyplot as plt

from biomass.plotting import *

from .observable import Observable


class Visualization(Observable):
    """
    Plotting parameters for customizing figure properties.

    Attributes
    ----------
    cm : matplotlib.colors.ListedColormap (default: ``plt.cm.get_cmap('tab10')``)
        Choosing colormaps for ``cmap``.
    single_observable_options : list of SingleObservable
        Visualization options for time-course simulation (single-observable).
    multiple_observables_options : MultipleObservables
        Visualization options for time-course simulation (multi-observables).
    sensitivity_options : SensitivityOptions
        Visualization options for sensitivity analysis results.
    """

    def __init__(self):
        super().__init__()

        self.cm = mpl.colormaps["tab10"]
        self.single_observable_options = [
            SingleObservable(self.cm, obs_name) for obs_name in self.obs_names
        ]
        self.multiple_observables_options = MultipleObservables(self.cm)
        self.sensitivity_options = SensitivityOptions(self.cm)

    def get_single_observable_options(self):
        for i, _ in enumerate(self.obs_names):
            self.single_observable_options[i].divided_by = 60
            self.single_observable_options[i].xlim = (0, 150)
            self.single_observable_options[i].xticks = [20 * i for i in range(8)]
            self.single_observable_options[i].xlabel = "TIME (min)"
            self.single_observable_options[i].ylim = (0, 310)
            self.single_observable_options[i].yticks = [50 * i for i in range(7)]
            self.single_observable_options[i].exp_data = False
            self.single_observable_options[
                self.obs_names.index("biphosphorylated_MAPK")
            ].ylabel = "ERK-PP"
            self.single_observable_options[
                self.obs_names.index("unphosphorylated_MAPK")
            ].ylabel = "ERK"

        return self.single_observable_options

    def get_multiple_observables_options(self):
        self.multiple_observables_options.fname = "oscillations_in_the_MAPK_concentrations"
        self.multiple_observables_options.observables = [
            "biphosphorylated_MAPK",
            "unphosphorylated_MAPK",
        ]
        self.multiple_observables_options.condition = "control"
        self.multiple_observables_options.xlim = (0, 150)
        self.multiple_observables_options.xticks = [20 * i for i in range(8)]
        self.multiple_observables_options.xlabel = "TIME (min)"
        self.multiple_observables_options.ylim = (0, 310)
        self.multiple_observables_options.yticks = [50 * i for i in range(7)]
        self.multiple_observables_options.ylabel = "MAPK concentrations (nM)"

        return self.multiple_observables_options

    def get_sensitivity_options(self):

        return self.sensitivity_options

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
        # plt.rcParams["savefig.format"] = "pdf"
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
