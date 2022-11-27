import matplotlib as mpl
import matplotlib.pyplot as plt

from biomass.plotting import *

from .observable import Observable


class Visualization(Observable):
    """
    Plotting parameters for customizing figure properties.

    Attributes
    ----------
    cm : matplotlib.colors.ListedColormap (default: ``mpl.colormaps['tab10']``)
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
            self.single_observable_options[i].xlim = (0, 72)
            self.single_observable_options[i].xticks = [0, 12, 24, 36, 48, 60, 72]
            self.single_observable_options[i].xlabel = "Time (h)"

        return self.single_observable_options

    def get_multiple_observables_options(self):
        self.multiple_observables_options.fname = "circadian_oscillations_in_continuous_darkness"
        self.multiple_observables_options.condition = "DD"
        self.multiple_observables_options.observables = [
            "Per_mRNA",
            "Cry_mRNA",
            "Bmal1_mRNA",
        ]
        self.multiple_observables_options.xlim = (0, 72)
        self.multiple_observables_options.xticks = [0, 12, 24, 36, 48, 60, 72]
        self.multiple_observables_options.xlabel = "Time (h)"
        self.multiple_observables_options.ylim = (0, 10)

        return self.multiple_observables_options

    def get_sensitivity_options(self):

        return self.sensitivity_options

    @staticmethod
    def set_timecourse_rcParams():
        """figure/simulation"""
        plt.rcParams["font.size"] = 18
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
        if name == 'MP':
            return 'Per mRNA'
        elif name == 'MC':
            return 'Cry mRNA'
        elif name == 'MB':
            return 'Bmal1 mRNA'
        """
        return name
