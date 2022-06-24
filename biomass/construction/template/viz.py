from typing import List

from biomass.plotting import *
from matplotlib import pyplot as plt

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

    def __init__(self) -> None:
        super().__init__()

        self.cm = plt.cm.get_cmap("tab10")
        self.single_observable_options = [
            SingleObservable(self.cm, obs_name) for obs_name in self.obs_names
        ]
        self.multiple_observables_options = MultipleObservables(self.cm)
        self.sensitivity_options = SensitivityOptions(self.cm)

    def get_single_observable_options(self) -> List[SingleObservable]:

        return self.single_observable_options

    def get_multiple_observables_options(self) -> MultipleObservables:

        return self.multiple_observables_options

    def get_sensitivity_options(self) -> SensitivityOptions:

        return self.sensitivity_options

    @staticmethod
    def set_timecourse_rcParams() -> None:
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
    def set_sensitivity_rcParams() -> None:
        """figure/sensitivity"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        plt.rcParams["savefig.bbox"] = "tight"
        # plt.rcParams["savefig.format"] = "pdf"
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name: str) -> str:
        """figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions.
        """
        return name
