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
            self.single_observable_options[i].divided_by = 60  # sec. -> min.
            self.single_observable_options[i].xlim = (-5, 95)
            self.single_observable_options[i].xticks = [0, 30, 60, 90]
            self.single_observable_options[i].xlabel = "Time (min)"
            self.single_observable_options[i].ylim = (-0.1, 1.3)
            self.single_observable_options[i].yticks = [0.0, 0.3, 0.6, 0.9, 1.2]
            self.single_observable_options[i].cmap = ["mediumblue", "red"]
            self.single_observable_options[i].shape = ["D", "s"]
            self.single_observable_options[i].dont_show = []

        self.single_observable_options[
            self.obs_names.index("Phosphorylated_MEKc")
        ].ylabel = "Phosphorylated MEK\n(cytoplasm)"

        self.single_observable_options[
            self.obs_names.index("Phosphorylated_ERKc")
        ].ylabel = "Phosphorylated ERK\n(cytoplasm)"

        self.single_observable_options[
            self.obs_names.index("Phosphorylated_RSKw")
        ].ylabel = "Phosphorylated RSK\n(whole cell)"

        self.single_observable_options[
            self.obs_names.index("Phosphorylated_CREBw")
        ].ylabel = "Phosphorylated CREB\n(whole cell)"

        self.single_observable_options[self.obs_names.index("dusp_mRNA")].ylabel = (
            r"$\it{dusp}$" + " mRNA\nexpression"
        )

        self.single_observable_options[self.obs_names.index("cfos_mRNA")].ylabel = (
            r"$\it{c}$" + "-" + r"$\it{fos}$" + " mRNA\nexpression"
        )

        self.single_observable_options[
            self.obs_names.index("cFos_Protein")
        ].ylabel = "c-Fos Protein\nexpression"

        self.single_observable_options[
            self.obs_names.index("Phosphorylated_cFos")
        ].ylabel = "Phosphorylated c-Fos\nProtein expression"

        return self.single_observable_options

    def get_multiple_observables_options(self):
        """
        self.multiple_observables_options.observables = [
            'Phosphorylated_ERKc',
            'Phosphorylated_CREBw',
            'Phosphorylated_cFos',
            'dusp_mRNA'
        ]
        self.multiple_observables_options.condition = 'EGF'
        self.multiple_observables_options.xlim = (-5, 95)
        self.multiple_observables_options.xticks = [0, 30, 60, 90]
        self.multiple_observables_options.xlabel = 'Time (min)'
        self.multiple_observables_options.ylim = (-0.1, 1.3)
        self.multiple_observables_options.yticks = [0.0, 0.3, 0.6, 0.9, 1.2]
        self.multiple_observables_options.ylabel = 'Intensity (a.u.)'
        """
        return self.multiple_observables_options

    def get_sensitivity_options(self):

        return self.sensitivity_options

    @staticmethod
    def set_timecourse_rcParams():
        """figure/simulation"""
        plt.rcParams["font.size"] = 20
        plt.rcParams["axes.linewidth"] = 1.5
        plt.rcParams["xtick.major.width"] = 1.5
        plt.rcParams["ytick.major.width"] = 1.5
        plt.rcParams["lines.linewidth"] = 1.8
        plt.rcParams["lines.markersize"] = 12
        plt.rcParams["savefig.bbox"] = "tight"
        # plt.rcParams["savefig.format"] = "pdf"
        # plt.rcParams["savefig.transparent"] = False
        # plt.rcParams["savefig.dpi"] = 1200
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
