import numpy as np

from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .ode import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of model observables.

    t : range
        Simulation time span.

    conditions : list of strings
        Experimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If :obj:`None`, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in ``sim.conditions`` will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names: list = [
            "Ski",
            "Skil",
            "Dnmt3a",
            "Sox4",
            "Jun",
            "Smad7",
            "Klf10",
            "Bmp4",
            "Cxcl15",
            "Dusp5",
            "Tgfa",
            "Pdk4",
        ]
        self.t: range = range(600 + 1)  # min
        self.conditions: list = ["WT", "Smad2OE", "Smad3OE", "Smad4OE"]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    @staticmethod
    def _set_gene_param(gene_name, x):
        if gene_name == "Ski":
            x[C.gene_turn] = x[C.Ski_turn]
            x[C.gene_act1] = x[C.Ski_act1]
            x[C.gene_act2] = x[C.Ski_act2]
            x[C.gene_act3] = x[C.Ski_act3]
            x[C.gene_inh1] = x[C.Ski_inh1]
            x[C.gene_inh2] = x[C.Ski_inh2]
            x[C.gene_inh3] = x[C.Ski_inh3]
        elif gene_name == "Skil":
            x[C.gene_turn] = x[C.Skil_turn]
            x[C.gene_act1] = x[C.Skil_act1]
            x[C.gene_act2] = x[C.Skil_act2]
            x[C.gene_act3] = x[C.Skil_act3]
            x[C.gene_inh1] = x[C.Skil_inh1]
            x[C.gene_inh2] = x[C.Skil_inh2]
            x[C.gene_inh3] = x[C.Skil_inh3]
        elif gene_name == "Dnmt3a":
            x[C.gene_turn] = x[C.Dnmt3a_turn]
            x[C.gene_act1] = x[C.Dnmt3a_act1]
            x[C.gene_act2] = x[C.Dnmt3a_act2]
            x[C.gene_act3] = x[C.Dnmt3a_act3]
            x[C.gene_inh1] = x[C.Dnmt3a_inh1]
            x[C.gene_inh2] = x[C.Dnmt3a_inh2]
            x[C.gene_inh3] = x[C.Dnmt3a_inh3]
        elif gene_name == "Sox4":
            x[C.gene_turn] = x[C.Sox4_turn]
            x[C.gene_act1] = x[C.Sox4_act1]
            x[C.gene_act2] = x[C.Sox4_act2]
            x[C.gene_act3] = x[C.Sox4_act3]
            x[C.gene_inh1] = x[C.Sox4_inh1]
            x[C.gene_inh2] = x[C.Sox4_inh2]
            x[C.gene_inh3] = x[C.Sox4_inh3]
        elif gene_name == "Jun":
            x[C.gene_turn] = x[C.Jun_turn]
            x[C.gene_act1] = x[C.Jun_act1]
            x[C.gene_act2] = x[C.Jun_act2]
            x[C.gene_act3] = x[C.Jun_act3]
            x[C.gene_inh1] = x[C.Jun_inh1]
            x[C.gene_inh2] = x[C.Jun_inh2]
            x[C.gene_inh3] = x[C.Jun_inh3]
        elif gene_name == "Smad7":
            x[C.gene_turn] = x[C.Smad7_turn]
            x[C.gene_act1] = x[C.Smad7_act1]
            x[C.gene_act2] = x[C.Smad7_act2]
            x[C.gene_act3] = x[C.Smad7_act3]
            x[C.gene_inh1] = x[C.Smad7_inh1]
            x[C.gene_inh2] = x[C.Smad7_inh2]
            x[C.gene_inh3] = x[C.Smad7_inh3]
        elif gene_name == "Klf10":
            x[C.gene_turn] = x[C.Klf10_turn]
            x[C.gene_act1] = x[C.Klf10_act1]
            x[C.gene_act2] = x[C.Klf10_act2]
            x[C.gene_act3] = x[C.Klf10_act3]
            x[C.gene_inh1] = x[C.Klf10_inh1]
            x[C.gene_inh2] = x[C.Klf10_inh2]
            x[C.gene_inh3] = x[C.Klf10_inh3]
        elif gene_name == "Bmp4":
            x[C.gene_turn] = x[C.Bmp4_turn]
            x[C.gene_act1] = x[C.Bmp4_act1]
            x[C.gene_act2] = x[C.Bmp4_act2]
            x[C.gene_act3] = x[C.Bmp4_act3]
            x[C.gene_inh1] = x[C.Bmp4_inh1]
            x[C.gene_inh2] = x[C.Bmp4_inh2]
            x[C.gene_inh3] = x[C.Bmp4_inh3]
        elif gene_name == "Cxcl15":
            x[C.gene_turn] = x[C.Cxcl15_turn]
            x[C.gene_act1] = x[C.Cxcl15_act1]
            x[C.gene_act2] = x[C.Cxcl15_act2]
            x[C.gene_act3] = x[C.Cxcl15_act3]
            x[C.gene_inh1] = x[C.Cxcl15_inh1]
            x[C.gene_inh2] = x[C.Cxcl15_inh2]
            x[C.gene_inh3] = x[C.Cxcl15_inh3]
        elif gene_name == "Dusp5":
            x[C.gene_turn] = x[C.Dusp5_turn]
            x[C.gene_act1] = x[C.Dusp5_act1]
            x[C.gene_act2] = x[C.Dusp5_act2]
            x[C.gene_act3] = x[C.Dusp5_act3]
            x[C.gene_inh1] = x[C.Dusp5_inh1]
            x[C.gene_inh2] = x[C.Dusp5_inh2]
            x[C.gene_inh3] = x[C.Dusp5_inh3]
        elif gene_name == "Tgfa":
            x[C.gene_turn] = x[C.Tgfa_turn]
            x[C.gene_act1] = x[C.Tgfa_act1]
            x[C.gene_act2] = x[C.Tgfa_act2]
            x[C.gene_act3] = x[C.Tgfa_act3]
            x[C.gene_inh1] = x[C.Tgfa_inh1]
            x[C.gene_inh2] = x[C.Tgfa_inh2]
            x[C.gene_inh3] = x[C.Tgfa_inh3]
        elif gene_name == "Pdk4":
            x[C.gene_turn] = x[C.Pdk4_turn]
            x[C.gene_act1] = x[C.Pdk4_act1]
            x[C.gene_act2] = x[C.Pdk4_act2]
            x[C.gene_act3] = x[C.Pdk4_act3]
            x[C.gene_inh1] = x[C.Pdk4_inh1]
            x[C.gene_inh2] = x[C.Pdk4_inh2]
            x[C.gene_inh3] = x[C.Pdk4_inh3]
        return x

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation

        for gene_name in self.obs_names:
            x = self._set_gene_param(gene_name, x)
            for i, condition in enumerate(self.conditions):
                if condition == "WT":
                    y0[V.S2] = x[C.S2tot]
                    y0[V.S3] = x[C.S3tot]
                    y0[V.S4] = x[C.S4tot]
                elif condition == "Smad2OE":
                    y0[V.S2] = 2 * x[C.S2tot]
                    y0[V.S3] = x[C.S3tot]
                    y0[V.S4] = x[C.S4tot]
                elif condition == "Smad3OE":
                    y0[V.S2] = x[C.S2tot]
                    y0[V.S3] = 16 * x[C.S3tot]
                    y0[V.S4] = x[C.S4tot]
                elif condition == "Smad4OE":
                    y0[V.S2] = x[C.S2tot]
                    y0[V.S3] = x[C.S3tot]
                    y0[V.S4] = 3 * x[C.S3tot]

                sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

                if sol is None:
                    return False
                else:
                    self.simulations[self.obs_names.index(gene_name), i] = np.log2(sol.y[V.gene])
        return None

    def set_data(self):
        pass

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return []
        assert False
