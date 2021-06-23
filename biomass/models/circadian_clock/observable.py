import numpy as np

from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .set_model import DifferentialEquation


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
        Expetimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If None, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in sim.conditions will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names: list = [
            "Per_mRNA",
            "Cry_mRNA",
            "Bmal1_mRNA",
        ]
        self.t: range = range(72 + 1)
        self.conditions: list = ["DD"]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.t), len(self.conditions))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        for i, condition in enumerate(self.conditions):
            if condition == "DD":
                pass

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("Per_mRNA"), :, i] = sol.y[V.MP, :]
                self.simulations[self.obs_names.index("Cry_mRNA"), :, i] = sol.y[V.MC, :]
                self.simulations[self.obs_names.index("Bmal1_mRNA"), :, i] = sol.y[V.MB, :]

    def set_data(self):
        pass

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return []
