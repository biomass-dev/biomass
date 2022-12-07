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
        Names of self.obs_names.

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
        self.obs_names = ["CycA", "CycA_tot", "CycE", "CycE_tot", "p27_tot"]
        self.t: range = range(901)
        self.conditions: list = ["control"]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation=None):
        if _perturbation is not None:
            self.perturbation = _perturbation
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == "control":
                pass

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("CycA_tot"), i] = sol.y[V.CycAT]
                self.simulations[self.obs_names.index("CycA"), i] = (
                    sol.y[V.CycAT] - sol.y[V.CycAp27]
                )
                self.simulations[self.obs_names.index("CycE_tot"), i] = sol.y[V.CycET]
                self.simulations[self.obs_names.index("CycE"), i] = (
                    sol.y[V.CycET] - sol.y[V.CycEp27]
                )
                self.simulations[self.obs_names.index("p27_tot"), i] = sol.y[V.p27T]
        return None

    def set_data(self):
        pass

    def get_timepoint(self, obs_name):
        if obs_name in self.obs_names:
            return []
        assert False
