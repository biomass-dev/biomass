import numpy as np

from biomass.dynamics.solver import solve_ode

from .name2idx import C, V
from .set_model import DifferentialEquation

observables = [
    "Per_mRNA",
    "Cry_mRNA",
    "Bmal1_mRNA",
]


class NumericalSimulation(DifferentialEquation):
    """Simulate a model using scipy.integrate.solve_ivp

    Attributes
    ----------
    normalization : nested dict
        Keys for each observable
        ------------------------
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If None, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in sim.conditions will be used.

    """

    def __init__(self):
        super().__init__(perturbation={})
        self.normalization = {}

    t = range(72 + 1)

    # Experimental conditions
    conditions = ["DD"]

    simulations = np.empty((len(observables), len(t), len(conditions)))

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
                self.simulations[observables.index("Per_mRNA"), :, i] = sol.y[V.MP, :]
                self.simulations[observables.index("Cry_mRNA"), :, i] = sol.y[V.MC, :]
                self.simulations[observables.index("Bmal1_mRNA"), :, i] = sol.y[V.MB, :]


class ExperimentalData(object):
    """
    Set experimental data.

    Attributes
    ----------
    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """

    def __init__(self):
        self.experiments = [None] * len(observables)
        self.error_bars = [None] * len(observables)

    def set_data(self):
        pass

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return []
