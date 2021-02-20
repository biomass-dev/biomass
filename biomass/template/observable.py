import numpy as np

from biomass.dynamics.solver import *

from .name2idx import C, V
from .set_model import DifferentialEquation

observables = []


class NumericalSimulation(DifferentialEquation):
    """Simulate a model using scipy.integrate.ode

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

    t = range(101)  # 0, 1, 2, ..., 100

    # Experimental conditions
    conditions = []

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x: list, y0: list, _perturbation: dict = {}):
        if _perturbation:
            self.perturbation = _perturbation
        # unperturbed steady state

        for i, condition in enumerate(self.conditions):

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))

            if sol is None:
                return False
            else:
                pass


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
    def get_timepoint(obs_name: str):
        if obs_name in observables:
            return []
