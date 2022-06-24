from typing import List

import numpy as np
from biomass.dynamics.solver import *

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
        self.obs_names: list = []
        self.t: range = range(101)
        self.conditions: list = []
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.conditions), len(self.t))
        )
        self.normalization: dict = {}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

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

    def set_data(self) -> None:
        pass

    def get_timepoint(self, obs_name: str) -> List[int]:
        if obs_name in self.obs_names:
            return []
        assert False
