import numpy as np
from scipy.spatial.distance import cosine

from .observable import Observable
from .search_param import SearchParam

class OptimizationProblem(Observable, SearchParam):
    def __init__(self):
        super(OptimizationProblem, self).__init__()

    @property
    def bounds(self):
        """
        Lower and upper bounds on independent variables.
        """
        search_region: np.ndarray = self.get_region()
        lb = 10 ** search_region[0]
        ub = 10 ** search_region[1]
        return tuple(zip(lb, ub))

    @staticmethod
    def _compute_objval_rss(sim_data, exp_data):
        """Return Residual Sum of Squares"""
        return np.dot((sim_data - exp_data), (sim_data - exp_data))

    @staticmethod
    def _diff_sim_and_exp
