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
    def _compute_objval_cos(sim_data, exp_data):
        """Return Cosine distance"""
        return cosine(sim_data, exp_data)

    @staticmethod
    def _diff_sim_and_exp(sim_matrix, exp_dict, exp_timepoint, conditions, sim_norm_max):
        sim_val = []
        exp_val = []

        for idx, condition in enumerate(conditions):
            if condition in exp_dict.keys():
                sim_val.extend(sim_matrix[idx, list(map(int, exp_timepoint[condition]))])
                exp_val.extend(exp_dict[condition])

        return np.array(sim_val) / sim_norm_max, np.array(exp_val)

    def objective(self, indiv, *args):
        """Define an objective function to be minimized"""
        if len(args) == 0:
            (x, y0) = self.update(indiv)
        elif len(args) == 1:
            raise ValueError("not enough values to unpack (expected 2, got 1)")
        elif len(args) == 2:
            (x, y0) = args
        else:
            raise ValueError("too many values to unpack (expected 2)")

        self.set_data()

        if self.simulate(x, y0) is None:
            error = np.zeros(len(self.obs_names))
            for i, obs_name in enumerate(self.obs_names):
                if self.experiments[i] is not None:
                    error[i] = self._compute_objval_rss(
                        *self._diff_sim_and_exp(
                            self.simulations[i],
                            self.experiments[i],
                            self.get_timepoint(obs_name),
                            self.conditions,
                            sim_norm_max=(
                                1
                                if not self.normalization
                                else (
                                    np.max(
                                        self.simulations[
                                            self.obs_names.index(obs_name),
                                            [
                                                self.conditions.index(c)
                                                for c in (
                                                    self.normalization[obs_name]["condition"]
                                                    if self.normalization[obs_name]["condition"]
                                                    else self.conditions
                                                )
                                            ],
                                            self.normalization[obs_name]["timepoint"],
                                        ]
                                    )
                                    if self.normalization[obs_name]["timepoint"] is not None
                                    else np.max(
                                        self.simulations[
                                            self.obs_names.index(obs_name),
                                            [
                                                self.conditions.index(c)
                                                for c in (
                                                    self.normalization[obs_name]["condition"]
                                                    if self.normalization[obs_name]["condition"]
                                                    else self.conditions
                                                )
                                            ],
                                        ]
                                    )
                                )
                            ),
                        )
                    )
            """
            error = np.zeros(16)

            norm_max = np.max(self.simulations[self.obs_names.index('Phosphorylated_MEKc')])
            error[0] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_MEKc'), self.t2, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_MEKc')]['EGF']
            )
            error[1] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_MEKc'), self.t2, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_MEKc')]['HRG']
            )

            norm_max = np.max(self.simulations[self.obs_names.index('Phosphorylated_ERKc')])
            error[2] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_ERKc'), self.t2, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_ERKc')]['EGF']
            )
            error[3] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_ERKc'), self.t2, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_ERKc')]['HRG']
            )

            norm_max = np.max(self.simulations[self.obs_names.index('Phosphorylated_RSKw')])
            error[4] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_RSKw'), self.t2, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_RSKw')]['EGF']
            )
            error[5] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_RSKw'), self.t2, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_RSKw')]['HRG']
            )

            norm_max = np.max(self.simulations[self.obs_names.index('Phosphorylated_CREBw')])
            error[6] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_CREBw'), self.t3, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_CREBw')]['EGF']
            )
            error[7] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_CREBw'), self.t3, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_CREBw')]['HRG']
            )

            norm_max = np.max(self.simulations[self.obs_names.index('dusp_mRNA')])
            error[8] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('dusp_mRNA'), self.t5, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('dusp_mRNA')]['EGF']
            )
            error[9] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('dusp_mRNA'), self.t5, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('dusp_mRNA')]['HRG']
            )

            norm_max = np.max(self.simulations[self.obs_names.index('cfos_mRNA')])
            error[10] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('cfos_mRNA'), self.t4, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('cfos_mRNA')]['EGF']
            )
            error[11] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('cfos_mRNA'), self.t4, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('cfos_mRNA')]['HRG']
            )

            norm_max = np.max(self.simulations[self.obs_names.index('cFos_Protein')])
            error[12] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('cFos_Protein'), self.t5, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('cFos_Protein')]['EGF']
            )
            error[13] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('cFos_Protein'), self.t5, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('cFos_Protein')]['HRG']
            )

            norm_max = np.max(self.simulations[self.obs_names.index('Phosphorylated_cFos')])
            error[14] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_cFos'), self.t2, self.conditions.index('EGF')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_cFos')]['EGF']
            )
            error[15] = self._compute_objval_rss(
                self.simulations[self.obs_names.index('Phosphorylated_cFos'), self.t2, self.conditions.index('HRG')]/norm_max,
                self.experiments[self.obs_names.index('Phosphorylated_cFos')]['HRG']
            )
            """
            return np.sum(error)  # < 1e12
        else:
            return 1e12
