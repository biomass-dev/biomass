import numpy as np
from scipy.spatial.distance import cosine

from biomass.model import C, V, f_params, initial_values
from biomass.observable import *
from .search_parameter import update_param
from .ga.converter import decode_gene2variable


def _compute_objval_rss(sim_data, exp_data):
    """Return Residual Sum of Squares
    """
    return np.dot((sim_data-exp_data), (sim_data-exp_data))


def _compute_objval_cos(sim_data, exp_data):
    """Return Cosine distance
    """
    return cosine(sim_data, exp_data)


def _diff_sim_and_exp(sim_matrix, exp_dict, exp_timepoint, conditions, 
                        sim_norm_max=1, exp_norm_max=1):
    sim_val = []
    exp_val = []

    for condition in conditions:
        if condition in exp_dict.keys():
            sim_val.extend(
                sim_matrix[exp_timepoint, conditions.index(condition)]
            )
            exp_val.extend(exp_dict[condition])

    return np.array(sim_val)/sim_norm_max, np.array(exp_val)/exp_norm_max


def objective(individual_gene, search_idx, search_region):
    """Define an objective function to be minimized
    """
    x = f_params()
    y0 = initial_values()

    indiv = decode_gene2variable(individual_gene, search_region)

    (x, y0) = update_param(indiv, x, y0)

    exp = ExperimentalData()
    sim = NumericalSimulation()

    if sim.simulate(x, y0) is None:
        error = np.zeros(len(observables))
        for i, _ in enumerate(observables):
            exp_t = exp.get_timepoint(i)
            norm_max = np.max(sim.simulations[i])
            if exp.experiments[i] is not None:
                error[i] = _compute_objval_rss(
                    *_diff_sim_and_exp(
                        sim.simulations[i], exp.experiments[i], exp_t, sim.conditions, 
                        sim_norm_max=norm_max, 
                    )
                )
        '''
        error = np.zeros(16)

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_MEKc')])
        error[0] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_MEKc'), exp.t2, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_MEKc')]['EGF']
        )
        error[1] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_MEKc'), exp.t2, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_MEKc')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_ERKc')])
        error[2] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_ERKc'), exp.t2, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_ERKc')]['EGF']
        )
        error[3] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_ERKc'), exp.t2, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_ERKc')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_RSKw')])
        error[4] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_RSKw'), exp.t2, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_RSKw')]['EGF']
        )
        error[5] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_RSKw'), exp.t2, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_RSKw')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_CREBw')])
        error[6] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_CREBw'), exp.t3, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_CREBw')]['EGF']
        )
        error[7] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_CREBw'), exp.t3, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_CREBw')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('dusp_mRNA')])
        error[8] = _compute_objval_rss(
            sim.simulations[observables.index('dusp_mRNA'), exp.t5, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('dusp_mRNA')]['EGF']
        )
        error[9] = _compute_objval_rss(
            sim.simulations[observables.index('dusp_mRNA'), exp.t5, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('dusp_mRNA')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('cfos_mRNA')])
        error[10] = _compute_objval_rss(
            sim.simulations[observables.index('cfos_mRNA'), exp.t4, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('cfos_mRNA')]['EGF']
        )
        error[11] = _compute_objval_rss(
            sim.simulations[observables.index('cfos_mRNA'), exp.t4, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('cfos_mRNA')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('cFos_Protein')])
        error[12] = _compute_objval_rss(
            sim.simulations[observables.index('cFos_Protein'), exp.t5, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('cFos_Protein')]['EGF']
        )
        error[13] = _compute_objval_rss(
            sim.simulations[observables.index('cFos_Protein'), exp.t5, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('cFos_Protein')]['HRG']
        )

        norm_max = np.max(sim.simulations[observables.index('Phosphorylated_cFos')])
        error[14] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_cFos'), exp.t2, sim.conditions.index('EGF')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_cFos')]['EGF']
        )
        error[15] = _compute_objval_rss(
            sim.simulations[observables.index('Phosphorylated_cFos'), exp.t2, sim.conditions.index('HRG')]/norm_max, 
            exp.experiments[observables.index('Phosphorylated_cFos')]['HRG']
        )
        '''
        return np.sum(error)
    else:
        return np.inf