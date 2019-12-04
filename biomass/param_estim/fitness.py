import numpy as np
from scipy.spatial.distance import cosine

from biomass.model.name2idx import parameters as C
from biomass.model.name2idx import variables as V
from biomass.model.param_const import f_params
from biomass.model.initial_condition import initial_values
from biomass.observable import *
from .ga.converter import decode_gene2variable

def compute_objval_abs(sim_data,exp_data):  # Residual Sum of Squares

    return np.dot((sim_data-exp_data),(sim_data-exp_data))


def compute_objval_cs(sim_data,exp_data):  # Cosine similarity

    return cosine(sim_data,exp_data)


# Define an objective function to be minimized.
def objective(individual_gene,search_idx,search_region):
    # update_param
    x = f_params()
    y0 = initial_values()

    indiv = decode_gene2variable(individual_gene,search_region)

    for i,j in enumerate(search_idx[0]):
        x[j] = indiv[i]
    for i,j in enumerate(search_idx[1]):
        y0[j] = indiv[i+len(search_idx[0])]

    # constraints --------------------------------------------------------------
    x[C.V6] = x[C.V5]
    x[C.Km6] = x[C.Km5]
    x[C.KimpDUSP] = x[C.KimDUSP]
    x[C.KexpDUSP] = x[C.KexDUSP]
    x[C.KimpcFOS] = x[C.KimFOS]
    x[C.KexpcFOS] = x[C.KexFOS]
    x[C.p52] = x[C.p47]
    x[C.m52] = x[C.m47]
    x[C.p53] = x[C.p48]
    x[C.p54] = x[C.p49]
    x[C.m54] = x[C.m49]
    x[C.p55] = x[C.p50]
    x[C.p56] = x[C.p51]
    x[C.m56] = x[C.m51]
    # --------------------------------------------------------------------------

    exp = ExperimentalData()
    sim = NumericalSimulation()

    if sim.simulate(x,y0) is None:
        error = np.zeros(len(species))
        for i,_ in enumerate(species):
            exp_t = exp.get_timepoint(i)
            norm_max = np.max(sim.simulations[i])
            if exp.experiments[i] is not None:
                error[i] = compute_objval_abs(
                    *diff_sim_and_exp(
                        sim.simulations[i], exp.experiments[i],
                        exp_timepoint=exp_t, conditions=sim.conditions,
                        norm_max_sim=norm_max
                    )
                )
        '''
        error = np.zeros(16)

        norm_max = np.max(sim.simulations[species.index('Phosphorylated_MEKc')])
        error[0] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_MEKc'), exp.t2, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('Phosphorylated_MEKc')]['EGF']
        )
        error[1] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_MEKc'), exp.t2, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('Phosphorylated_MEKc')]['HRG']
        )

        norm_max = np.max(sim.simulations[species.index('Phosphorylated_ERKc')])
        error[2] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_ERKc'), exp.t2, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('Phosphorylated_ERKc')]['EGF']
        )
        error[3] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_ERKc'), exp.t2, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('Phosphorylated_ERKc')]['HRG']
        )

        norm_max = np.max(sim.simulations[species.index('Phosphorylated_RSKw')])
        error[4] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_RSKw'), exp.t2, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('Phosphorylated_RSKw')]['EGF']
        )
        error[5] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_RSKw'), exp.t2, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('Phosphorylated_RSKw')]['HRG']
        )

        norm_max = np.max(sim.simulations[species.index('Phosphorylated_CREBw')])
        error[6] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_CREBw'), exp.t3, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('Phosphorylated_CREBw')]['EGF']
        )
        error[7] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_CREBw'), exp.t3, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('Phosphorylated_CREBw')]['HRG']
        )

        norm_max = np.max(sim.simulations[species.index('dusp_mRNA')])
        error[8] = compute_objval_abs(
            sim.simulations[species.index('dusp_mRNA'), exp.t5, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('dusp_mRNA')]['EGF']
        )
        error[9] = compute_objval_abs(
            sim.simulations[species.index('dusp_mRNA'), exp.t5, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('dusp_mRNA')]['HRG']
        )

        norm_max = np.max(sim.simulations[species.index('cfos_mRNA')])
        error[10] = compute_objval_abs(
            sim.simulations[species.index('cfos_mRNA'), exp.t4, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('cfos_mRNA')]['EGF']
        )
        error[11] = compute_objval_abs(
            sim.simulations[species.index('cfos_mRNA'), exp.t4, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('cfos_mRNA')]['HRG']
        )

        norm_max = np.max(sim.simulations[species.index('cFos_Protein')])
        error[12] = compute_objval_abs(
            sim.simulations[species.index('cFos_Protein'), exp.t5, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('cFos_Protein')]['EGF']
        )
        error[13] = compute_objval_abs(
            sim.simulations[species.index('cFos_Protein'), exp.t5, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('cFos_Protein')]['HRG']
        )

        norm_max = np.max(sim.simulations[species.index('Phosphorylated_cFos')])
        error[14] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_cFos'), exp.t2, sim.conditions.index('EGF')]/norm_max,
            exp.experiments[species.index('Phosphorylated_cFos')]['EGF']
        )
        error[15] = compute_objval_abs(
            sim.simulations[species.index('Phosphorylated_cFos'), exp.t2, sim.conditions.index('HRG')]/norm_max,
            exp.experiments[species.index('Phosphorylated_cFos')]['HRG']
        )
        '''
        return np.sum(error)
    else:
        return np.inf