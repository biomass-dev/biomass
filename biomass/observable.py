import numpy as np
    
observable_names = [
    'Phosphorylated_MEKc',
    'Phosphorylated_ERKc',
    'Phosphorylated_RSKw',
    'Phosphorylated_CREBw',
    'dusp_mRNA',
    'cfos_mRNA',
    'cFos_Protein',
    'Phosphorylated_cFos'\
]

num_observables = len(observable_names)
species = dict(zip(observable_names,range(num_observables)))


def diff_sim_and_exp(
    sim_matrix,exp_dict,exp_timepoint,num_condition,
    norm_max_sim=1,norm_max_exp=1
    ):
        
    return (
        np.r_[
            [sim_matrix[exp_timepoint,i] for i in range(num_condition)]
        ].flatten()/norm_max_sim,
        np.r_[
            list(exp_dict.values())
        ].flatten()/norm_max_exp
    )