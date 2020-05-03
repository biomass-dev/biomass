# Num. of Rate Equations (len(v) in diffeq in model/set_model.py)
num_reaction = 64

def get_reaction_module():
    n_module = 16
    reaction_module = [None]*n_module
    # ERK_activation
    reaction_module[0] = [i for i in range(1,7)]
    # ERK_dephosphorylation_by_DUSP
    reaction_module[1] = [i for i in range(47,57)]
    # ERK_transport
    reaction_module[2] = [i for i in range(7,10)]
    # RSK_activation
    reaction_module[3] = [24,25]
    # RSK_transport
    reaction_module[4] = [26]
    # Elk1_activation
    reaction_module[5] = [29,30]
    # CREB_activation
    reaction_module[6] = [27,28]
    # dusp_production_etc
    reaction_module[7] = [i for i in range(10,14)]
    # DUSP_transport
    reaction_module[8] = [18,19]
    # DUSP_stabilization
    reaction_module[9] = [14,15,20,21]
    # DUSP_degradation
    reaction_module[10] = [16,17,22,23]
    # cfos_production_etc
    reaction_module[11] = [i for i in range(31,35)]
    # cFos_transport
    reaction_module[12] = [40,41]
    # cFos_stabilization
    reaction_module[13] = [35,36,37,42,43,44]
    # cFos_degradation
    reaction_module[14] = [38,39,45,46]
    # Feedback_from_F
    reaction_module[15] = [i for i in range(57,64)]

    return reaction_module


def get_sort_idx():
    reaction_module = get_reaction_module()
    sort_idx = [0] * num_reaction
    left_end = 0
    for i, ith_module in enumerate(reaction_module):
        for j, k in enumerate(ith_module):
            if i != 0 and j == 0:
                left_end += len(reaction_module[i-1])
            sort_idx[left_end+j] = k
            
    return sort_idx