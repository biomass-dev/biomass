def _group():
    """
    Group reactions according to biological processes
    """
    process = []

    # ERK_activation
    process.append([i for i in range(1, 7)])

    # ERK_dephosphorylation_by_DUSP
    process.append([i for i in range(47, 57)])

    # ERK_transport
    process.append([i for i in range(7, 10)])

    # RSK_activation
    process.append([24, 25])

    # RSK_transport
    process.append([26])

    # Elk1_activation
    process.append([29, 30])

    # CREB_activation
    process.append([27, 28])

    # dusp_production_etc
    process.append([i for i in range(10, 14)])

    # DUSP_transport
    process.append([18, 19])

    # DUSP_stabilization
    process.append([14, 15, 20, 21])

    # DUSP_degradation
    process.append([16, 17, 22, 23])

    # cfos_production_etc
    process.append([i for i in range(31, 35)])

    # cFos_transport
    process.append([40, 41])

    # cFos_stabilization
    process.append([35, 36, 37, 42, 43, 44])

    # cFos_degradation
    process.append([38, 39, 45, 46])

    # Feedback_from_F
    process.append([i for i in range(57, 64)])

    return process


def _is_duplicate(biological_processes):
    all_proc = sum(biological_processes, [])
    duplicate_reaction = [i for i in set(all_proc) if all_proc.count(i) > 1]
    if not duplicate_reaction:
        return False
    else:
        which_process = []
        for i in duplicate_reaction:
            for j, process in enumerate(biological_processes):
                if i in process:
                    which_process.append(j)
            raise ValueError(
                'Duplicate reaction:{:d} found in process:{}.'.format(
                    i, which_process
                )
            )

def rxn2proc():
    """
    Get n_reaction and indices
    """
    biological_processes = _group()
    if not _is_duplicate(biological_processes):
        n_reaction = 1
        for rxn in biological_processes:
            n_reaction += len(rxn)
        sort_idx = [0] * n_reaction
        left_end = 0
        for i, process in enumerate(biological_processes):
            for j, k in enumerate(process):
                if i != 0 and j == 0:
                    left_end += len(biological_processes[i-1])
                sort_idx[left_end+j] = k
                
        return biological_processes, n_reaction, sort_idx