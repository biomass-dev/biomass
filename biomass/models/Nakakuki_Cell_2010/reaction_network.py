import numpy as np


class ReactionNetwork(object):
    reactions = {
        'ERK_activation': [i for i in range(1, 7)],
        'ERK_dephosphorylation_by_DUSP': [i for i in range(47, 57)],
        'ERK_transport': [i for i in range(7, 10)],
        'RSK_activation': [24, 25],
        'RSK_transport': [26],
        'Elk1_activation': [29, 30],
        'CREB_activation': [27, 28],
        'dusp_production_etc': [i for i in range(10, 14)],
        'DUSP_transport': [18, 19],
        'DUSP_stabilization': [14, 15, 20, 21],
        'DUSP_degradation': [16, 17, 22, 23],
        'cfos_production_etc': [i for i in range(31, 35)],
        'cFos_transport': [40, 41],
        'cFos_stabilization': [35, 36, 37, 42, 43, 44],
        'cFos_degradation': [38, 39, 45, 46],
        'Feedback_from_F': [i for i in range(57, 64)],
    }

    def _is_duplicate(self, biological_processes):
        reaction_indices = np.sum(biological_processes, axis=0) \
            if len(self.reactions) > 1 else biological_processes[0]
        duplicate_reaction = \
            [i for i in set(reaction_indices) if reaction_indices.count(i) > 1]
        if not duplicate_reaction:
            return False
        else:
            which_process = []
            for reaction_index in duplicate_reaction:
                for process, indices in self.reactions.items():
                    if reaction_index in indices:
                        which_process.append(process)
            raise ValueError(
                'Duplicate reaction: {:d} found in {}.'.format(
                    reaction_index, which_process
                )
            )

    def group(self):
        """
        Group reactions according to biological processes
        """
        for process, indices in self.reactions.items():
            if not isinstance(indices, list):
                raise TypeError(
                    'Use list for reaction indices in {}'.format(process)
                )
        biological_processes = []
        for process, indices in self.reactions.items():
            biological_processes.append(indices)

        if not self._is_duplicate(biological_processes):
            return biological_processes