import numpy as np


class ReactionNetwork(object):
    """
    Reaction indices grouped according to biological processes.
    This is used for sensitivity analysis (target='reaction').
    """

    reactions = {}

    def _is_duplicate(self, biological_processes):
        reaction_indices = np.sum(biological_processes, axis=0) if len(self.reactions) > 1 else biological_processes[0]
        duplicate_reaction = [i for i in set(reaction_indices) if reaction_indices.count(i) > 1]
        if not duplicate_reaction:
            return False
        else:
            which_process = []
            for reaction_index in duplicate_reaction:
                for process, indices in self.reactions.items():
                    if reaction_index in indices:
                        which_process.append(process)
            raise ValueError("Duplicate reaction: {:d} found in {}.".format(reaction_index, ", ".join(which_process)))

    def group(self):
        """
        Group reactions according to biological processes
        """
        for process, indices in self.reactions.items():
            if not isinstance(indices, list):
                raise TypeError("Use list for reaction indices in {}".format(process))
        biological_processes = []
        for process, indices in self.reactions.items():
            biological_processes.append(indices)

        if not self._is_duplicate(biological_processes):
            return biological_processes