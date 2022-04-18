from typing import Dict, List

from biomass.analysis.reaction import is_duplicate

from .name2idx import C, V


class ReactionNetwork(object):
    def __init__(self) -> None:
        """
        Reaction indices grouped according to biological processes.
        This is used for sensitivity analysis (target='reaction').
        """
        super(ReactionNetwork, self).__init__()

        self.reactions: Dict[str, List[int]] = {
            "all": [i for i in range(1, 114)],
        }

    def group(self):
        """
        Group reactions according to biological processes
        """
        biological_processes = []
        for process, indices in self.reactions.items():
            if not isinstance(indices, list):
                raise TypeError("Use list for reaction indices in {}".format(process))
            biological_processes.append(indices)

        if not is_duplicate(self.reactions, biological_processes):
            return biological_processes

    @staticmethod
    def flux(t, y, x) -> dict:
        """
        Rate equations in the model.

        Parameters
        ----------
        t : float
            Time point.
        y : ndarray
            Concentration vector.
        x : tuple
            Parameter values.

        Returns
        -------
        v : dict
            Flux vector.
        """
        v = {}

        return v
