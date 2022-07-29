from typing import Dict, List

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
