from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    """
    Reaction indices grouped according to biological processes.
    This is used for sensitivity analysis (``target``='reaction').
    """

    def __init__(self) -> None:
        super(ReactionNetwork, self).__init__()
        self.reactions: Dict[str, List[int]] = {}

    @staticmethod
    def flux(t, y, x):
        """
        Flux vector.
        """

        v = {}

        return v
