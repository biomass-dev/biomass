from typing import Dict, List

from biomass.analysis.reaction import is_duplicate


class ReactionNetwork(object):
    """
    Reaction indices grouped according to biological processes.
    This is used for sensitivity analysis (target='reaction').
    """

    def __init__(self) -> None:
        self.reactions: Dict[str, List[int]] = {}

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
