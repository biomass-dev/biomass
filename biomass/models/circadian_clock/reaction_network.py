from typing import Dict, List


class ReactionNetwork(object):
    """
    Reaction indices grouped according to biological processes.
    This is used for sensitivity analysis (target='reaction').
    """

    def __init__(self) -> None:
        self.reactions: Dict[str, List[int]] = {}
