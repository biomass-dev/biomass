from typing import List


NAMES: List[str] = []

for idx, name in enumerate(NAMES):
    exec("{} = {:d}".format(name, idx))

NUM: int = len(NAMES)
