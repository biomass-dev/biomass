from typing import List


NAMES: List[str] = []

for idx, name in enumerate(NAMES):
    exec(f"{name} = {idx:d}")

NUM: int = len(NAMES)