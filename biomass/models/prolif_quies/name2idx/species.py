from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "uRb",
    "tE2f",
    "E2f",
    "tE1",
    "tP21",
    "tCe",
    "tCa",
    "CeP21",
    "CaP21",
    "C1",
    "E1C1",
    "aPcna",
    "iPcna",
    "Rc",
    "pRc",
    "aRc",
    "iRc",
    "Dna",
    "P53",
    "Dam",
    "Pr",
]

NUM: int = len(NAMES)

Species = make_dataclass(
    cls_name="Species",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

V = Species(**name2idx)

del name2idx
