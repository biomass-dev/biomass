from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "MP",
    "MC",
    "MB",
    "PC",
    "CC",
    "PCP",
    "CCP",
    "PCC",
    "PCN",
    "PCCP",
    "PCNP",
    "BC",
    "BCP",
    "BN",
    "BNP",
    "IN",
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
