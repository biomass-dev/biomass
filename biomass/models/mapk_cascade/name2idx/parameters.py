from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "V1",
    "n",
    "KI",
    "K1",
    "V2",
    "K2",
    "k3",
    "K3",
    "k4",
    "K4",
    "V5",
    "K5",
    "V6",
    "K6",
    "k7",
    "K7",
    "k8",
    "K8",
    "V9",
    "K9",
    "V10",
    "K10",
]

NUM: int = len(NAMES)

Parameters = make_dataclass(
    cls_name="Parameters",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

C = Parameters(**name2idx)

del name2idx
