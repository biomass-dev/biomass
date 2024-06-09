from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "IL13",
    "Rec",
    "Rec_i",
    "IL13_Rec",
    "p_IL13_Rec",
    "p_IL13_Rec_i",
    "JAK2",
    "pJAK2",
    "SHP1",
    "STAT5",
    "pSTAT5",
    "SOCS3mRNA",
    "DecoyR",
    "IL13_DecoyR",
    "SOCS3",
    "CD274mRNA",
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
