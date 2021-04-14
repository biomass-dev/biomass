from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "ppMEKc",
    "CREBn",
    "pCREBn",
    "ERKc",
    "ERKn",
    "pERKc",
    "pERKn",
    "ppERKc",
    "ppERKn",
    "Elk1n",
    "pElk1n",
    "cFOSc",
    "cFOSn",
    "pcFOSc",
    "pcFOSn",
    "DUSPc",
    "DUSPn",
    "pDUSPc",
    "pDUSPn",
    "DUSPn_ERKn",
    "DUSPn_pERKn",
    "DUSPn_ppERKn",
    "pDUSPn_ERKn",
    "pDUSPn_pERKn",
    "pDUSPn_ppERKn",
    "RSKc",
    "pRSKc",
    "pRSKn",
    "PrecfosmRNAn",
    "PreduspmRNAn",
    "cfosmRNAc",
    "duspmRNAc",
    "Fc",
    "Fn",
    "FmRNAc",
    "PreFmRNAn",
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
