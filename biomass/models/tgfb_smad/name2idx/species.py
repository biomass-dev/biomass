from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "TGFb",
    "Rec",
    "TGFb_pRec",
    "S2",
    "S3",
    "S4",
    "ppS2_ppS2_ppS2",
    "ppS3_ppS3_ppS3",
    "S4_S4_S4",
    "pS2",
    "pS3",
    "ppS2",
    "ppS3",
    "ppS2_ppS2_S4",
    "ppS2_ppS2_ppS3",
    "ppS2_ppS3_ppS3",
    "ppS3_ppS3_S4",
    "ppS2_ppS3_S4",
    "ppS3_S4_S4",
    "ppS2_S4_S4",
    "gene",
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
