from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "Ins",
    "pro_IRcom",
    "IRcom",
    "p1IRcom",
    "p2IRcom",
    "p1p2IRcom",
    "iAKT",
    "pAKT",
    "imTOR",
    "pmTOR",
    "iX",
    "pX",
    "iS6K",
    "pS6K",
    "iGSK3B",
    "pGSK3B",
    "iFoxO1",
    "pFoxO1",
    "G6Pase",
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
