from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "k1_synthesis",
    "k1_InsIRcom",
    "k2_InsIRcom",
    "k1_p1IRcomDeg",
    "k1_p1IRcomPhos",
    "k1_p1p2IRcomdePhos",
    "k1_IRcomPhos",
    "k1_p2IRcomdePhos",
    "k1_p2IRcomDeg",
    "k1_Insp2IRcom",
    "k2_Insp2IRcom",
    "k1_p1p2IRcomDeg",
    "k1_AKTPhos",
    "k1_pAKTdePhos",
    "k1_mTORPhos",
    "k1_pmTORdePhos",
    "k1_S6KPhos",
    "k1_pS6KdePhos",
    "k1_XPhos",
    "k1_pXdePhos",
    "k1_GSK3BPhos",
    "k1_pGSK3BdePhos",
    "k1_FoxO1Phos",
    "k1_pFoxO1dePhos",
    "k1_G6PaseSynthesis",
    "k1_G6PaseDeg",
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
