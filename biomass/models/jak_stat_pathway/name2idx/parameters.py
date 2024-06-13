from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "Kon_IL13Rec",
    "Rec_phosphorylation",
    "pRec_intern",
    "pRec_degradation",
    "Rec_intern",
    "Rec_recycle",
    "JAK2_phosphorylation",
    "pJAK2_dephosphorylation",
    "STAT5_phosphorylation",
    "pSTAT5_dephosphorylation",
    "SOCS3mRNA_production",
    "DecoyR_binding",
    "JAK2_p_inhibition",
    "SOCS3_translation",
    "SOCS3_accumulation",
    "SOCS3_degradation",
    "CD274mRNA_production",
    "scale_pJAK2",
    "scale_pIL4Ra",
    "scale_IL13_cell",
    "scale_SOCS3mRNA",
    "scale_CD274mRNA",
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
