from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "uptake",
    "TNF",
    "trigger_iIkk",
    "deact_TNFR",
    "deact_ppIkk",
    "deact_pnNfk",
    "act_Ikk_by_TNF",
    "act_pIkk",
    "act_Ikb_by_Ikk",
    "act_Nfk_by_Ikk",
    "act_Nfk_by_Ikk_complex",
    "act_Ikb_complex",
    "form_complex",
    "form_complex_nuc",
    "ext_nNfkIkb",
    "Vnuc",
    "split_NfkpIkb",
    "split_NfkIkb",
    "int_Nfk",
    "int_Ikb",
    "eta_int_pNfk",
    "degrad_Ikb",
    "degrad_mIkb",
    "degrad_RnaA20",
    "degrad_A20",
    "prod_Ikb",
    "prod_mIkb_by_nNfk",
    "build_RnaA20",
    "build_A20",
    "shuttle_RnaA20",
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
