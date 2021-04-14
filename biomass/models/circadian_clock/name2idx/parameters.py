from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "ck1",
    "ck2",
    "ck3",
    "ck4",
    "ck5",
    "ck6",
    "ck7",
    "ck8",
    "cKAP",
    "cKAC",
    "cKIB",
    "ckdmb",
    "ckdmc",
    "ckdmp",
    "ckdn",
    "ckdnc",
    "cKd",
    "cKdp",
    "cKp",
    "cKmB",
    "cKmC",
    "cKmP",
    "cksB",
    "cksC",
    "cksP",
    "cm",
    "cn",
    "cV1B",
    "cV1C",
    "cV1P",
    "cV1PC",
    "cV2B",
    "cV2C",
    "cV2P",
    "cV2PC",
    "cV3B",
    "cV3PC",
    "cV4B",
    "cV4PC",
    "cVphos",
    "cvdBC",
    "cvdBN",
    "cvdCC",
    "cvdIN",
    "cvdPC",
    "cvdPCC",
    "cvdPCN",
    "cvmB",
    "cvmC",
    "cvmP",
    "cvsB",
    "cvsC",
    "cvsP",
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
