from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    ## CYCE SYNTHESISDEGRADATION AND P27 BINDING/DISSOCIATION:
    "kscyce",
    "kdcyce",
    "kdcycee",
    "kdcycea",
    "kasse",
    "kdise",
    ## CYCA SYNTHESISDEGRADATION AND P27 BINDING/DISSOCIATION:
    "kscyca",
    "kdcyca",
    "kdcycac1",
    "kassa",
    "kdisa",
    ## P27 SYNTHESIS AND DEGRADATION:
    "ks27",
    "kd27",
    "kd27e",
    "kd27a",
    ## EMI1 SYNTHESIS AND DEGRADATION:
    "ksemi1",
    "kdemi1",
    ## CDH1 REGULATION:
    "Cdh1T",
    "kacdh1",
    "kicdh1e",
    "kicdh1a",
    "kasec",
    "kdiec",
    ## SKP2 SYNTHESIS AND DEGRADATION:
    "ksskp2",
    "kdskp2",
    "kdskp2c1",
    ## CDK INHIBITOR
    "Inhibitor",
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
