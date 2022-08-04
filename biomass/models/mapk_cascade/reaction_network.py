from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    def __init__(self) -> None:
        super(ReactionNetwork, self).__init__()
        self.reactions: Dict[str, List[int]] = {
            "all": [i for i in range(1, 11)],
        }

    @staticmethod
    def flux(t, y, x) -> dict:

        v = {}
        v[1] = (
            x[C.V1]
            * y[V.MKKK]
            / ((1 + (y[V.MAPK_PP] / x[C.KI]) ** x[C.n]) * (x[C.K1] + y[V.MKKK]))
        )
        v[2] = x[C.V2] * y[V.MKKK_P] / (x[C.K2] + y[V.MKKK_P])
        v[3] = x[C.k3] * y[V.MKKK_P] * y[V.MKK] / (x[C.K3] + y[V.MKK])
        v[4] = x[C.k4] * y[V.MKKK_P] * y[V.MKK_P] / (x[C.K3] + y[V.MKK_P])
        v[5] = x[C.V5] * y[V.MKK_PP] / (x[C.K5] + y[V.MKK_PP])
        v[6] = x[C.V6] * y[V.MKK_P] / (x[C.K6] + y[V.MKK_P])
        v[7] = x[C.k7] * y[V.MKK_PP] * y[V.MAPK] / (x[C.K7] + y[V.MAPK])
        v[8] = x[C.k8] * y[V.MKK_PP] * y[V.MAPK_P] / (x[C.K8] + y[V.MAPK_P])
        v[9] = x[C.V9] * y[V.MAPK_PP] / (x[C.K9] + y[V.MAPK_PP])
        v[10] = x[C.V10] * y[V.MAPK_P] / (x[C.K10] + y[V.MAPK_P])

        return v
