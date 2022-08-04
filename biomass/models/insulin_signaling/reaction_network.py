from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    def __init__(self) -> None:
        super(ReactionNetwork, self).__init__()
        self.reactions: Dict[str, List[int]] = {}

    @staticmethod
    def flux(t, y, x):

        v = {}

        v[1] = x[C.k1_synthesis] * (y[V.pro_IRcom] - y[V.IRcom])
        v[2] = x[C.k1_InsIRcom] * y[V.Ins] * y[V.IRcom] - x[C.k2_InsIRcom] * y[V.p1IRcom]
        v[3] = x[C.k1_p1IRcomDeg] * y[V.p1IRcom]
        v[4] = x[C.k1_p1IRcomPhos] * y[V.pmTOR] * y[V.p1IRcom]
        v[5] = x[C.k1_p1p2IRcomdePhos] * y[V.p1p2IRcom]
        v[6] = x[C.k1_IRcomPhos] * y[V.pmTOR] * y[V.IRcom]
        v[7] = x[C.k1_p2IRcomdePhos] * y[V.p2IRcom]
        v[8] = x[C.k1_p2IRcomDeg] * y[V.p2IRcom]
        v[9] = x[C.k1_Insp2IRcom] * y[V.Ins] * y[V.p2IRcom] - x[C.k2_Insp2IRcom] * y[V.p1p2IRcom]
        v[10] = x[C.k1_p1p2IRcomDeg] * y[V.p1p2IRcom]
        v[11] = x[C.k1_AKTPhos] * (y[V.iAKT] - y[V.pAKT]) * y[V.p1IRcom]
        v[12] = x[C.k1_pAKTdePhos] * y[V.pAKT]
        v[13] = x[C.k1_mTORPhos] * (y[V.imTOR] - y[V.pmTOR]) * y[V.pAKT]
        v[14] = x[C.k1_pmTORdePhos] * y[V.pmTOR]
        v[15] = x[C.k1_S6KPhos] * (y[V.iS6K] - y[V.pS6K]) * y[V.pmTOR]
        v[16] = x[C.k1_pS6KdePhos] * y[V.pS6K] * y[V.pX]
        v[17] = x[C.k1_XPhos] * (y[V.iX] - y[V.pX]) * y[V.pmTOR]
        v[18] = x[C.k1_pXdePhos] * y[V.pX]
        v[19] = x[C.k1_GSK3BPhos] * (y[V.iGSK3B] - y[V.pGSK3B]) * y[V.pAKT]
        v[20] = x[C.k1_pGSK3BdePhos] * y[V.pGSK3B]
        v[21] = x[C.k1_FoxO1Phos] * (y[V.iFoxO1] - y[V.pFoxO1]) * y[V.pAKT]
        v[22] = x[C.k1_pFoxO1dePhos] * y[V.pFoxO1]
        v[23] = x[C.k1_G6PaseSynthesis] * (y[V.iFoxO1] - y[V.pFoxO1])
        v[24] = x[C.k1_G6PaseDeg] * y[V.G6Pase]

        return v
