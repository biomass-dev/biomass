from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    def __init__(self) -> None:
        super(ReactionNetwork, self).__init__()
        self.reactions: Dict[str, List[int]] = {}

    @staticmethod
    def flux(t, y, x) -> dict:

        v = {}

        v[1] = y[V.Rec] * x[C.Rec_act] * y[V.TGFb]
        v[2] = y[V.TGFb_pRec] * x[C.pRec_degind]
        v[3] = x[C.k_on_222] * (y[V.ppS2] ** 3)
        v[4] = 3 * x[C.S_dephosphos] * y[V.ppS2_ppS2_ppS2]
        v[5] = x[C.k_on_333] * y[V.ppS3] ** 3
        v[6] = 3 * x[C.S_dephosphos] * y[V.ppS3_ppS3_ppS3]
        v[7] = (y[V.S4] ** 3) * x[C.k_on_444]
        v[8] = y[V.S4_S4_S4] * x[C.kdiss_SS]
        v[9] = y[V.S2] * x[C.S_phos] * y[V.TGFb_pRec]
        v[10] = x[C.S_dephosphos] * y[V.ppS2]
        v[11] = x[C.S_dephos] * y[V.pS2]
        v[12] = y[V.S3] * x[C.S_phos] * y[V.TGFb_pRec]
        v[13] = x[C.S_dephosphos] * y[V.ppS3]
        v[14] = x[C.S_dephos] * y[V.pS3]
        v[15] = x[C.k_on_223] * (y[V.ppS2] ** 2) * y[V.ppS3]
        v[16] = 2 * x[C.S_dephosphos] * y[V.ppS2_ppS2_ppS3]
        v[17] = x[C.S_dephosphos] * y[V.ppS2_ppS2_ppS3]
        v[18] = y[V.S4] * x[C.k_on_224] * (y[V.ppS2] ** 2)
        v[19] = 2 * x[C.S_dephosphos] * y[V.ppS2_ppS2_S4]
        v[20] = x[C.k_on_233] * y[V.ppS2] * (y[V.ppS3] ** 2)
        v[21] = x[C.S_dephosphos] * y[V.ppS2_ppS3_ppS3]
        v[22] = 2 * x[C.S_dephosphos] * y[V.ppS2_ppS3_ppS3]
        v[23] = y[V.S4] * x[C.k_on_334] * (y[V.ppS3] ** 2)
        v[24] = 2 * x[C.S_dephosphos] * y[V.ppS3_ppS3_S4]
        v[25] = (y[V.S4] ** 2) * x[C.k_on_244] * y[V.ppS2]
        v[26] = x[C.S_dephosphos] * y[V.ppS2_S4_S4]
        v[27] = (y[V.S4] ** 2) * x[C.k_on_344] * y[V.ppS3]
        v[28] = x[C.S_dephosphos] * y[V.ppS3_S4_S4]
        v[29] = y[V.S4] * x[C.k_on_234] * y[V.ppS2] * y[V.ppS3]
        v[30] = x[C.S_dephosphos] * y[V.ppS2_ppS3_S4]
        v[31] = x[C.S_dephosphos] * y[V.ppS2_ppS3_S4]
        v[32] = (
            x[C.gene_turn]
            + x[C.gene_act1] * y[V.ppS2_ppS3_ppS3]
            + x[C.gene_act2] * y[V.ppS2_S4_S4]
            + x[C.gene_act3] * y[V.ppS2_ppS3_S4]
        ) / (
            x[C.gene_inh1] * y[V.ppS2_ppS3_ppS3]
            + x[C.gene_inh2] * y[V.ppS2_S4_S4]
            + x[C.gene_inh3] * y[V.ppS2_ppS3_S4]
            + 1
        )
        v[33] = y[V.gene] * x[C.gene_turn]

        return v
