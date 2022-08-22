from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    def __init__(self) -> None:
        """
        Reaction indices grouped according to biological processes.
        This is used for sensitivity analysis (target='reaction').
        """
        super(ReactionNetwork, self).__init__()

        self.reactions: Dict[str, List[int]] = {
            "ERK_activation": [i for i in range(1, 7)],
            "ERK_dephosphorylation_by_DUSP": [i for i in range(47, 57)],
            "ERK_transport": [i for i in range(7, 10)],
            "RSK_activation": [24, 25],
            "RSK_transport": [26],
            "Elk1_activation": [29, 30],
            "CREB_activation": [27, 28],
            "dusp_production_etc": [i for i in range(10, 14)],
            "DUSP_transport": [18, 19],
            "DUSP_stabilization": [14, 15, 20, 21],
            "DUSP_degradation": [16, 17, 22, 23],
            "cfos_production_etc": [i for i in range(31, 35)],
            "cFos_transport": [40, 41],
            "cFos_stabilization": [35, 36, 37, 42, 43, 44],
            "cFos_degradation": [38, 39, 45, 46],
            "Feedback_from_F": [i for i in range(57, 64)],
        }

    @staticmethod
    def flux(t, y, x) -> dict:
        """
        Rate equations in the model.

        Parameters
        ----------
        t : float
            Time point.
        y : ndarray
            Concentration vector.
        x : tuple
            Parameter values.

        Returns
        -------
        v : dict
            Flux vector.
        """
        v = {}
        # Rate equations
        v[1] = (
            x[C.V1]
            * x[C.a]
            * y[V.ppMEKc]
            * y[V.ERKc]
            / (x[C.Km1] * (1 + y[V.pERKc] / x[C.Km2]) + y[V.ERKc])
        )
        v[2] = (
            x[C.V2]
            * x[C.a]
            * y[V.ppMEKc]
            * y[V.pERKc]
            / (x[C.Km2] * (1 + y[V.ERKc] / x[C.Km1]) + y[V.pERKc])
        )
        v[3] = x[C.V3] * y[V.pERKc] / (x[C.Km3] * (1 + y[V.ppERKc] / x[C.Km4]) + y[V.pERKc])
        v[4] = x[C.V4] * y[V.ppERKc] / (x[C.Km4] * (1 + y[V.pERKc] / x[C.Km3]) + y[V.ppERKc])
        v[5] = x[C.V5] * y[V.pERKn] / (x[C.Km5] * (1 + y[V.ppERKn] / x[C.Km6]) + y[V.pERKn])
        v[6] = x[C.V6] * y[V.ppERKn] / (x[C.Km6] * (1 + y[V.pERKn] / x[C.Km5]) + y[V.ppERKn])
        v[7] = x[C.KimERK] * y[V.ERKc] - x[C.KexERK] * (x[C.Vn] / x[C.Vc]) * y[V.ERKn]
        v[8] = x[C.KimpERK] * y[V.pERKc] - x[C.KexpERK] * (x[C.Vn] / x[C.Vc]) * y[V.pERKn]
        v[9] = x[C.KimppERK] * y[V.ppERKc] - x[C.KexppERK] * (x[C.Vn] / x[C.Vc]) * y[V.ppERKn]
        v[10] = (
            x[C.V10] * y[V.ppERKn] ** x[C.n10] / (x[C.Km10] ** x[C.n10] + y[V.ppERKn] ** x[C.n10])
        )
        v[11] = x[C.p11] * y[V.PreduspmRNAn]
        v[12] = x[C.p12] * y[V.duspmRNAc]
        v[13] = x[C.p13] * y[V.duspmRNAc]
        v[14] = x[C.V14] * y[V.ppERKc] * y[V.DUSPc] / (x[C.Km14] + y[V.DUSPc])
        v[15] = x[C.V15] * y[V.pDUSPc] / (x[C.Km15] + y[V.pDUSPc])
        v[16] = x[C.p16] * y[V.DUSPc]
        v[17] = x[C.p17] * y[V.pDUSPc]
        v[18] = x[C.KimDUSP] * y[V.DUSPc] - x[C.KexDUSP] * (x[C.Vn] / x[C.Vc]) * y[V.DUSPn]
        v[19] = x[C.KimpDUSP] * y[V.pDUSPc] - x[C.KexpDUSP] * (x[C.Vn] / x[C.Vc]) * y[V.pDUSPn]
        v[20] = x[C.V20] * y[V.ppERKn] * y[V.DUSPn] / (x[C.Km20] + y[V.DUSPn])
        v[21] = x[C.V21] * y[V.pDUSPn] / (x[C.Km21] + y[V.pDUSPn])
        v[22] = x[C.p22] * y[V.DUSPn]
        v[23] = x[C.p23] * y[V.pDUSPn]
        v[24] = x[C.V24] * y[V.ppERKc] * y[V.RSKc] / (x[C.Km24] + y[V.RSKc])
        v[25] = x[C.V25] * y[V.pRSKc] / (x[C.Km25] + y[V.pRSKc])
        v[26] = x[C.KimRSK] * y[V.pRSKc] - x[C.KexRSK] * (x[C.Vn] / x[C.Vc]) * y[V.pRSKn]
        v[27] = x[C.V27] * y[V.pRSKn] * y[V.CREBn] / (x[C.Km27] + y[V.CREBn])
        v[28] = x[C.V28] * y[V.pCREBn] / (x[C.Km28] + y[V.pCREBn])
        v[29] = x[C.V29] * y[V.ppERKn] * y[V.Elk1n] / (x[C.Km29] + y[V.Elk1n])
        v[30] = x[C.V30] * y[V.pElk1n] / (x[C.Km30] + y[V.pElk1n])
        v[31] = (
            x[C.V31]
            * (y[V.pCREBn] * y[V.pElk1n]) ** x[C.n31]
            / (
                x[C.Km31] ** x[C.n31]
                + (y[V.pCREBn] * y[V.pElk1n]) ** x[C.n31]
                + (y[V.Fn] / x[C.KF31]) ** x[C.nF31]
            )
        )
        v[32] = x[C.p32] * y[V.PrecfosmRNAn]
        v[33] = x[C.p33] * y[V.cfosmRNAc]
        v[34] = x[C.p34] * y[V.cfosmRNAc]
        v[35] = x[C.V35] * y[V.ppERKc] * y[V.cFOSc] / (x[C.Km35] + y[V.cFOSc])
        v[36] = x[C.V36] * y[V.pRSKc] * y[V.cFOSc] / (x[C.Km36] + y[V.cFOSc])
        v[37] = x[C.V37] * y[V.pcFOSc] / (x[C.Km37] + y[V.pcFOSc])
        v[38] = x[C.p38] * y[V.cFOSc]
        v[39] = x[C.p39] * y[V.pcFOSc]
        v[40] = x[C.KimFOS] * y[V.cFOSc] - x[C.KexFOS] * (x[C.Vn] / x[C.Vc]) * y[V.cFOSn]
        v[41] = x[C.KimpcFOS] * y[V.pcFOSc] - x[C.KexpcFOS] * (x[C.Vn] / x[C.Vc]) * y[V.pcFOSn]
        v[42] = x[C.V42] * y[V.ppERKn] * y[V.cFOSn] / (x[C.Km42] + y[V.cFOSn])
        v[43] = x[C.V43] * y[V.pRSKn] * y[V.cFOSn] / (x[C.Km43] + y[V.cFOSn])
        v[44] = x[C.V44] * y[V.pcFOSn] / (x[C.Km44] + y[V.pcFOSn])
        v[45] = x[C.p45] * y[V.cFOSn]
        v[46] = x[C.p46] * y[V.pcFOSn]
        v[47] = x[C.p47] * y[V.DUSPn] * y[V.ppERKn] - x[C.m47] * y[V.DUSPn_ppERKn]
        v[48] = x[C.p48] * y[V.DUSPn_ppERKn]
        v[49] = x[C.p49] * y[V.DUSPn] * y[V.pERKn] - x[C.m49] * y[V.DUSPn_pERKn]
        v[50] = x[C.p50] * y[V.DUSPn_pERKn]
        v[51] = x[C.p51] * y[V.DUSPn] * y[V.ERKn] - x[C.m51] * y[V.DUSPn_ERKn]
        v[52] = x[C.p52] * y[V.pDUSPn] * y[V.ppERKn] - x[C.m52] * y[V.pDUSPn_ppERKn]
        v[53] = x[C.p53] * y[V.pDUSPn_ppERKn]
        v[54] = x[C.p54] * y[V.pDUSPn] * y[V.pERKn] - x[C.m54] * y[V.pDUSPn_pERKn]
        v[55] = x[C.p55] * y[V.pDUSPn_pERKn]
        v[56] = x[C.p56] * y[V.pDUSPn] * y[V.ERKn] - x[C.m56] * y[V.pDUSPn_ERKn]
        v[57] = (
            x[C.V57] * y[V.pcFOSn] ** x[C.n57] / (x[C.Km57] ** x[C.n57] + y[V.pcFOSn] ** x[C.n57])
        )
        v[58] = x[C.p58] * y[V.PreFmRNAn]
        v[59] = x[C.p59] * y[V.FmRNAc]
        v[60] = x[C.p60] * y[V.FmRNAc]
        v[61] = x[C.p61] * y[V.Fc]
        v[62] = x[C.KimF] * y[V.Fc] - x[C.KexF] * (x[C.Vn] / x[C.Vc]) * y[V.Fn]
        v[63] = x[C.p63] * y[V.Fn]

        return v
