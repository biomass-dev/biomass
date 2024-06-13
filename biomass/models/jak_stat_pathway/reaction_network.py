from typing import Dict, List

from .name2idx import C, V


class ReactionNetwork(object):
    def __init__(self) -> None:
        """
        Reaction indices grouped according to biological processes.
        This is used for sensitivity analysis (target='reaction').
        """
        super(ReactionNetwork, self).__init__()
        self.reactions: Dict[str, List[int]] = {}

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
        IL13_scale = 2.265
        v[1] = x[C.Kon_IL13Rec] * y[V.IL13] * y[V.Rec] * IL13_scale
        v[2] = x[C.Rec_intern] * y[V.Rec] - x[C.Rec_recycle] * y[V.Rec_i]
        v[3] = x[C.Rec_phosphorylation] * y[V.IL13_Rec] * y[V.pJAK2]
        v[4] = x[C.pRec_intern] * y[V.p_IL13_Rec]
        v[5] = x[C.pRec_degradation] * y[V.p_IL13_Rec_i]
        v[6] = (
            x[C.JAK2_phosphorylation]
            * y[V.IL13_Rec]
            * y[V.JAK2]
            / (1 + x[C.JAK2_p_inhibition] * y[V.SOCS3])
        )
        v[7] = (
            x[C.JAK2_phosphorylation]
            * y[V.p_IL13_Rec]
            * y[V.JAK2]
            / (1 + x[C.JAK2_p_inhibition] * y[V.SOCS3])
        )
        v[8] = x[C.pJAK2_dephosphorylation] * y[V.pJAK2] * y[V.SHP1]
        v[9] = x[C.STAT5_phosphorylation] * y[V.STAT5] * y[V.pJAK2]
        v[10] = x[C.pSTAT5_dephosphorylation] * y[V.pSTAT5] * y[V.SHP1]
        v[11] = x[C.SOCS3mRNA_production] * y[V.pSTAT5]
        v[12] = x[C.DecoyR_binding] * y[V.IL13] * y[V.DecoyR] * IL13_scale
        v[13] = (
            x[C.SOCS3_translation] * y[V.SOCS3mRNA] / (x[C.SOCS3_accumulation] + y[V.SOCS3mRNA])
        )
        v[14] = x[C.SOCS3_degradation] * y[V.SOCS3]
        v[15] = x[C.CD274mRNA_production] * y[V.pSTAT5]

        return v
