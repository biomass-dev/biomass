from .name2idx import C, V
from .reaction_network import ReactionNetwork


class DifferentialEquation(ReactionNetwork):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        """Kinetic equations"""
        v = self.flux(t, y, x)

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM

        dydt[V.IL13] = 0.0
        dydt[V.Rec] = -v[1] - v[2]
        dydt[V.Rec_i] = v[2]
        dydt[V.IL13_Rec] = v[1] - v[3]
        dydt[V.p_IL13_Rec] = v[3] - v[4]
        dydt[V.p_IL13_Rec_i] = v[4] - v[5]
        dydt[V.JAK2] = -v[6] - v[7] + v[8]
        dydt[V.pJAK2] = v[6] + v[7] - v[8]
        dydt[V.SHP1] = 0.0
        dydt[V.STAT5] = -v[9] + v[10]
        dydt[V.pSTAT5] = v[9] - v[10]
        dydt[V.SOCS3mRNA] = v[11]
        dydt[V.DecoyR] = -v[12]
        dydt[V.IL13_DecoyR] = v[12]
        dydt[V.SOCS3] = v[13] - v[14]
        dydt[V.CD274mRNA] = v[15]

        return dydt


def param_values():
    x = [0] * C.NUM

    x[C.Kon_IL13Rec] = 0.00342
    x[C.Rec_phosphorylation] = 999.63100
    x[C.pRec_intern] = 0.15254
    x[C.pRec_degradation] = 0.17292
    x[C.Rec_intern] = 0.10335
    x[C.Rec_recycle] = 0.00136
    x[C.JAK2_phosphorylation] = 0.15706
    x[C.pJAK2_dephosphorylation] = 0.00062
    x[C.STAT5_phosphorylation] = 0.03826
    x[C.pSTAT5_dephosphorylation] = 0.00034
    x[C.SOCS3mRNA_production] = 0.00216
    x[C.DecoyR_binding] = 0.00012
    x[C.JAK2_p_inhibition] = 0.01683
    x[C.SOCS3_translation] = 11.90860
    x[C.SOCS3_accumulation] = 3.70803
    x[C.SOCS3_degradation] = 0.04292
    x[C.CD274mRNA_production] = 0.00008

    x[C.scale_pJAK2] = 1.39040
    x[C.scale_pIL4Ra] = 1.88700
    x[C.scale_IL13_cell] = 5.56750
    x[C.scale_SOCS3mRNA] = 17.66990
    x[C.scale_CD274mRNA] = 2.48547

    return x


def initial_values():
    y0 = [0] * V.NUM

    y0[V.Rec] = 1.3
    y0[V.JAK2] = 2.8
    y0[V.SHP1] = 91.0
    y0[V.STAT5] = 165.0
    y0[V.DecoyR] = 0.34
    y0[V.Rec_i] = 113.19400

    return y0
