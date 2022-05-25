from .name2idx import C, V
from .reaction_network import ReactionNetwork


class DifferentialEquation(ReactionNetwork):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    # Refined Model
    def diffeq(self, t, y, *x):

        v = self.flux(t, y, x)

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM

        CycE = y[V.CycET] - y[V.CycEp27]
        CycA = y[V.CycAT] - y[V.CycAp27]

        Vdp27 = x[C.kd27] + (x[C.kd27e] * CycE + x[C.kd27a] * CycA) * y[V.Skp2]
        Vdcyce = (
            x[C.kdcyce]
            + x[C.kdcycee] * CycE / (1.0 + x[C.Inhibitor])
            + x[C.kdcycea] * CycA / (1.0 + x[C.Inhibitor])
        )
        Vdcyca = x[C.kdcyca] + x[C.kdcycac1] * y[V.Cdh1]
        Vdskp2 = x[C.kdskp2] + x[C.kdskp2c1] * y[V.Cdh1]

        Vicdh1 = x[C.kicdh1e] * CycE / (1.0 + x[C.Inhibitor]) + x[C.kicdh1a] * CycA / (
            1.0 + x[C.Inhibitor]
        )

        # protein
        dydt[V.p27T] = x[C.ks27] - Vdp27 * y[V.p27T]
        dydt[V.Skp2] = x[C.ksskp2] - Vdskp2 * y[V.Skp2]
        dydt[V.CycET] = x[C.kscyce] - Vdcyce * y[V.CycET]
        dydt[V.CycAT] = x[C.kscyca] - Vdcyca * y[V.CycAT]
        dydt[V.Emi1T] = x[C.ksemi1] - x[C.kdemi1] * y[V.Emi1T]

        dydt[V.CycEp27] = (
            x[C.kasse] * (y[V.CycET] - y[V.CycEp27]) * (y[V.p27T] - y[V.CycAp27] - y[V.CycEp27])
            - (x[C.kdise] + Vdp27 + Vdcyce) * y[V.CycEp27]
        )

        dydt[V.CycAp27] = (
            x[C.kassa] * (y[V.CycAT] - y[V.CycAp27]) * (y[V.p27T] - y[V.CycAp27] - y[V.CycEp27])
            - (x[C.kdisa] + Vdp27 + Vdcyca) * y[V.CycAp27]
        )

        dydt[V.EmiC] = (
            x[C.kasec] * (x[C.Cdh1T] - y[V.EmiC]) * (y[V.Emi1T] - y[V.EmiC])
            - (x[C.kdiec] + x[C.kdemi1]) * y[V.EmiC]
        )

        dydt[V.Cdh1dp] = x[C.kacdh1] * (x[C.Cdh1T] - y[V.Cdh1dp]) - Vicdh1 * y[V.Cdh1dp]

        dydt[V.Cdh1] = (
            (x[C.kdiec] + x[C.kdemi1]) * (y[V.Cdh1dp] - y[V.Cdh1])
            - x[C.kasec] * y[V.Cdh1] * (y[V.Emi1T] - y[V.EmiC])
            + x[C.kacdh1] * (x[C.Cdh1T] - y[V.EmiC] - y[V.Cdh1])
            - Vicdh1 * y[V.Cdh1]
        )

        return dydt


def param_values():

    x = [0] * C.NUM

    ## CYCE SYNTHESISx[C.DEGRADATION AND P27 BINDING/DISSOCIATION:
    x[C.kscyce] = 0.003
    x[C.kdcyce] = 0.001
    x[C.kdcycee] = 0.0001
    x[C.kdcycea] = 0.03
    x[C.kasse] = 1
    x[C.kdise] = 0.02
    ## CYCA SYNTHESISx[C.DEGRADATION AND P27 BINDING/DISSOCIATION:
    x[C.kscyca] = 0.0025
    x[C.kdcyca] = 0.002
    x[C.kdcycac1] = 0.4
    x[C.kassa] = 1
    x[C.kdisa] = 0.02
    ## P27 SYNTHESIS AND DEGRADATION:
    x[C.ks27] = 0.008
    x[C.kd27] = 0.004
    x[C.kd27e] = 2
    x[C.kd27a] = 2
    ## EMI1 SYNTHESIS AND DEGRADATION:
    x[C.ksemi1] = 0.003
    x[C.kdemi1] = 0.001
    ## CDH1 REGULATION:
    x[C.Cdh1T] = 1
    x[C.kacdh1] = 0.02
    x[C.kicdh1e] = 0.07
    x[C.kicdh1a] = 0.2
    x[C.kasec] = 2
    x[C.kdiec] = 0.02
    ## SKP2 SYNTHESIS AND DEGRADATION:
    x[C.ksskp2] = 0.004
    x[C.kdskp2] = 0.002
    x[C.kdskp2c1] = 0.2
    ## CDK INHIBITOR
    x[C.Inhibitor] = 0

    return x


def initial_values():

    y0 = [0] * V.NUM

    y0[V.Cdh1dp] = 1.0
    y0[V.Cdh1] = 1.0

    return y0
