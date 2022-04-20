from .name2idx import C, V
from .reaction_network import ReactionNetwork


class DifferentialEquation(ReactionNetwork):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    @staticmethod
    def _heaviside(x):

        return 1 * (x > 0)

    # Refined Model
    def diffeq(self, t, y, *x):

        v = self.flux(t, y, x)

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        Cd = 0.65
        tRb = 5
        tC1 = 1
        Skp2 = 1
        Cdt2 = 1

        ### ALGEBRAIC EQUATIONS
        Rb = y[V.uRb] - y[V.tE2f] + y[V.E2f]
        pRb = tRb - y[V.uRb]
        RbE2f = y[V.tE2f] - y[V.E2f]

        Ce = y[V.tCe] - y[V.CeP21]
        Ca = y[V.tCa] - y[V.CaP21]
        E1 = y[V.tE1] - y[V.E1C1]
        pC1 = tC1 - y[V.C1] - y[V.E1C1]
        P21 = y[V.tP21] - y[V.CeP21] - y[V.CaP21] - y[V.iPcna] - y[V.iRc]
        # tPcna = y[V.aPcna] + y[V.iPcna] + y[V.aRc] + y[V.iRc]

        ### REACTION RATES
        rSyE2f = x[C.kSyE2f] + x[C.kSyE2fE2f] * y[V.E2f] / (x[C.jSyE2f] + y[V.E2f])
        rPhRb = x[C.kPhRbCd] * Cd + x[C.kPhRbCe] * Ce + x[C.kPhRbCa] * Ca
        rPhRc = x[C.kPhRc] * (Ce + Ca) ** x[C.n] / (x[C.jCy] ** x[C.n] + (Ce + Ca) ** x[C.n])
        rPhC1 = x[C.kPhC1] + x[C.kPhC1Ce] * Ce + x[C.kPhC1Ca] * Ca

        rDeP21 = x[C.kDeP21] + x[C.kDeP21Cy] * Skp2 * (Ce + Ca) + x[C.kDeP21aRc] * Cdt2 * y[V.aRc]
        rDeCe = x[C.kDeCe] + x[C.kDeCeCa] * Ca
        rDeCa = x[C.kDeCa] + x[C.kDeCaC1] * y[V.C1]

        rDeP53 = x[C.kDeP53] / (x[C.jP53] + y[V.Dam])
        rReDam = x[C.kReDam] + x[C.kReDamP53] * y[V.P53] / (x[C.jDam] + y[V.Dam])

        rDsRc = self._heaviside(y[V.Dna] - 1)

        ### MODEL STATES
        dydt = [0] * V.NUM

        dydt[V.uRb] = -rPhRb * y[V.uRb] + x[C.kDpRb] * pRb
        dydt[V.tE2f] = rSyE2f - x[C.kDeE2f] * y[V.tE2f]
        dydt[V.E2f] = (
            rSyE2f
            - x[C.kDeE2f] * y[V.E2f]
            - x[C.kAsRbE2f] * Rb * y[V.E2f]
            + (x[C.kDsRbE2f] + rPhRb) * RbE2f
        )

        dydt[V.tE1] = x[C.kSyE1] * y[V.E2f] - x[C.kDeE1C1] * y[V.E1C1] - x[C.kDeE1] * E1
        dydt[V.tP21] = x[C.kSyP21] + x[C.kSyP21P53] * y[V.P53] - rDeP21 * y[V.tP21]
        dydt[V.tCe] = x[C.kSyCe] * y[V.E2f] - rDeCe * y[V.tCe]
        dydt[V.tCa] = x[C.kSyCa] * y[V.E2f] - rDeCa * y[V.tCa]

        dydt[V.CeP21] = x[C.kAsCyP21] * P21 * Ce - (x[C.kDsCyP21] + rDeCe + rDeP21) * y[V.CeP21]
        dydt[V.CaP21] = x[C.kAsCyP21] * P21 * Ca - (x[C.kDsCyP21] + rDeCa + rDeP21) * y[V.CaP21]

        dydt[V.C1] = (
            -rPhC1 * y[V.C1]
            + x[C.kDpC1] * pC1
            - x[C.kAsE1C1] * E1 * y[V.C1]
            + (x[C.kDsE1C1] + x[C.kDeE1C1]) * y[V.E1C1]
        )
        dydt[V.E1C1] = (
            -rPhC1 * y[V.E1C1]
            + x[C.kAsE1C1] * E1 * y[V.C1]
            - (x[C.kDsE1C1] + x[C.kDeE1C1]) * y[V.E1C1]
        )

        dydt[V.aPcna] = (
            x[C.kImPc]
            - x[C.kAsPcP21] * P21 * y[V.aPcna]
            + (x[C.kDsPcP21] + rDeP21) * y[V.iPcna]
            - (x[C.kExPc] + x[C.kAsRcPc] * y[V.pRc]) * y[V.aPcna]
            + (x[C.kDsRcPc] + rDsRc) * y[V.aRc]
        )
        dydt[V.iPcna] = (
            x[C.kAsPcP21] * P21 * y[V.aPcna]
            - (x[C.kDsPcP21] + rDeP21) * y[V.iPcna]
            - (x[C.kExPc] + x[C.kAsRcPc] * y[V.pRc]) * y[V.iPcna]
            + (x[C.kDsRcPc] + rDsRc) * y[V.iRc]
        )

        dydt[V.Rc] = -rPhRc * y[V.Rc] + x[C.kDpRc] * y[V.pRc] - rDsRc * y[V.Rc]
        dydt[V.pRc] = (
            rPhRc * y[V.Rc]
            - x[C.kDpRc] * y[V.pRc]
            - x[C.kAsRcPc] * (y[V.aPcna] + y[V.iPcna]) * y[V.pRc]
            + x[C.kDsRcPc] * (y[V.aRc] + y[V.iRc])
            - rDsRc * y[V.pRc]
        )
        dydt[V.aRc] = (
            -x[C.kAsPcP21] * P21 * y[V.aRc]
            + (x[C.kDsPcP21] + rDeP21) * y[V.iRc]
            + x[C.kAsRcPc] * y[V.aPcna] * y[V.pRc]
            - x[C.kDsRcPc] * y[V.aRc]
            - rDsRc * y[V.aRc]
        )
        dydt[V.iRc] = (
            x[C.kAsPcP21] * P21 * y[V.aRc]
            - (x[C.kDsPcP21] + rDeP21) * y[V.iRc]
            + x[C.kAsRcPc] * y[V.iPcna] * y[V.pRc]
            - x[C.kDsRcPc] * y[V.iRc]
            - rDsRc * y[V.iRc]
        )

        dydt[V.Dna] = x[C.kSyDna] * y[V.aRc]

        dydt[V.P53] = x[C.kSyP53] - rDeP53 * y[V.P53]
        dydt[V.Dam] = x[C.kGeDam] + x[C.kGeDamArc] * y[V.aRc] - rReDam * y[V.Dam]

        dydt[V.Pr] = x[C.kSyPr] - (x[C.kDePr] + x[C.kDeCaC1] * y[V.C1]) * y[V.Pr]

        return dydt


def param_values():

    x = [0] * C.NUM

    x[C.Bg] = 0.05
    x[C.kSyE2f] = 0.03
    x[C.kSyE2fE2f] = 0.04
    x[C.jSyE2f] = 0.2
    x[C.kAsRbE2f] = 5.0
    x[C.kDsRbE2f] = 0.005
    x[C.kDeE2f] = 0.05
    x[C.kPhRbCd] = 0.2
    x[C.kPhRbCe] = 0.3
    x[C.kPhRbCa] = 0.3
    x[C.kDpRb] = 0.05
    x[C.kSyE1] = 0.005
    x[C.kDeE1C1] = 0.005
    x[C.kDeE1] = 5.0e-4
    x[C.kPhC1] = 0.0
    x[C.kPhC1Ce] = 0.01
    x[C.kPhC1Ca] = 1.0
    x[C.kDpC1] = 0.05
    x[C.kAsE1C1] = 10.0
    x[C.kDsE1C1] = 0.01
    x[C.kSyP21] = 0.002
    x[C.kSyP21P53] = 0.008
    x[C.kDeP21] = 0.0025
    x[C.kDeP21Cy] = 0.007
    x[C.kDeP21aRc] = 1.0
    x[C.kSyCe] = 0.01
    x[C.kSyCa] = 0.02
    x[C.kAsCyP21] = 1.0
    x[C.kDsCyP21] = 0.05
    x[C.kDeCe] = 0.004
    x[C.kDeCa] = 0.01
    x[C.kDeCeCa] = 0.015
    x[C.kDeCaC1] = 2.0
    x[C.kImPc] = 0.003
    x[C.kExPc] = 0.006
    x[C.kPhRc] = 0.1
    x[C.kDpRc] = 0.05
    x[C.jCy] = 1.8
    x[C.n] = 6.0
    x[C.kAsRcPc] = 0.01
    x[C.kDsRcPc] = 0.001
    x[C.kAsPcP21] = 100.0
    x[C.kDsPcP21] = 0.01
    x[C.kSyDna] = 0.0093
    x[C.kSyP53] = 0.05
    x[C.kDeP53] = 0.05
    x[C.jP53] = 0.01
    x[C.kGeDam] = 0.001
    x[C.kGeDamArc] = 0.012
    x[C.kReDam] = 0.001
    x[C.kReDamP53] = 0.005
    x[C.jDam] = 0.5
    x[C.kSyPr] = 0.01
    x[C.kDePr] = 1.0e-4

    return x


def initial_values():

    y0 = [0] * V.NUM

    y0[V.tP21] = 0.6
    y0[V.aPcna] = 0.5
    y0[V.Rc] = 1
    y0[V.tCe] = 0.5
    y0[V.tCa] = 1.2
    y0[V.Pr] = 0.5

    return y0
