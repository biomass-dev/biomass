from .name2idx import C, V
from .reaction_network import ReactionNetwork


class DifferentialEquation(ReactionNetwork):
    """Kinetic equations comprising the computational model of the MAPK cascade."""

    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        # Rate equation
        v = self.flux(t, y, x)

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM
        dydt[V.MKKK] = v[2] - v[1]
        dydt[V.MKKK_P] = v[1] - v[2]
        dydt[V.MKK] = v[6] - v[3]
        dydt[V.MKK_P] = v[3] + v[5] - v[4] - v[6]
        dydt[V.MKK_PP] = v[4] - v[5]
        dydt[V.MAPK] = v[10] - v[7]
        dydt[V.MAPK_P] = v[7] + v[9] - v[8] - v[10]
        dydt[V.MAPK_PP] = v[8] - v[9]

        return dydt


def param_values():
    # Parameter values
    x = [0] * C.NUM
    x[C.V1] = 2.5
    x[C.n] = 1
    x[C.KI] = 9
    x[C.K1] = 10
    x[C.V2] = 0.25
    x[C.K2] = 8
    x[C.k3] = 0.025
    x[C.K3] = 15
    x[C.k4] = 0.025
    x[C.K4] = 15
    x[C.V5] = 0.75
    x[C.K5] = 15
    x[C.V6] = 0.75
    x[C.K6] = 15
    x[C.k7] = 0.025
    x[C.K7] = 15
    x[C.k8] = 0.025
    x[C.K8] = 15
    x[C.V9] = 0.5
    x[C.K9] = 15
    x[C.V10] = 0.5
    x[C.K10] = 15

    return x


def initial_values():
    # Value of the initial condition
    y0 = [0] * V.NUM

    y0[V.MKKK] = 90
    y0[V.MKKK_P] = 10
    y0[V.MKK] = 280
    y0[V.MKK_P] = 10
    y0[V.MKK_PP] = 10
    y0[V.MAPK] = 280
    y0[V.MAPK_P] = 10
    y0[V.MAPK_PP] = 10

    return y0
