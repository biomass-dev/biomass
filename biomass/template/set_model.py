from .name2idx import C, V


class DifferentialEquation(object):
    def __init__(self, perturbation):
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        """Kinetic equations"""
        # v : flux vector
        # v = {}

        # if self.perturbation:
        #    for i, dv in self.perturbation.items():
        #        v[i] = v[i] * dv

        dydt = [0] * V.NUM

        return dydt


def param_values():
    """Parameter values"""
    x = [0] * C.NUM

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM

    return y0