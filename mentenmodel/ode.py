from re import X
from .name2idx import C, V

class DifferentialEquation(object):
    def __init__(self, pertubation) -> None:
        super(DifferentialEquation, self).__init__()
        self.pertubation = pertubation

    def diffeq(self, t, y, *x):
        dydt = [0] * V.NUM
        dydt[V.S] = - x[C.kf] * y[V.S] * y[V.E] + x[C.kr] * y[V.ES]
        dydt[V.E] = - x[C.kf] * y[V.S] * y[V.E] + x[C.kr] * y[V.ES] + x[C.kcat] * y[V.ES]
        dydt[V.ES] = x[C.kf] * y[V.S] * y[V.E] - x[C.kr] * y[V.ES] - x[C.kcat] * y[V.ES]
        dydt[V.P] = x[C.kcat] * y[V.ES]
        return dydt

def param_values():
    x = [0] * C.NUM
    x[C.kf] = 0.1
    x[C.kr] = 0.1
    x[C.kcat] = 0.1

    return x

def initial_values():
    y0 = [0] * V.NUM

    y0[V.S] = 100
    y0[V.E] = 1

    return y0
