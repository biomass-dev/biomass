from typing import NamedTuple

import numpy as np

from .name2idx import C, V
from .reaction_network import ReactionNetwork


class KineticInfo(NamedTuple):
    reactants: tuple
    products: tuple
    modifiers: tuple
    rate: str


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

        return dydt


def param_values():
    """Parameter values"""
    x = [1] * C.NUM

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM

    return y0


def kinetics():
    """Kinetics information"""
    kinetics = None

    return kinetics


def stoichiometric_matrix():
    """Stoichiometric matrix"""
    stoichiometry = None

    return stoichiometry
