"""
References
----------
- DIDC:
    Kimura, S. & Konagaya, A. A Genetic Algorithm with Distance
    Independent Diversity Control for High Dimensional Function
    Optimization. J. Japanese Soc. Artif. Intell. 18, 193–202 (2003).

- AGLSDC:
    Kimura S, Nakakuki T, Kirita S, Okada M. AGLSDC: A Genetic Local Search
    Suitable for Parallel Computation. SICE J Control Meas Syst Integr 2012;
    4: 105–113.
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from typing import Callable
from dataclasses import dataclass, field


@dataclass(frozen=True)
class RealCodedGeneticAlgorithm(object):
    obj_func: Callable[[np.ndarray], float]
    n_population: int
    n_gene: int
    n_children: int
    maxiter: int
    workers: int
    n_children_for_endx: int = field(default=10, init=False)

    def _xover(self, parents: np.ndarray) -> np.ndarray:
        """Extended Normal Distribution Xover"""
        ALPHA = (1.0 - 2 * 0.35 ** 2) ** 0.5 / 2.0
        BETA = 0.35 / (self.n_gene - 1) ** 0.5

        child = np.empty(self.n_gene + 1)

        t1 = (parents[1, : self.n_gene] - parents[0, : self.n_gene]) / 2.0
        t2 = np.random.normal(scale=ALPHA) * (parents[1, : self.n_gene] - parents[0, : self.n_gene])
        t3 = np.sum(
            np.random.normal(scale=BETA, size=self.n_gene)[:, np.newaxis]
            * (parents[2:, : self.n_gene] - (np.sum(parents[2:, : self.n_gene], axis=0) / self.n_gene)),
            axis=0,
        )
        child[: self.n_gene] = t1 + t2 + t3
        child[: self.n_gene] = np.clip(child[: self.n_gene], 0.0, 1.0)
        child[-1] = 1e12  # assigns the worst objective value to the children.

        return child

    def converging(self, ip: np.ndarray, population: np.ndarray) -> np.ndarray:
        children = np.empty((self.n_children_for_endx, self.n_gene + 1))

        for i in range(self.n_children_for_endx):
            ip[2:] = np.random.choice(self.n_population, self.n_gene, replace=False)
            children[i, :] = self._xover(population[ip, :])

        family = np.empty((self.n_children_for_endx + 2, self.n_gene + 1))
        family[: self.n_children_for_endx, :] = children
        family[-2, :] = population[ip[0], :]
        family[-1, :] = population[ip[1], :]

        family = family[np.argsort(family[:, -1]), :]
        # Best, either of parents
        population[ip[0], :] = family[0, :]
        # Random
        random_child_idx = np.random.randint(low=1, high=self.n_children_for_endx + 2, dtype=np.int)
        population[ip[1], :] = family[random_child_idx, :]

        if 1e12 <= population[ip[1], -1] or np.isnan(population[ip[1], -1]):
            population[ip[1], -1] = self.obj_func(population[ip[1], : self.n_gene])

        population = population[np.argsort(population[:, -1]), :]

        return population

    def _mutation(self, parents: np.ndarray) -> np.ndarray:
        """Normal Distribution Mutation"""
        GAMMA = 0.35 / self.n_gene ** 0.5

        child = np.empty(self.n_gene + 1)

        t2 = np.sum(
            np.random.normal(scale=GAMMA, size=self.n_gene + 1)[:, np.newaxis]
            * (parents[1:, : self.n_gene] - (np.sum(parents[1:, : self.n_gene], axis=0) / (self.n_gene + 1))),
            axis=0,
        )
        child[: self.n_gene] = parents[0, : self.n_gene] + t2
        child[: self.n_gene] = np.clip(child[: self.n_gene], 0.0, 1.0)
        child[-1] = self.obj_func(child[: self.n_gene])

        return child

    def local_search(self, ip: np.ndarray, population: np.ndarray, method: str) -> np.ndarray:
        """
        Apply the local search method to the best individual in a family
        consisting of the two parents, i.e., p1 and p2, and their children.
        """
        if method == "mutation":
            idx = [True] * self.n_population
            idx[ip[0]] = False
            children = np.empty((self.n_children, self.n_gene + 1))
            for i in range(self.n_children):
                ip[1:] = np.random.choice(np.arange(self.n_population)[idx], self.n_gene + 1, replace=False)
                children[i, :] = self._mutation(population[ip, :])
            family = np.empty((self.n_children + 1, self.n_gene + 1))
            family[: self.n_children, :] = children
            family[-1, :] = population[ip[0], :]
            family = family[np.argsort(family[:, -1]), :]
            population[ip[0], :] = family[0, :]  # Elite
        elif method == "powell":
            lower = np.min(population[ip, : self.n_gene], axis=0)
            upper = np.max(population[ip, : self.n_gene], axis=0)
            direc = np.identity(self.n_gene) * 0.3 * (upper - lower)
            res = minimize(
                self.obj_func,
                population[ip[0], : self.n_gene],
                method="Powell",
                bounds=tuple(zip(lower, upper)),
                callback=lambda xk: True
                if self.obj_func(xk) < self.obj_func(population[ip[0], : self.n_gene])
                else False,
                options={
                    "xtol": 1.0,
                    "ftol": 1.0,
                    "maxiter": self.maxiter,
                    "maxfev": 100 * self.n_gene,
                    "direc": direc,
                },
            )
            obj_val = self.obj_func(res.x)
            if obj_val < self.obj_func(population[ip[0], : self.n_gene]):
                population[ip[0], : self.n_gene] = res.x
                population[ip[0], -1] = obj_val
        elif method == "de":
            res = differential_evolution(
                self.obj_func,
                ((0.0, 1.0),) * self.n_gene,
                strategy="best2bin",
                mutation=0.1,
                recombination=0.9,
                maxiter=self.maxiter,
                popsize=1,
                polish=False,
                init=population[ip, : self.n_gene],
                updating="immediate" if self.workers == 1 else "deferred",
                workers=self.workers,
            )
            obj_val = self.obj_func(res.x)
            if obj_val < self.obj_func(population[ip[0], : self.n_gene]):
                population[ip[0], : self.n_gene] = res.x
                population[ip[0], -1] = obj_val

        population = population[np.argsort(population[:, -1]), :]

        return population