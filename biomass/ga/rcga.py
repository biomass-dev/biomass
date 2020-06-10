import numpy as np


class UnimodalNormalDistributionXover(object):
    """ UNDX + MGG
    """
    def __init__(self, objective):
        self.objective = objective

    @staticmethod
    def _undx(parents, n_gene):
        """Unimodal Normal Distribution Xover
        """
        child = np.empty(n_gene+1)
        ALPHA = 0.5
        BETA = 0.35 / (n_gene**0.5)

        d1 = np.linalg.norm(parents[1, :n_gene] - parents[0, :n_gene])
        d2 = np.linalg.norm(
            (
                parents[2, :n_gene] - parents[0, :n_gene]
            ) - (
                np.dot(
                    (
                        parents[2, :n_gene]-parents[0, :n_gene]
                    ), (
                        parents[1, :n_gene]-parents[0, :n_gene]
                    )
                ) / (
                    d1 ** 2
                )
            ) * (
                parents[1, :n_gene]-parents[0, :n_gene]
            )
        )
        e1 = parents[0, :n_gene]/d1

        t = np.random.normal(scale=BETA, size=n_gene) * d2
        t = t - np.dot(t, e1) * e1
        t = t + np.random.normal(scale=ALPHA) * d1 * e1

        child[:n_gene] = t + (parents[0, :n_gene] + parents[1, :n_gene]) / 2.

        return child

    def _get_new_child(self, parents, n_gene):
        """
        MAXITER = 100
        for _ in range(MAXITER):
            child = _undx(parents, n_gene)
            if 0. <= np.min(child[:n_gene]) and np.max(child[:n_gene]) <= 1.:
                break
        else:
            child[:n_gene] = np.clip(child[:n_gene], 0., 1.)
        """
        child = self._undx(parents, n_gene)
        child[:n_gene] = np.clip(child[:n_gene], 0., 1.)
        child[-1] = self.objective(child[:n_gene])

        return child

    @staticmethod
    def _rank_selection(n_family):
        ranking = np.repeat(
            np.arange(1, n_family), np.arange(1, n_family)[-1::-1]
        )
        # np.random.shuffle(ranking)
        idx = np.random.randint(len(ranking))

        return ranking[idx]

    def mgg_alternation(self, population, n_population, n_children, n_gene):
        ip = [None for _ in range(3)]
        ip[:2] = np.random.choice(n_population, 2, replace=False)
        idx = [True] * n_population
        idx[ip[0]] = False
        idx[ip[1]] = False

        children = np.empty((n_children, n_gene+1))

        for i in range(n_children):
            ip[2] = np.random.choice(np.arange(n_population)[idx])
            children[i, :] = self._get_new_child(
                population[ip, :], n_gene
            )
        family = np.empty((n_children+2, n_gene+1))
        family[:n_children, :] = children
        family[-2, :] = population[ip[0], :]
        family[-1, :] = population[ip[1], :]

        family = family[np.argsort(family[:, -1]), :]
        # Elite
        population[ip[0], :] = family[0, :]
        # Rank-based Roulette Selection
        ic1 = self._rank_selection(n_children+2)
        population[ip[1], :] = family[ic1, :]

        population = population[np.argsort(population[:, -1]), :]

        return population


class DistanceIndependentDiversityControl(object):
    """Modified
    """
    def __init__(self, objective):
        self.objective = objective

    @staticmethod
    def _endx(parents, n_gene):
        """Extended Normal Distribution Xover
        """
        ALPHA = (1.-2*0.35**2)**0.5/2.
        BETA = 0.35/(n_gene-1)**0.5

        child = np.empty(n_gene+1)

        t1 = (parents[1, :n_gene]-parents[0, :n_gene]) / 2.
        t2 = np.random.normal(scale=ALPHA) * (
            parents[1, :n_gene] - parents[0, :n_gene]
        )
        t3 = np.sum(
            np.random.normal(scale=BETA, size=n_gene)[:, np.newaxis]
            * (
                parents[2:, :n_gene] - (
                    np.sum(parents[2:, :n_gene], axis=0) / n_gene
                )
            ), axis=0
        )
        child[:n_gene] = t1 + t2 + t3

        return child

    def _xover(self, parents, n_gene):
        """
        MAXITER = 100
        for _ in range(MAXITER):
            child = _endx(parents, n_gene)
            if 0. <= np.min(child[:n_gene]) and np.max(child[:n_gene]) <= 1.:
                break
        else:
            child[:n_gene] = np.clip(child[:n_gene], 0., 1.)
        """
        child = self._endx(parents, n_gene)
        child[:n_gene] = np.clip(child[:n_gene], 0., 1.)
        child[-1] = np.inf  # assigns the worst objective value to the children.

        return child

    def converging(self, ip, population, n_population, n_gene):
        n_children = 10
        children = np.empty((n_children, n_gene+1))

        for i in range(n_children):
            ip[2:] = np.random.choice(n_population, n_gene, replace=False)
            children[i, :] = self._xover(population[ip, :], n_gene)

        family = np.empty((n_children+2, n_gene+1))
        family[:n_children, :] = children
        family[-2, :] = population[ip[0], :]
        family[-1, :] = population[ip[1], :]

        family = family[np.argsort(family[:, -1]), :]
        # Best, either of parents
        population[ip[0], :] = family[0, :]
        # Random
        population[ip[1], :] = \
            family[np.random.randint(low=1, high=n_children+2, dtype=np.int), :]

        if not np.isfinite(population[ip[1], -1]):
            population[ip[1], -1] = self.objective(population[ip[1], :n_gene])

        population = population[np.argsort(population[:, -1]), :]

        return population

    @staticmethod
    def _ndm(parents, n_gene):
        """Normal Distribution Mutation
        """
        GAMMA = 0.35/n_gene**0.5

        child = np.empty(n_gene+1)

        t2 = np.sum(
            np.random.normal(scale=GAMMA, size=n_gene+1)[:, np.newaxis]
            * (
                parents[1:, :n_gene] - (
                    np.sum(parents[1:, :n_gene], axis=0) / (n_gene+1)
                )
            ), axis=0
        )
        child[:n_gene] = parents[0, :n_gene] + t2

        return child

    def _mutation(self, parents, n_gene):
        """
        MAXITER = 100
        for _ in range(MAXITER):
            child = _ndm(parents, n_gene)
            if 0. <= np.min(child[:n_gene]) and np.max(child[:n_gene]) <= 1.:
                break
        else:
            child[:n_gene] = np.clip(child[:n_gene], 0., 1.)
        """
        child = self._ndm(parents, n_gene)
        child[:n_gene] = np.clip(child[:n_gene], 0., 1.)

        child[-1] = self.objective(child[:n_gene])

        return child

    def local_search(self, ip, population, n_population, n_children, n_gene):
        idx = [True] * n_population
        idx[ip[0]] = False

        children = np.empty((n_children, n_gene+1))

        for i in range(n_children):
            ip[1:] = np.random.choice(
                np.arange(n_population)[idx], n_gene+1, replace=False
            )
            children[i, :] = self._mutation(
                population[ip, :], n_gene
            )
        family = np.empty((n_children+1, n_gene+1))
        family[:n_children, :] = children
        family[-1, :] = population[ip[0], :]

        family = family[np.argsort(family[:, -1]), :]

        population[ip[0], :] = family[0, :]  # Elite

        population = population[np.argsort(population[:, -1]), :]

        return population
