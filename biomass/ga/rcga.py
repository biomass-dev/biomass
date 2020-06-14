import numpy as np


class UnimodalNormalDistributionXover(object):
    """ - UNDX: 
            Ono I., Kita H., Kobayashi S.. A robust real-coded genetic 
            algorithm using unimodal normal distribution crossover augmented by 
            uniform crossover: effects of self-adaptation of crossover 
            probabilities. in Proceedings of the 1st Annual Conference on 
            Genetic and Evolutionary Computation 1, 496–503 (1999).
        - MGG:
            Sato, H., Ono, I. & Kobayashi, S. A new generation alternation model 
            of genetic algorithms and its assesment. J. Jpn Soc. Artif. Intell. 
            12, 734–744 (1997).
    """
    def __init__(self, obj_func, n_population, n_children, n_gene):
        self.obj_func = obj_func
        self.n_population = n_population
        self.n_children = n_children
        self.n_gene = n_gene

    def _undx(self, parents):
        """Unimodal Normal Distribution Xover
        """
        child = np.empty(self.n_gene+1)
        ALPHA = 0.5
        BETA = 0.35 / (self.n_gene**0.5)

        d1 = np.linalg.norm(
            parents[1, :self.n_gene] - parents[0, :self.n_gene]
        )
        d2 = np.linalg.norm(
            (
                parents[2, :self.n_gene] - parents[0, :self.n_gene]
            ) - (
                np.dot(
                    (
                        parents[2, :self.n_gene]-parents[0, :self.n_gene]
                    ), (
                        parents[1, :self.n_gene]-parents[0, :self.n_gene]
                    )
                ) / (
                    d1 ** 2
                )
            ) * (
                parents[1, :self.n_gene]-parents[0, :self.n_gene]
            )
        )
        e1 = parents[0, :self.n_gene] / d1

        t = np.random.normal(scale=BETA, size=self.n_gene) * d2
        t = t - np.dot(t, e1) * e1
        t = t + np.random.normal(scale=ALPHA) * d1 * e1

        child[:self.n_gene] = t + \
            (parents[0, :self.n_gene] + parents[1, :self.n_gene]) / 2.

        return child

    def _get_new_child(self, parents):
        """
        MAXITER = 100
        for _ in range(MAXITER):
            child = _undx(parents, self.n_gene)
            if 0. <= np.min(child[:self.n_gene]) and \
                    np.max(child[:self.n_gene]) <= 1.:
                break
        else:
            child[:self.n_gene] = np.clip(child[:self.n_gene], 0., 1.)
        """
        child = self._undx(parents)
        child[:self.n_gene] = np.clip(child[:self.n_gene], 0., 1.)
        child[-1] = self.obj_func(child[:self.n_gene])

        return child

    @staticmethod
    def _rank_selection(n_family):
        ranking = np.repeat(
            np.arange(1, n_family), np.arange(1, n_family)[-1::-1]
        )
        # np.random.shuffle(ranking)
        idx = np.random.randint(len(ranking))

        return ranking[idx]

    def mgg_alternation(self, population):
        ip = [None for _ in range(3)]
        ip[:2] = np.random.choice(self.n_population, 2, replace=False)
        idx = [True] * self.n_population
        idx[ip[0]] = False
        idx[ip[1]] = False

        children = np.empty((self.n_children, self.n_gene+1))

        for i in range(self.n_children):
            ip[2] = np.random.choice(np.arange(self.n_population)[idx])
            children[i, :] = self._get_new_child(population[ip, :])
        family = np.empty((self.n_children+2, self.n_gene+1))
        family[:self.n_children, :] = children
        family[-2, :] = population[ip[0], :]
        family[-1, :] = population[ip[1], :]
        family = family[np.argsort(family[:, -1]), :]
        # Elite
        population[ip[0], :] = family[0, :]
        # Rank-based Roulette Selection
        ic1 = self._rank_selection(self.n_children+2)
        population[ip[1], :] = family[ic1, :]

        population = population[np.argsort(population[:, -1]), :]

        return population


class DistanceIndependentDiversityControl(object):
    """ DIDC:
            Kimura, S. & Konagaya, A. A Genetic Algorithm with Distance 
            Independent Diversity Control for High Dimensional Function 
            Optimization. J. Japanese Soc. Artif. Intell. 18, 193–202 (2003).
    """
    def __init__(self, obj_func, n_population, n_children, n_gene):
        self.obj_func = obj_func
        self.n_population = n_population
        self.n_children = n_children
        self.n_gene = n_gene
        #
        self.n_children_for_endx = 10

    def _endx(self, parents):
        """Extended Normal Distribution Xover
        """
        ALPHA = (1.-2*0.35**2)**0.5/2.
        BETA = 0.35/(self.n_gene-1)**0.5

        child = np.empty(self.n_gene+1)

        t1 = (parents[1, :self.n_gene]-parents[0, :self.n_gene]) / 2.
        t2 = np.random.normal(scale=ALPHA) * (
            parents[1, :self.n_gene] - parents[0, :self.n_gene]
        )
        t3 = np.sum(
            np.random.normal(scale=BETA, size=self.n_gene)[:, np.newaxis]
            * (
                parents[2:, :self.n_gene] - (
                    np.sum(parents[2:, :self.n_gene], axis=0) / self.n_gene
                )
            ), axis=0
        )
        child[:self.n_gene] = t1 + t2 + t3

        return child

    def _xover(self, parents):
        """
        MAXITER = 100
        for _ in range(MAXITER):
            child = _endx(parents, self.n_gene)
            if 0. <= np.min(child[:self.n_gene]) and \
                    np.max(child[:self.n_gene]) <= 1.:
                break
        else:
            child[:self.n_gene] = np.clip(child[:self.n_gene], 0., 1.)
        """
        child = self._endx(parents)
        child[:self.n_gene] = np.clip(child[:self.n_gene], 0., 1.)
        child[-1] = np.inf  # assigns the worst objective value to the children.

        return child

    def converging(self, ip, population):
        children = np.empty((self.n_children_for_endx, self.n_gene+1))

        for i in range(self.n_children_for_endx):
            ip[2:] = np.random.choice(
                self.n_population, self.n_gene, replace=False
            )
            children[i, :] = self._xover(population[ip, :])

        family = np.empty((self.n_children_for_endx+2, self.n_gene+1))
        family[:self.n_children_for_endx, :] = children
        family[-2, :] = population[ip[0], :]
        family[-1, :] = population[ip[1], :]

        family = family[np.argsort(family[:, -1]), :]
        # Best, either of parents
        population[ip[0], :] = family[0, :]
        # Random
        population[ip[1], :] = \
            family[
                np.random.randint(
                    low=1, high=self.n_children_for_endx+2, dtype=np.int
                ), :]

        if not np.isfinite(population[ip[1], -1]):
            population[ip[1], -1] = \
                self.obj_func(population[ip[1], :self.n_gene])
        population = population[np.argsort(population[:, -1]), :]

        return population

    def _ndm(self, parents):
        """Normal Distribution Mutation
        """
        GAMMA = 0.35/self.n_gene**0.5

        child = np.empty(self.n_gene+1)

        t2 = np.sum(
            np.random.normal(scale=GAMMA, size=self.n_gene+1)[:, np.newaxis]
            * (
                parents[1:, :self.n_gene] - (
                    np.sum(parents[1:, :self.n_gene], axis=0) / (self.n_gene+1)
                )
            ), axis=0
        )
        child[:self.n_gene] = parents[0, :self.n_gene] + t2

        return child

    def _mutation(self, parents):
        """
        MAXITER = 100
        for _ in range(MAXITER):
            child = _ndm(parents, self.n_gene)
            if 0. <= np.min(child[:self.n_gene]) and \
                    np.max(child[:self.n_gene]) <= 1.:
                break
        else:
            child[:self.n_gene] = np.clip(child[:self.n_gene], 0., 1.)
        """
        child = self._ndm(parents)
        child[:self.n_gene] = np.clip(child[:self.n_gene], 0., 1.)
        child[-1] = self.obj_func(child[:self.n_gene])

        return child

    def local_search(self, ip, population):
        idx = [True] * self.n_population
        idx[ip[0]] = False

        children = np.empty((self.n_children, self.n_gene+1))

        for i in range(self.n_children):
            ip[1:] = np.random.choice(
                np.arange(self.n_population)[idx], self.n_gene+1, replace=False
            )
            children[i, :] = self._mutation(population[ip, :])
        family = np.empty((self.n_children+1, self.n_gene+1))
        family[:self.n_children, :] = children
        family[-1, :] = population[ip[0], :]
        family = family[np.argsort(family[:, -1]), :]
        population[ip[0], :] = family[0, :]  # Elite
        population = population[np.argsort(population[:, -1]), :]

        return population
