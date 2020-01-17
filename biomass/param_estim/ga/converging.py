import numpy as np

from biomass.param_estim.fitness import objective


def converging(ip, population, n_population, n_gene, search_idx, search_region):
    n_children = 10
    children = np.empty((n_children, n_gene+1))

    for i in range(n_children):
        ip[2:] = np.random.choice(n_population, n_gene, replace=False)
        children[i, :] = xover(population[ip, :], n_gene)

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

    if np.isinf(population[ip[1], -1]):
        population[ip[1], -1] = objective(
            population[ip[1], :n_gene], search_idx, search_region
        )
    population = population[np.argsort(population[:, -1]), :]

    return ip, population


def xover(parents, n_gene):
    MAXITER = 100
    for _ in range(MAXITER):
        child = endx(parents, n_gene)
        if 0. <= np.min(child[:n_gene]) and np.max(child[:n_gene]) <= 1.:
            break
    else:
        child[:n_gene] = np.clip(child[:n_gene], 0., 1.)

    child[-1] = np.inf  # assigns the worst objective value to the children.

    return child


def endx(parents, n_gene):
    """Extended Normal Distribution Xover
    """
    ALPHA = (1.-2*0.35**2)**0.5/2.
    BETA = 0.35/(n_gene-1)**0.5

    child = np.empty(n_gene+1)

    t1 = (parents[1, :n_gene]-parents[0, :n_gene])/2.
    t2 = np.random.normal(scale=ALPHA) * (parents[1, :n_gene] - parents[0, :n_gene])
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
