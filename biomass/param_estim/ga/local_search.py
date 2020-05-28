import numpy as np

from biomass.models import objective


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


def _mutation(parents, n_gene):
    """
    MAXITER = 100
    for _ in range(MAXITER):
        child = _ndm(parents, n_gene)
        if 0. <= np.min(child[:n_gene]) and np.max(child[:n_gene]) <= 1.:
            break
    else:
        child[:n_gene] = np.clip(child[:n_gene], 0., 1.)
    """
    child = _ndm(parents, n_gene)
    child[:n_gene] = np.clip(child[:n_gene], 0., 1.)

    child[-1] = objective(child[:n_gene])

    return child


def local_search(ip, population, n_population, n_children, n_gene):
    idx = [True]*n_population
    idx[ip[0]] = False

    children = np.empty((n_children, n_gene+1))

    for i in range(n_children):
        ip[1:] = np.random.choice(
            np.arange(n_population)[idx], n_gene+1, replace=False
        )
        children[i, :] = _mutation(
            population[ip, :], n_gene
        )
    family = np.empty((n_children+1, n_gene+1))
    family[:n_children, :] = children
    family[-1, :] = population[ip[0], :]

    family = family[np.argsort(family[:, -1]), :]

    population[ip[0], :] = family[0, :]  # Elite

    population = population[np.argsort(population[:, -1]), :]

    return ip, population