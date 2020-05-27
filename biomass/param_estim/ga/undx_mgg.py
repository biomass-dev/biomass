import numpy as np

from biomass.models import objective


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

    child[:n_gene] = t + (parents[0, :n_gene]+parents[1, :n_gene])/2.

    return child


def _get_new_child(parents, n_gene, search_rgn):
    """
    MAXITER = 100
    for _ in range(MAXITER):
        child = _undx(parents, n_gene)
        if 0. <= np.min(child[:n_gene]) and np.max(child[:n_gene]) <= 1.:
            break
    else:
        child[:n_gene] = np.clip(child[:n_gene], 0., 1.)
    """
    child = _undx(parents, n_gene)
    child[:n_gene] = np.clip(child[:n_gene], 0., 1.)

    child[-1] = objective(child[:n_gene], search_rgn)

    return child


def _rank_selection(n_family):
    ranking = np.repeat(
        np.arange(1, n_family), np.arange(1, n_family)[-1::-1]
    )
    # np.random.shuffle(ranking)
    idx = np.random.randint(len(ranking))

    return ranking[idx]


def mgg_alternation(population, n_population, n_children, n_gene, search_rgn):
    ip = [None for _ in range(3)]
    ip[:2] = np.random.choice(n_population, 2, replace=False)
    idx = [True] * n_population
    idx[ip[0]] = False
    idx[ip[1]] = False

    children = np.empty((n_children, n_gene+1))

    for i in range(n_children):
        ip[2] = np.random.choice(np.arange(n_population)[idx])
        children[i, :] = _get_new_child(
            population[ip, :], n_gene, search_rgn
        )
    family = np.empty((n_children+2, n_gene+1))
    family[:n_children, :] = children
    family[-2, :] = population[ip[0], :]
    family[-1, :] = population[ip[1], :]

    family = family[np.argsort(family[:, -1]), :]
    # Elite
    population[ip[0], :] = family[0, :]
    # Rank-based Roulette Selection
    ic1 = _rank_selection(n_children+2)
    population[ip[1], :] = family[ic1, :]

    population = population[np.argsort(population[:, -1]), :]

    return population
