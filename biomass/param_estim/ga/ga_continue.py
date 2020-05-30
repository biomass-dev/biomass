import time
import numpy as np

from biomass.models import objective, get_search_region, decode_gene2val
from .undx_mgg import mgg_alternation
from .converging import converging
from .local_search import local_search


def optimize_continue(nth_paramset):

    np.random.seed(
        time.time_ns()*nth_paramset % 2**32
    )
    search_rgn = get_search_region()

    max_generation = 10000
    n_population = int(5*search_rgn.shape[1])
    n_children = 50
    n_gene = search_rgn.shape[1]
    allowable_error = 0.35
    p0_bounds = [0.1, 10.0]  # [lower_bounds, upper bounds]

    (best_indiv, best_fitness) = ga_v2_continue(
        nth_paramset,
        max_generation,
        n_population,
        n_children,
        n_gene,
        allowable_error,
        p0_bounds
    )


def _encode_val2gene(indiv):
    search_rgn = get_search_region()
    indiv_gene = (
        np.log10(indiv) - search_rgn[0, :]
    ) / (
        search_rgn[1, :] - search_rgn[0, :]
    )

    return indiv_gene


def _encode_bestIndivVal2randGene(best_indiv, p0_bounds):
    search_rgn = get_search_region()
    rand_gene = (
        np.log10(
            best_indiv * 10**(
                np.random.rand(len(best_indiv))
                * np.log10(p0_bounds[1]/p0_bounds[0])
                + np.log10(p0_bounds[0])
            )
        ) - search_rgn[0, :]
    ) / (search_rgn[1, :] - search_rgn[0, :])

    return rand_gene


def _get_initial_population_continue(nth_paramset, n_population, n_gene, p0_bounds):
    best_generation = np.load(
        './out/{:d}/generation.npy'.format(
            nth_paramset
        )
    )
    best_indiv = np.load(
        './out/{:d}/fit_param{:d}.npy'.format(
            nth_paramset, int(best_generation)
        )
    )
    population = np.full((n_population, n_gene+1), np.inf)

    with open('./out/{:d}/initpop.log'.format(nth_paramset), mode='w') as f:
        f.write(
            'Generating the initial population. . .\n'
        )
    for i in range(n_population):
        while not np.isfinite(population[i, -1]):
            population[i, :n_gene] = _encode_bestIndivVal2randGene(
                best_indiv, p0_bounds
            )
            population[i, :n_gene] = np.clip(population[i, :n_gene], 0., 1.)
            population[i, -1] = objective(
                population[i, :n_gene]
            )
        with open('./out/{:d}/initpop.log'.format(nth_paramset), mode='a') as f:
            f.write(
                '{:d} / {:d}\n'.format(
                    i + 1, n_population
                )
            )
    population = population[np.argsort(population[:, -1]), :]

    return population


def ga_v1_continue(nth_paramset, max_generation, n_population, n_children, n_gene,
                    allowable_error, p0_bounds):
    count_num = np.load(
        './out/{:d}/count_num.npy'.format(
            nth_paramset
        )
    )
    best_generation = np.load(
        './out/{:d}/generation.npy'.format(
            nth_paramset
        )
    )
    best_indiv = np.load(
        './out/{:d}/fit_param{:d}.npy'.format(
            nth_paramset, int(best_generation)
        )
    )
    best_fitness = objective(
        _encode_val2gene(best_indiv)
    )
    population = _get_initial_population_continue(
        nth_paramset, n_population, n_gene, p0_bounds
    )
    if best_fitness < population[0, -1]:
        population[0, :n_gene] = (
            _encode_val2gene(best_indiv)
        )
        population[0, -1] = best_fitness
    else:
        best_indiv = decode_gene2val(
            population[0, :n_gene]
        )
        best_fitness = population[0, -1]
        np.save(
            './out/{:d}/fit_param{:d}.npy'.format(
                nth_paramset, int(count_num) + 1
            ), best_indiv
        )
    with open('./out/{:d}/out.log'.format(nth_paramset), mode='a') as f:
        f.write(
            '\n----------------------------------------\n\n' +
            'Generation{:d}: Best Fitness = {:e}\n'.format(
                int(count_num) + 1, best_fitness
            )
        )
    if population[0, -1] <= allowable_error:
        best_indiv = decode_gene2val(
            population[0, :n_gene]
        )
        best_fitness = population[0, -1]
        return best_indiv, best_fitness

    generation = 1
    while generation < max_generation:
        population = mgg_alternation(
            population, n_population, n_children, n_gene
        )
        best_indiv = decode_gene2val(
            population[0, :n_gene]
        )
        if population[0, -1] < best_fitness:
            np.save(
                './out/{:d}/fit_param{:d}.npy'.format(
                    nth_paramset, generation + int(count_num) + 1
                ), best_indiv
            )
            np.save(
                './out/{:d}/generation.npy'.format(
                    nth_paramset
                ), generation + int(count_num) + 1
            )
            np.save(
                './out/{:d}/best_fitness'.format(
                    nth_paramset
                ), best_fitness
            )
        best_fitness = population[0, -1]

        np.save(
            './out/{:d}/count_num.npy'.format(
                nth_paramset
            ), generation + int(count_num) + 1
        )
        with open('./out/{:d}/out.log'.format(nth_paramset), mode='a') as f:
            f.write(
                'Generation{:d}: Best Fitness = {:e}\n'.format(
                    generation + int(count_num) + 1, best_fitness
                )
            )
        if population[0, -1] <= allowable_error:
            best_indiv = decode_gene2val(
                population[0, :n_gene]
            )
            best_fitness = population[0, -1]
            return best_indiv, best_fitness

        generation += 1

    best_indiv = decode_gene2val(
        population[0, :n_gene]
    )
    best_fitness = population[0, -1]

    return best_indiv, best_fitness


def ga_v2_continue(nth_paramset, max_generation, n_population, n_children, n_gene,
                    allowable_error, p0_bounds):
    if n_population < n_gene + 2:
        raise ValueError(
            'n_population must be larger than {:d}'.format(
                n_gene + 2
            )
        )
    n_iter = 1
    n0 = np.empty(3*n_population)

    count_num = np.load(
        './out/{:d}/count_num.npy'.format(
            nth_paramset
        )
    )
    best_generation = np.load(
        './out/{:d}/generation.npy'.format(
            nth_paramset
        )
    )
    best_indiv = np.load(
        './out/{:d}/fit_param{:d}.npy'.format(
            nth_paramset, int(best_generation)
        )
    )
    best_fitness = objective(
        _encode_val2gene(best_indiv)
    )
    population = _get_initial_population_continue(
        nth_paramset, n_population, n_gene, p0_bounds
    )
    if best_fitness < population[0, -1]:
        population[0, :n_gene] = (
            _encode_val2gene(best_indiv)
        )
        population[0, -1] = best_fitness
    else:
        best_indiv = decode_gene2val(
            population[0, :n_gene]
        )
        best_fitness = population[0, -1]
        np.save(
            './out/{:d}/fit_param{:d}.npy'.format(
                nth_paramset, int(count_num) + 1
            ), best_indiv
        )
    with open('./out/{:d}/out.log'.format(nth_paramset), mode='a') as f:
        f.write(
            '\n----------------------------------------\n\n' +
            'Generation{:d}: Best Fitness = {:e}\n'.format(
                int(count_num) + 1, best_fitness
            )
        )
    n0[0] = population[0, -1]

    if population[0, -1] <= allowable_error:
        best_indiv = decode_gene2val(
            population[0, :n_gene]
        )
        best_fitness = population[0, -1]
        return best_indiv, best_fitness

    generation = 1
    while generation < max_generation:
        ip = np.random.choice(n_population, n_gene+2, replace=False)
        ip, population = converging(
            ip, population, n_population, n_gene
        )
        ip, population = local_search(
            ip, population, n_population, n_children, n_gene
        )
        for _ in range(n_iter-1):
            ip = np.random.choice(n_population, n_gene+2, replace=False)
            ip, population = converging(
                ip, population, n_population, n_gene
            )
        if generation % len(n0) == len(n0) - 1:
            n0[-1] = population[0, -1]
            if n0[0] == n0[-1]:
                n_iter *= 2
            else:
                n_iter = 1
        else:
            n0[generation % len(n0)] = population[0, -1]

        best_indiv = decode_gene2val(
            population[0, :n_gene]
        )
        if population[0, -1] < best_fitness:
            np.save(
                './out/{:d}/generation.npy'.format(
                    nth_paramset
                ), generation + int(count_num) + 1
            )
            np.save(
                './out/{:d}/fit_param{:d}.npy'.format(
                    nth_paramset, generation + int(count_num) + 1
                ), best_indiv
            )
            np.save(
                './out/{:d}/best_fitness'.format(
                    nth_paramset
                ), best_fitness
            )
        best_fitness = population[0, -1]

        np.save(
            './out/{:d}/count_num.npy'.format(
                nth_paramset
            ), generation + int(count_num) + 1
        )
        with open('./out/{:d}/out.log'.format(nth_paramset), mode='a') as f:
            f.write(
                'Generation{:d}: Best Fitness = {:e}\n'.format(
                    generation + int(count_num) + 1, best_fitness
                )
            )
        if population[0, -1] <= allowable_error:
            best_indiv = decode_gene2val(
                population[0, :n_gene]
            )
            best_fitness = population[0, -1]
            return best_indiv, best_fitness

        generation += 1

    best_indiv = decode_gene2val(
        population[0, :n_gene]
    )
    best_fitness = population[0, -1]

    return best_indiv, best_fitness