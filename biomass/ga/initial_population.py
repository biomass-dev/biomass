import numpy as np


class InitialPopulation(object):
    def __init__(self, get_region, objective):
        self.get_region = get_region
        self.objective = objective

    def set_initial(self, nth_paramset, n_population, n_gene):
        population = np.full((n_population, n_gene+1), np.inf)
        with open('./out/{:d}/initpop.log'.format(nth_paramset), mode='w') as f:
            f.write(
                'Generating the initial population. . .\n'
            )
        for i in range(n_population):
            while not np.isfinite(population[i, -1]):
                population[i, :n_gene] = np.random.rand(n_gene)
                population[i, -1] = self.objective(population[i, :n_gene])
            with open('./out/{:d}/initpop.log'.format(nth_paramset), mode='a') as f:
                f.write(
                    '{:d} / {:d}\n'.format(i + 1, n_population)
                )
        population = population[np.argsort(population[:, -1]), :]

        return population

    def set_continue(self, nth_paramset, n_population, n_gene, p0_bounds):
        best_generation = np.load(
            './out/{:d}/generation.npy'.format(nth_paramset)
        )
        best_indiv = np.load(
            './out/{:d}/fit_param{:d}.npy'.format(nth_paramset, int(best_generation))
        )
        population = np.full((n_population, n_gene+1), np.inf)

        with open('./out/{:d}/initpop.log'.format(nth_paramset), mode='w') as f:
            f.write(
                'Generating the initial population. . .\n'
            )
        for i in range(n_population):
            while not np.isfinite(population[i, -1]):
                population[i, :n_gene] = self._encode_bestIndivVal2randGene(
                    best_indiv, p0_bounds
                )
                population[i, :n_gene] = np.clip(population[i, :n_gene], 0., 1.)
                population[i, -1] = self.objective(population[i, :n_gene])
            with open('./out/{:d}/initpop.log'.format(nth_paramset), mode='a') as f:
                f.write(
                    '{:d} / {:d}\n'.format(i + 1, n_population)
                )
        population = population[np.argsort(population[:, -1]), :]

        return population

    def _encode_bestIndivVal2randGene(self, best_indiv, p0_bounds):
        search_rgn = self.get_region()
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