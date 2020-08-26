import time
import numpy as np

from biomass.exec_model import ExecModel
from .rcga import (UnimodalNormalDistributionXover,
                   DistanceIndependentDiversityControl)

class GeneticAlgorithmContinue(ExecModel):
    def __init__(self, model, max_generation, allowable_error, p0_bounds):
        super().__init__(model)
        self.search_rgn = self.sp.get_region()
        self.n_population = int(5*self.search_rgn.shape[1])
        self.n_children = 50
        self.n_gene = self.search_rgn.shape[1]
        self.max_generation = max_generation
        self.allowable_error = allowable_error
        self.p0_bounds = p0_bounds

        if self.n_population < self.n_gene + 2:
            raise ValueError(
                'self.n_population must be larger than {:d}'.format(
                    self.n_gene + 2
                )
            )

    def run(self, nth_paramset):

        np.random.seed(
            time.time_ns()*nth_paramset % 2**32
        )
        (best_indiv, best_fitness) = self._ga_v2_continue(nth_paramset)

    def _set_continue(self, nth_paramset):
        best_generation = np.load(
            self.model_path + '/out/{:d}/generation.npy'.format(nth_paramset)
        )
        best_indiv = np.load(
            self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                nth_paramset, int(best_generation)
            )
        )
        population = np.full((self.n_population, self.n_gene+1), np.inf)

        with open(
                self.model_path + '/out/{:d}/'
                'optimization.log'.format(nth_paramset), mode='a') as f:
            f.write(
                '\n########################################' +
                '\n############### Continue ###############' +
                '\n########################################' +
                '\nGenerating the initial population. . .\n'
            )
        for i in range(self.n_population):
            while not np.isfinite(population[i, -1]):
                population[i, :self.n_gene] = \
                    self._encode_bestIndivVal2randGene(best_indiv)
                population[i, :self.n_gene] = \
                    np.clip(population[i, :self.n_gene], 0., 1.)
                population[i, -1] = self.obj_func(population[i, :self.n_gene])
            with open(
                    self.model_path + '/out/{:d}/'
                    'optimization.log'.format(nth_paramset), mode='a') as f:
                f.write(
                    '{:d} / {:d}\n'.format(i + 1, self.n_population)
                )
        population = population[np.argsort(population[:, -1]), :]

        return population

    def _encode_bestIndivVal2randGene(self, best_indiv):
        rand_gene = (
            np.log10(
                best_indiv * 10**(
                    np.random.rand(len(best_indiv))
                    * np.log10(self.p0_bounds[1] / self.p0_bounds[0])
                    + np.log10(self.p0_bounds[0])
                )
            ) - self.search_rgn[0, :]
        ) / (self.search_rgn[1, :] - self.search_rgn[0, :])

        return rand_gene

    def _ga_v1_continue(self, nth_paramset):
        undx = UnimodalNormalDistributionXover(
            self.obj_func, self.n_population, self.n_children, self.n_gene
        )
        count_num = np.load(
            self.model_path + '/out/{:d}/count_num.npy'.format(nth_paramset)
        )
        best_generation = np.load(
            self.model_path + '/out/{:d}/generation.npy'.format(nth_paramset)
        )
        best_indiv = np.load(
            self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                nth_paramset, int(best_generation)
            )
        )
        best_indiv_gene = self.sp.val2gene(best_indiv)
        best_fitness = self.obj_func(best_indiv_gene)

        if self.max_generation <= count_num:
            raise ValueError(
                'max_generation should be larger than {:d}'.format(
                    int(count_num)
                )
            )
        population = self._set_continue(nth_paramset)
        if best_fitness < population[0, -1]:
            population[0, :self.n_gene] = best_indiv_gene
            population[0, -1] = best_fitness
        else:
            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            best_fitness = population[0, -1]
            np.save(
                self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                    nth_paramset, int(count_num) + 1
                ), best_indiv
            )
        with open(
                self.model_path + '/out/{:d}/'
                'optimization.log'.format(nth_paramset), mode='a') as f:
            f.write(
                '\n----------------------------------------\n\n' +
                'Generation{:d}: Best Fitness = {:e}\n'.format(
                    int(count_num) + 1, best_fitness
                )
            )
        if population[0, -1] <= self.allowable_error:
            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            best_fitness = population[0, -1]

            return best_indiv, best_fitness

        generation = 1 + int(count_num)
        while generation < self.max_generation:
            population = undx.mgg_alternation(population)
            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            if population[0, -1] < best_fitness:
                np.save(
                    self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                        nth_paramset, generation + 1
                    ), best_indiv
                )
                np.save(
                    self.model_path + '/out/{:d}/generation.npy'.format(
                        nth_paramset
                    ), generation + 1
                )
                np.save(
                    self.model_path + '/out/{:d}/best_fitness'.format(
                        nth_paramset
                    ), best_fitness
                )
            best_fitness = population[0, -1]

            np.save(
                self.model_path + '/out/{:d}/count_num.npy'.format(
                    nth_paramset
                ), generation + 1
            )
            with open(
                    self.model_path + '/out/{:d}/'
                    'optimization.log'.format(nth_paramset), mode='a') as f:
                f.write(
                    'Generation{:d}: Best Fitness = {:e}\n'.format(
                        generation + 1, best_fitness
                    )
                )
            if population[0, -1] <= self.allowable_error:
                best_indiv = self.sp.gene2val(population[0, :self.n_gene])
                best_fitness = population[0, -1]

                return best_indiv, best_fitness

            generation += 1

        best_indiv = self.sp.gene2val(population[0, :self.n_gene])
        best_fitness = population[0, -1]

        return best_indiv, best_fitness

    def _ga_v2_continue(self, nth_paramset):
        didc = DistanceIndependentDiversityControl(
            self.obj_func, self.n_population, self.n_children, self.n_gene
        )
        n_iter = 1
        n0 = np.empty(3*self.n_population)

        count_num = np.load(
            self.model_path + '/out/{:d}/count_num.npy'.format(nth_paramset)
        )
        best_generation = np.load(
            self.model_path + '/out/{:d}/generation.npy'.format(nth_paramset)
        )
        best_indiv = np.load(
            self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                nth_paramset, int(best_generation)
            )
        )
        best_indiv_gene = self.sp.val2gene(best_indiv)
        best_fitness = self.obj_func(best_indiv_gene)

        if self.max_generation <= count_num:
            raise ValueError(
                'max_generation should be larger than {:d}'.format(
                    int(count_num)
                )
            )

        population = self._set_continue(nth_paramset)
        if best_fitness < population[0, -1]:
            population[0, :self.n_gene] = best_indiv_gene
            population[0, -1] = best_fitness
        else:
            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            best_fitness = population[0, -1]
            np.save(
                self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                    nth_paramset, int(count_num) + 1
                ), best_indiv
            )
        with open(
                self.model_path + '/out/{:d}/'
                'optimization.log'.format(nth_paramset), mode='a') as f:
            f.write(
                '\n----------------------------------------\n\n' +
                'Generation{:d}: Best Fitness = {:e}\n'.format(
                    int(count_num) + 1, best_fitness
                )
            )
        n0[0] = population[0, -1]

        if population[0, -1] <= self.allowable_error:
            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            best_fitness = population[0, -1]

            return best_indiv, best_fitness

        generation = 1 + int(count_num)
        while generation < self.max_generation:
            ip = np.random.choice(
                self.n_population, self.n_gene+2, replace=False
            )
            population = didc.converging(ip, population)
            population = didc.local_search(ip, population)
            for _ in range(n_iter-1):
                ip = np.random.choice(
                    self.n_population, self.n_gene+2, replace=False
                )
                population = didc.converging(ip, population)
            if generation % len(n0) == len(n0) - 1:
                n0[-1] = population[0, -1]
                if n0[0] == n0[-1]:
                    n_iter *= 2
                else:
                    n_iter = 1
            else:
                n0[generation % len(n0)] = population[0, -1]

            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            if population[0, -1] < best_fitness:
                np.save(
                    self.model_path + '/out/{:d}/generation.npy'.format(
                        nth_paramset
                    ), generation + 1
                )
                np.save(
                    self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                        nth_paramset, generation + 1
                    ), best_indiv
                )
                np.save(
                    self.model_path + '/out/{:d}/best_fitness'.format(
                        nth_paramset
                    ), best_fitness
                )
            best_fitness = population[0, -1]
            np.save(
                self.model_path + '/out/{:d}/count_num.npy'.format(
                    nth_paramset
                ), generation + 1
            )
            with open(
                    self.model_path + '/out/{:d}/'
                    'optimization.log'.format(nth_paramset), mode='a') as f:
                f.write(
                    'Generation{:d}: Best Fitness = {:e}\n'.format(
                        generation + 1, best_fitness
                    )
                )
            if population[0, -1] <= self.allowable_error:
                best_indiv = self.sp.gene2val(population[0, :self.n_gene])
                best_fitness = population[0, -1]

                return best_indiv, best_fitness

            generation += 1

        best_indiv = self.sp.gene2val(population[0, :self.n_gene])
        best_fitness = population[0, -1]

        return best_indiv, best_fitness