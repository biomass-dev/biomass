import os
import time
import numpy as np

from biomass.exec_model import ExecModel
from .rcga import (UnimodalNormalDistributionXover,
                   DistanceIndependentDiversityControl)

class GeneticAlgorithmInit(ExecModel):
    def __init__(self, model, max_generation, allowable_error):
        super().__init__(model)
        self.search_rgn = self.sp.get_region()
        self.n_population = int(5*self.search_rgn.shape[1])
        self.n_children = 50
        self.n_gene = self.search_rgn.shape[1]
        self.max_generation = max_generation
        self.allowable_error = allowable_error
    
        if self.n_population < self.n_gene + 2:
            raise ValueError(
                'n_population must be larger than {:d}'.format(
                    self.n_gene + 2
                )
            )

    def run(self, nth_paramset):
        if not os.path.isdir(self.model_path + '/out'):
            os.mkdir(self.model_path + '/out')
        try:
            files = os.listdir(
                self.model_path + '/out/{:d}'.format(
                    nth_paramset
                )
            )
            for file in files:
                if any(map(file.__contains__, ('.npy', '.log'))):
                    os.remove(
                        self.model_path + '/out/{:d}/{}'.format(
                            nth_paramset, file
                        )
                    )
        except FileNotFoundError:
            os.mkdir(
                self.model_path + '/out/{:d}'.format(
                    nth_paramset
                )
            )
        np.random.seed(
            time.time_ns()*nth_paramset % 2**32
        )
        (best_indiv, best_fitness) = self._ga_v2(nth_paramset)

    def _set_initial(self, nth_paramset):
        population = np.full((self.n_population, self.n_gene+1), np.inf)
        with open(
                self.model_path + '/out/{:d}/'
                'optimization.log'.format(nth_paramset), mode='w') as f:
            f.write(
                'Generating the initial population. . .\n'
            )
        for i in range(self.n_population):
            while not np.isfinite(population[i, -1]):
                population[i, :self.n_gene] = np.random.rand(self.n_gene)
                population[i, -1] = self.obj_func(population[i, :self.n_gene])
            with open(
                    self.model_path + '/out/{:d}/'
                    'optimization.log'.format(nth_paramset), mode='a') as f:
                f.write(
                    '{:d} / {:d}\n'.format(i + 1, self.n_population)
                )
        population = population[np.argsort(population[:, -1]), :]

        return population

    def _ga_v1(self, nth_paramset):
        undx = UnimodalNormalDistributionXover(
            self.obj_func, self.n_population, self.n_children, self.n_gene
        )
        population = self._set_initial(nth_paramset)
        with open(
                self.model_path + '/out/{:d}/'
                'optimization.log'.format(nth_paramset), mode='a') as f:
            f.write(
                '\n----------------------------------------\n\n' +
                'Generation1: Best Fitness = {:e}\n'.format(population[0, -1])
            )
        best_indiv = self.sp.gene2val(population[0, :self.n_gene])
        best_fitness = population[0, -1]

        np.save(
            self.model_path + '/out/{:d}/generation.npy'.format(
                nth_paramset
            ), 1
        )
        np.save(
            self.model_path + '/out/{:d}/fit_param1'.format(
                nth_paramset
            ), best_indiv
        )
        np.save(
            self.model_path + '/out/{:d}/best_fitness.npy'.format(
                nth_paramset
            ), best_fitness
        )
        if population[0, -1] <= self.allowable_error:
            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            best_fitness = population[0, -1]

            return best_indiv, best_fitness

        generation = 1
        while generation < self.max_generation:
            population = undx.mgg_alternation(population)
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
            best_fitness = population[0, -1]
            np.save(
                self.model_path + '/out/{:d}/best_fitness.npy'.format(
                    nth_paramset
                ), best_fitness
            )
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

        best_indiv = self.sp.gene2val(
            population[0, :self.n_gene]
        )
        best_fitness = population[0, -1]

        return best_indiv, best_fitness

    def _ga_v2(self, nth_paramset):
        """ga_v2 optimizes an objective function through the following procedure.

        1. Initialization
            As an initial population, create np individuals randomly.
            ga_v2 also represents individuals as n-dimensional real number
            vectors, where n is the dimension of the search space. Set
            Generation to 0, and set the iteration number of converging
            operations Niter to 1.

        2. Selection for reproduction
            As parents for the recombination operator, ENDX, select m
            individuals, p1, p2, . . . ,pm, without replacement from the
            population.

        3. Generation of offsprings
            Generate Nc children by applying ENDX to the selected parents.
            This algorithm assigns the worst objective value to the children.

        4. Local Search (NDM/MGG)
            Apply the local search method to the best individual in a family
            consisting of the two parents, i.e., p1 and p2, and their children.
            Note here that the children are assumed to have the worst objective
            value. Thus, whenever the objective values of the two parents have
            been actually computed in previous generations, the algorithm
            applies the local search to either of the parents. When all of the
            individuals in the family have the same objective value, on the
            other hand, the local search is applied to a randomly selected
            individual from the family.

        5. Selection for survival
            Select two individuals from the family. The first selected
            individual should be the individual with the best objective value,
            and the second should be selected randomly. Then, replace the two
            parents (p1 and p2) with the selected individuals. Note that the
            individual to which the local search has been applied in the
            previous step is always selected as the best.

        6. Application of ENDX/MGG
            To achieve a good search performance, ga_v2 optimizes a function,
            gradually narrowing the search space. For this purpose, the
            converging phase slightly converges the population by repeating the
            following procedure Niter times.

            (i) Select m individuals without replacement from the population.
                The selected individuals, expressed here as p1, p2, . . . , pm,
                are used as the parents for an extended normal distribution
                crossover (ENDX) applied in the next step.

            (ii) Generate Nc children by applying ENDX to the parents selected
                in the previous step. To reduce the computational cost, ga_v2
                forgoes any computation of the objective values of the Nc
                individuals generated here. Instead, the algorithm assigns the
                newly generated children a single objective value, one which is
                inferior to the objective values of any of the possible
                candidate solutions.

            (iii) Select two individuals from a family containing the two
                parents, i.e., p1 and p2, and their children. The first
                selected individual should be the one with the best objective
                value, and the second should be selected randomly. Then,
                replace the two parents with the selected individuals.

        7. Adaptation of Niter
            If the best individual has not improved during the last np
            generations, Niter <- 2 * Niter. Otherwise, set Niter to 1.

        8. Termination
            Stop if the halting criteria are satisfied.
            Otherwise, Generation <- Generation + 1, and return to the step 2.
        """
        didc = DistanceIndependentDiversityControl(
            self.obj_func, self.n_population, self.n_children, self.n_gene
        )
        n_iter = 1
        n0 = np.empty(3*self.n_population)

        population = self._set_initial(nth_paramset)
        n0[0] = population[0, -1]

        with open(
                self.model_path + '/out/{:d}/'
                'optimization.log'.format(nth_paramset), mode='a') as f:
            f.write(
                '\n----------------------------------------\n\n' +
                'Generation1: Best Fitness = {:e}\n'.format(population[0, -1])
            )
        best_indiv = self.sp.gene2val(population[0, :self.n_gene])
        best_fitness = population[0, -1]

        np.save(
            self.model_path + '/out/{:d}/generation.npy'.format(
                nth_paramset
            ), 1
        )
        np.save(
            self.model_path + '/out/{:d}/fit_param1.npy'.format(
                nth_paramset
            ), best_indiv
        )
        np.save(
            self.model_path + '/out/{:d}/best_fitness.npy'.format(
                nth_paramset
            ), best_fitness
        )
        if population[0, -1] <= self.allowable_error:
            best_indiv = self.sp.gene2val(population[0, :self.n_gene])
            best_fitness = population[0, -1]

            return best_indiv, best_fitness

        generation = 1
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
            best_fitness = population[0, -1]
            np.save(
                self.model_path + '/out/{:d}/best_fitness.npy'.format(
                    nth_paramset
                ), best_fitness
            )
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
