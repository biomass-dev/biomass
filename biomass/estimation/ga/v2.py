import os
import warnings
import time
import numpy as np
from typing import Optional, NoReturn

from ...exec_model import BioMassModel, ExecModel
from .rcga import RealCodedGeneticAlgorithm


class OptimizeWarning(UserWarning):
    pass


def _check_unknown_options(unknown_options: dict) -> Optional[NoReturn]:
    if unknown_options:
        msg = ", ".join(map(str, unknown_options.keys()))
        warnings.warn(f"Unknown solver options: {msg}", OptimizeWarning)


class GeneticAlgorithmInit(ExecModel):
    def __init__(
        self,
        model: BioMassModel,
        popsize: int,
        max_generation: int,
        allowable_error: float,
        local_search_method: str,
        n_children: int,
        maxiter: int,
        workers: int,
        overwrite: bool,
        **unknown_options,
    ) -> None:
        super().__init__(model)
        self.search_rgn: np.ndarray = self.model.sp.get_region()
        self.n_population: int = int(popsize * self.search_rgn.shape[1])
        self.n_gene: int = self.search_rgn.shape[1]
        self.max_generation: int = max_generation
        self.allowable_error: float = allowable_error
        self.local_search_method: str = local_search_method.lower()
        self.n_children: int = n_children
        self.maxiter: int = maxiter
        self.workers: int = workers
        self.overwrite: bool = overwrite

        _check_unknown_options(unknown_options)
        if self.n_population < self.n_gene + 2:
            raise ValueError(f"n_population must be larger than {self.n_gene + 2:d}")

    def run(self, nth_paramset: int) -> None:
        os.makedirs(
            os.path.join(
                self.model.path,
                "out",
            ),
            exist_ok=True,
        )
        if not self.overwrite and os.path.isdir(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
            )
        ):
            raise FileExistsError(
                "Set options['overwrite']=True to overwrite "
                + os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                )
            )
        else:
            try:
                files = os.listdir(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                    )
                )
                for file in files:
                    if any(map(file.__contains__, (".npy", ".log"))):
                        os.remove(
                            os.path.join(
                                self.model.path,
                                "out",
                                f"{nth_paramset:d}",
                                f"{file}",
                            )
                        )
            except FileNotFoundError:
                os.mkdir(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                    )
                )
        np.random.seed(time.time_ns() * nth_paramset % 2 ** 32)
        warnings.filterwarnings("ignore")
        self._ga_v2(nth_paramset)

    def _set_initial(self, nth_paramset: int) -> np.ndarray:
        population = np.full((self.n_population, self.n_gene + 1), np.inf)
        with open(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "optimization.log",
            ),
            mode="w",
        ) as f:
            f.write("Generating the initial population. . .\n")
        for i in range(self.n_population):
            while 1e12 <= population[i, -1] or np.isnan(population[i, -1]):
                population[i, : self.n_gene] = np.random.rand(self.n_gene)
                population[i, -1] = self.model.obj_func(population[i, : self.n_gene])
            with open(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "optimization.log",
                ),
                mode="a",
            ) as f:
                f.write(f"{i + 1:d} / {self.n_population:d}\n")
        population = population[np.argsort(population[:, -1]), :]

        return population

    def _ga_v2(self, nth_paramset: int) -> None:
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

        4. Local Search
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
        rcga = RealCodedGeneticAlgorithm(
            self.model.obj_func,
            self.n_population,
            self.n_gene,
            self.n_children,
            self.maxiter,
            self.workers,
        )
        n_iter = 1
        n0 = np.empty(3 * self.n_population)

        population = self._set_initial(nth_paramset)
        n0[0] = population[0, -1]

        with open(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "optimization.log",
            ),
            mode="a",
        ) as f:
            f.write(
                "\n----------------------------------------\n\n"
                + f"Generation1: Best Fitness = {population[0, -1]:e}\n"
            )
        best_individual = self.model.sp.gene2val(population[0, : self.n_gene])
        best_fitness = population[0, -1]

        np.save(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "generation.npy",
            ),
            1,
        )
        np.save(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "fit_param1.npy",
            ),
            best_individual,
        )
        np.save(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "best_fitness.npy",
            ),
            best_fitness,
        )
        if population[0, -1] <= self.allowable_error:
            return

        generation = 1
        while generation < self.max_generation:
            ip = np.random.choice(self.n_population, self.n_gene + 2, replace=False)
            population = rcga.converging(ip, population)
            population = rcga.local_search(ip, population, self.local_search_method)
            for _ in range(n_iter - 1):
                ip = np.random.choice(self.n_population, self.n_gene + 2, replace=False)
                population = rcga.converging(ip, population)
            # Adaptation of n_iter
            if generation % len(n0) == len(n0) - 1:
                n0[-1] = population[0, -1]
                if n0[0] == n0[-1]:
                    n_iter *= 2
                else:
                    n_iter = 1
            else:
                n0[generation % len(n0)] = population[0, -1]

            best_individual = self.model.sp.gene2val(population[0, : self.n_gene])
            if population[0, -1] < best_fitness:
                np.save(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                        "generation.npy",
                    ),
                    generation + 1,
                )
                np.save(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                        f"fit_param{generation + 1:d}.npy",
                    ),
                    best_individual,
                )
            best_fitness = population[0, -1]
            np.save(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "best_fitness.npy",
                ),
                best_fitness,
            )
            np.save(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "count_num.npy",
                ),
                generation + 1,
            )
            with open(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "optimization.log",
                ),
                mode="a",
            ) as f:
                f.write(f"Generation{generation + 1:d}: " f"Best Fitness = {best_fitness:e}\n")
            if population[0, -1] <= self.allowable_error:
                break

            generation += 1

        return


class GeneticAlgorithmContinue(ExecModel):
    def __init__(
        self,
        model: BioMassModel,
        popsize: int,
        max_generation: int,
        allowable_error: float,
        local_search_method: str,
        n_children: int,
        maxiter: int,
        workers: int,
        p0_bounds: list,
        **unknown_options,
    ) -> None:
        super().__init__(model)
        self.search_rgn: np.ndarray = self.model.sp.get_region()
        self.n_population: int = int(popsize * self.search_rgn.shape[1])
        self.n_gene: int = self.search_rgn.shape[1]
        self.n_children: int = 50
        self.max_generation: int = max_generation
        self.allowable_error: float = allowable_error
        self.local_search_method: str = local_search_method.lower()
        self.n_children: int = n_children
        self.maxiter: int = maxiter
        self.workers: int = workers
        self.p0_bounds: list = p0_bounds

        _check_unknown_options(unknown_options)
        if self.n_population < self.n_gene + 2:
            raise ValueError(f"n_population must be larger than {self.n_gene + 2:d}")

    def run(self, nth_paramset: int) -> None:

        np.random.seed(time.time_ns() * nth_paramset % 2 ** 32)
        warnings.filterwarnings("ignore")
        self._my_ga_continue(nth_paramset)

    def _set_continue(self, nth_paramset: int) -> np.ndarray:
        best_generation = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "generation.npy",
            )
        )
        best_individual = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                f"fit_param{int(best_generation):d}.npy",
            )
        )
        population = np.full((self.n_population, self.n_gene + 1), np.inf)

        with open(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "optimization.log",
            ),
            mode="a",
        ) as f:
            f.write(
                "\n########################################"
                "\n############### Continue ###############"
                "\n########################################"
                "\nGenerating the initial population. . .\n"
            )
        for i in range(self.n_population):
            while 1e12 <= population[i, -1]:
                population[i, : self.n_gene] = self._encode_bestIndivVal2randGene(best_individual)
                population[i, : self.n_gene] = np.clip(population[i, : self.n_gene], 0.0, 1.0)
                population[i, -1] = self.model.obj_func(population[i, : self.n_gene])
            with open(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "optimization.log",
                ),
                mode="a",
            ) as f:
                f.write(f"{i + 1:d} / {self.n_population:d}\n")
        population = population[np.argsort(population[:, -1]), :]

        return population

    def _encode_bestIndivVal2randGene(self, best_individual: np.ndarray) -> np.ndarray:
        rand_gene = (
            np.log10(
                best_individual
                * 10
                ** (
                    np.random.rand(len(best_individual)) * np.log10(self.p0_bounds[1] / self.p0_bounds[0])
                    + np.log10(self.p0_bounds[0])
                )
            )
            - self.search_rgn[0, :]
        ) / (self.search_rgn[1, :] - self.search_rgn[0, :])

        return rand_gene

    def _my_ga_continue(self, nth_paramset: int) -> None:
        rcga = RealCodedGeneticAlgorithm(
            self.model.obj_func,
            self.n_population,
            self.n_gene,
            self.n_children,
            self.maxiter,
            self.workers,
        )
        n_iter = 1
        n0 = np.empty(3 * self.n_population)

        count_num = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "count_num.npy",
            )
        )
        best_generation = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "generation.npy",
            )
        )
        best_individual = np.load(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                f"fit_param{int(best_generation):d}.npy",
            )
        )
        best_individual_gene = self.model.sp.val2gene(best_individual)
        best_fitness = self.model.obj_func(best_individual_gene)

        if self.max_generation <= count_num:
            raise ValueError(f"max_generation should be larger than {int(count_num):d}")

        population = self._set_continue(nth_paramset)
        if best_fitness < population[0, -1]:
            population[0, : self.n_gene] = best_individual_gene
            population[0, -1] = best_fitness
        else:
            best_individual = self.model.sp.gene2val(population[0, : self.n_gene])
            best_fitness = population[0, -1]
            np.save(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    f"fit_param{int(count_num) + 1:d}.npy",
                ),
                best_individual,
            )
        with open(
            os.path.join(
                self.model.path,
                "out",
                f"{nth_paramset:d}",
                "optimization.log",
            ),
            mode="a",
        ) as f:
            f.write(
                "\n----------------------------------------\n\n"
                f"Generation{int(count_num) + 1:d}: "
                f"Best Fitness = {best_fitness:e}\n"
            )
        n0[0] = population[0, -1]

        if population[0, -1] <= self.allowable_error:
            return

        generation = 1 + int(count_num)
        while generation < self.max_generation:
            ip = np.random.choice(self.n_population, self.n_gene + 2, replace=False)
            population = rcga.converging(ip, population)
            population = rcga.local_search(ip, population, self.local_search_method)
            for _ in range(n_iter - 1):
                ip = np.random.choice(self.n_population, self.n_gene + 2, replace=False)
                population = rcga.converging(ip, population)
            if generation % len(n0) == len(n0) - 1:
                n0[-1] = population[0, -1]
                if n0[0] == n0[-1]:
                    n_iter *= 2
                else:
                    n_iter = 1
            else:
                n0[generation % len(n0)] = population[0, -1]

            best_individual = self.model.sp.gene2val(population[0, : self.n_gene])
            if population[0, -1] < best_fitness:
                np.save(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                        "generation.npy",
                    ),
                    generation + 1,
                )
                np.save(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                        f"fit_param{generation + 1:d}.npy",
                    ),
                    best_individual,
                )
                np.save(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                        "best_fitness",
                    ),
                    best_fitness,
                )
            best_fitness = population[0, -1]
            np.save(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "count_num.npy",
                ),
                generation + 1,
            )
            with open(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "optimization.log",
                ),
                mode="a",
            ) as f:
                f.write(f"Generation{generation + 1:d}: " f"Best Fitness = {best_fitness:e}\n")
            if population[0, -1] <= self.allowable_error:
                break

            generation += 1

        return