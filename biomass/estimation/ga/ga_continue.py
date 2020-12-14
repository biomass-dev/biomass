import time
import numpy as np

from biomass.exec_model import ExecModel
from .rcga import RealCodedGeneticAlgorithm


class GeneticAlgorithmContinue(ExecModel):
    def __init__(
        self,
        model,
        popsize,
        max_generation,
        allowable_error,
        local_search_method,
        n_children,
        workers,
        p0_bounds,
    ):
        super().__init__(model)
        self.search_rgn: np.ndarray = self.sp.get_region()
        self.n_population: int = int(popsize * self.search_rgn.shape[1])
        self.n_gene: int = self.search_rgn.shape[1]
        self.n_children: int = 50
        self.max_generation: int = max_generation
        self.allowable_error: float = allowable_error
        self.local_search_method: str = local_search_method.lower()
        self.n_children: int = n_children
        self.workers: int = workers
        self.p0_bounds: list = p0_bounds

        if self.n_population < self.n_gene + 2:
            raise ValueError(f"self.n_population must be larger than {self.n_gene + 2:d}")

        if self.local_search_method not in ["mutation", "powell", "de"]:
            raise ValueError(
                f"'{local_search_method}': Invalid local_search_method. "
                "Should be one of ['mutation', 'Powell', 'DE']"
            )

    def run(self, nth_paramset: int) -> None:

        np.random.seed(time.time_ns() * nth_paramset % 2 ** 32)
        self._my_ga_continue(nth_paramset)

    def _set_continue(self, nth_paramset: int) -> np.ndarray:
        best_generation = np.load(self.model_path + f"/out/{nth_paramset:d}/generation.npy")
        best_individual = np.load(self.model_path + f"/out/{nth_paramset:d}/fit_param{int(best_generation):d}.npy")
        population = np.full((self.n_population, self.n_gene + 1), np.inf)

        with open(self.model_path + f"/out/{nth_paramset:d}/optimization.log", mode="a") as f:
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
                population[i, -1] = self.obj_func(population[i, : self.n_gene])
            with open(self.model_path + f"/out/{nth_paramset:d}/optimization.log", mode="a") as f:
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
        rcga = RealCodedGeneticAlgorithm(self.obj_func, self.n_population, self.n_gene, self.n_children, self.workers)
        n_iter = 1
        n0 = np.empty(3 * self.n_population)

        count_num = np.load(self.model_path + f"/out/{nth_paramset:d}/count_num.npy")
        best_generation = np.load(self.model_path + f"/out/{nth_paramset:d}/generation.npy")
        best_individual = np.load(self.model_path + f"/out/{nth_paramset:d}/fit_param{int(best_generation):d}.npy")
        best_individual_gene = self.sp.val2gene(best_individual)
        best_fitness = self.obj_func(best_individual_gene)

        if self.max_generation <= count_num:
            raise ValueError(f"max_generation should be larger than {int(count_num):d}")

        population = self._set_continue(nth_paramset)
        if best_fitness < population[0, -1]:
            population[0, : self.n_gene] = best_individual_gene
            population[0, -1] = best_fitness
        else:
            best_individual = self.sp.gene2val(population[0, : self.n_gene])
            best_fitness = population[0, -1]
            np.save(
                self.model_path + f"/out/{nth_paramset:d}/fit_param{int(count_num) + 1:d}.npy",
                best_individual,
            )
        with open(self.model_path + f"/out/{nth_paramset:d}/optimization.log", mode="a") as f:
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

            best_individual = self.sp.gene2val(population[0, : self.n_gene])
            if population[0, -1] < best_fitness:
                np.save(
                    self.model_path + f"/out/{nth_paramset:d}/generation.npy",
                    generation + 1,
                )
                np.save(
                    self.model_path + f"/out/{nth_paramset:d}/fit_param{generation + 1:d}.npy",
                    best_individual,
                )
                np.save(
                    self.model_path + f"/out/{nth_paramset:d}/best_fitness",
                    best_fitness,
                )
            best_fitness = population[0, -1]
            np.save(self.model_path + f"/out/{nth_paramset:d}/count_num.npy", generation + 1)
            with open(self.model_path + f"/out/{nth_paramset:d}/optimization.log", mode="a") as f:
                f.write(f"Generation{generation + 1:d}: " f"Best Fitness = {best_fitness:e}\n")
            if population[0, -1] <= self.allowable_error:
                break

            generation += 1

        return