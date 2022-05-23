import multiprocessing
import os
import shutil
import warnings
from dataclasses import dataclass
from math import isfinite

import numpy as np
from tqdm import tqdm

from ..exec_model import ModelObject


@dataclass
class InitialPopulation(object):
    """
    Attributes
    ----------
    model : ModelObject
        BioMASS model object.
    popsize : int (default: 3)
        A multiplier for setting the total population size.
    threshold : float (default: 1e12)
        Allowable error for generating initial population.
        Default value is 1e12 (numerically solvable).
    """

    model: ModelObject
    popsize: int = 3
    threshold: float = 1e12

    def __post_init__(self) -> None:
        self.n_gene = len(self.model.problem.bounds)
        self.n_population = self.popsize * self.n_gene
        self.initpop_path = os.path.join(self.model.path, "_initpop")
        os.makedirs(self.initpop_path, exist_ok=True)

    def _set_gene_vector(self, pop_id: int) -> None:
        obj_val = np.inf
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            while self.threshold <= obj_val or not isfinite(obj_val):
                gene = np.random.rand(self.n_gene)
                obj_val = self.model.get_obj_val(gene)
            np.save(os.path.join(self.initpop_path, f"population_{pop_id}"), gene)

    def generate(self, n_proc: int = 1, progress: bool = False) -> np.ndarray:
        """
        Return initial population for ``optimizer_options['init'] `` in ``biomass.optimize`` function.

        Parameters
        ----------
        n_proc : int (default: 1)
            Number of processes to use (default: 1). Set to a larger number
            (e.g. the number of CPU cores available) for parallel execution of simulations.
        progress : bool (default: :obj:`False`)
            Whether the progress bar is animating or not.

        Returns
        -------
        population : numpy.ndarray
            Array specifying the initial population.
            The array should have shape (M, len(x)),
            where M is the total population size and len(x) is the number of parameters.

        Examples
        --------
        >>> from biomass import create_model, optimize
        >>> from biomass.estimation import InitialPopulation
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> model = create_model(Nakakuki_Cell_2010.__package__)
        >>> initpop = InitialPopulation(model).generate(n_proc=2, progress=True)
        >>> optimize(model, x_id=1, optimizer_options={"init": initpop})
        """
        p = multiprocessing.Pool(processes=n_proc)
        population = np.empty((self.n_population, self.n_gene))
        try:
            with tqdm(total=self.n_population, disable=not progress) as t:
                for _ in p.imap_unordered(self._set_gene_vector, range(self.n_population)):
                    t.update(1)
            for i in range(self.n_population):
                gene = np.load(os.path.join(self.initpop_path, f"population_{i}.npy"))
                population[i] = gene
            return population
        finally:
            p.close()
            shutil.rmtree(self.initpop_path)