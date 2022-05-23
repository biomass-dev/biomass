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
    """

    model: ModelObject

    def generate(
        self,
        popsize: int = 3,
        threshold: float = 1e12,
        progress: bool = False,
    ) -> np.ndarray:
        """
        Return initial population for ``optimizer_options['init'] `` in ``biomass.optimize`` function.

        Parameters
        ----------
        popsize : int (default: 3)
            A multiplier for setting the total population size.
        threshold : float (default: 1e12)
            Allowable error for generating initial population.
            Default value is 1e12 (numerically solvable).
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
        >>> initpop = InitialPopulation(model).generate(progress=True)
        >>> optimize(model, x_id=1, optimizer_options={"init": initpop})
        """
        n_gene = len(self.model.problem.bounds)
        n_population = popsize * n_gene
        population = np.full((n_population, n_gene + 1), np.inf)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in tqdm(range(n_population), disable=not progress):
                while threshold <= population[i, -1] or not isfinite(population[i, -1]):
                    population[i, : n_gene] = np.random.rand(n_gene)
                    population[i, -1] = self.model.get_obj_val(population[i, : n_gene])
        population = population[np.argsort(population[:, -1]), :]
        return population[:, : n_gene]