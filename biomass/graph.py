import os
import re
import warnings
from collections import defaultdict
from types import ModuleType
from typing import List, Literal

import pygraphviz as pgv


class NetworkGraph(object):
    """
    Visualization of the biological network as a graph.

    Parameters
    ----------
    path : str
        Path to a biomass model.
    biomass_model : Any
        A package containing biomass model properties.

    Attributes
    ----------
    path : str
        Path to the model.
    parameters : list of strings
        Names of model parameters.
    species : list of strings
        Names of model species.
    observables : list of strings
        Names of model observables.
    pval : Callable
        Numerical values of the parameters.
    ival : Callable
        Initial values.
    problem : OptimizationProblem
        The optimization problem.
    viz : Visualization
        Plotting parameters for customizing figure properties.
    rxn : ReactionNetwork
        Reaction indices grouped according to biological processes.
    """

    def __init__(self, path: str, biomass_model: ModuleType):
        self._path = path
        self._parameters = biomass_model.C.NAMES
        self._species = biomass_model.V.NAMES
        self.pval = biomass_model.param_values
        self.ival = biomass_model.initial_values
        self.problem = biomass_model.OptimizationProblem()
        self.viz = biomass_model.Visualization()
        self.rxn = biomass_model.ReactionNetwork()

    @property
    def path(self) -> str:
        return self._path

    @property
    def parameters(self) -> List[str]:
        return self._parameters

    @property
    def species(self) -> list:
        return self._species

    @property
    def observables(self) -> List[str]:
        duplicate = [
            name for name in set(self.problem.obs_names) if self.problem.obs_names.count(name) > 1
        ]
        if not duplicate:
            return self.problem.obs_names
        else:
            raise NameError(f"Duplicate observables: {', '.join(duplicate)}")

    @staticmethod
    def _extract_equation(dir: str, filepath: str, left_match: str, right_match: str) -> dict:
        """
        Matches species that interact with each other by reading python files as text and looking for equations.
        """
        data = {}
        with open(os.path.join(dir, filepath)) as f:
            was_warned = False
            lefthand = None
            equation_cache = ""
            for line in f.readlines():
                match = re.search(left_match, line)
                cache_full = bool(len(equation_cache))
                if "return" in line:
                    if cache_full:
                        equation_cache = "".join(equation_cache).split("=")[1]
                        data[lefthand] = [
                            found
                            for found in re.findall(right_match, equation_cache)
                            if found != lefthand
                        ]
                    equation_cache = ""
                    cache_full = False
                if match is not None:
                    if cache_full:
                        equation_cache = "".join(equation_cache).split("=")[1]
                        data[lefthand] = [
                            found
                            for found in re.findall(right_match, equation_cache)
                            if found != lefthand
                        ]
                        equation_cache = ""
                        cache_full = False
                    lefthand = match.group(0)
                    equation_cache += line
                elif cache_full:
                    equation_cache += line
                if len(equation_cache) == 0:
                    try:
                        warning_match = re.findall(r"\[V\.", line.split("=")[1])
                    except IndexError:
                        warning_match = None
                    if warning_match and not was_warned:
                        warnings.warn(
                            "Usage of species concentrations outside of rate equations detected."
                            + " Species concentrations outside of rate equations will not be considered in Graph generation."
                            + " Might lead to faulty graph."
                        )
                        was_warned = True
        return data

    def to_graph(
        self,
        file_name: str,
        gviz_args: str = "",
        gviz_prog: Literal["neato", "dot", "twopi", "circo", "fdp", "nop"] = "dot",
    ) -> None:
        """Constructs and draws a directed graph of the model.
        Using the pygraphviz library and graphviz a directed graph of the model is constructed by parsing the equations from
        ode.py/reaction_network.py. Equations will be split at the equal sign and an edge is added between the species on the
        lefthand side to all species on the righthand side. Self references will not be considered. Graphviz is then used to
        save the graph as an image.

        Parameters
        ----------
        file_name : string
                    Name as which the image of the graph will be stored.
        gviz_args : string, optional, default=""
                    Used to specify command line options for gviz, see https://graphviz.org/pdf/dot.1.pdf for available options.
        gviz_prog : {"neato", "dot", "twopi", "circo", "fdp", "nop"}, default="dot"
                    Layout engine with which the graph will be arranged. For details see https://graphviz.org/docs/layouts/ .

        Raises
        ------
        AssertionError
            If more species are part of the model, than are detected from ode.py/reaction_network.py or vice versa.

        Warns
        -----
        UserWarning
            If species equations are detected outside of the ODE section.

        Examples
        --------
        >>> model.create_graph("path/to/graph.png")
        Creates graph with dot layout and default options.
        >>> model.create_graph("path/to/graph.pdf, gviz_prog="-Nshape=box -Nstyle=filled -Nfillcolor="#ffe4c4" -Edir=none")
        Creates graph with dot layout in pdf file format. Nodes will be rectangular and colored bisque, edges will have no arrows indicating direction.
        """

        # Check gviz_prog
        if gviz_prog not in (available_layout := ["neato", "dot", "twopi", "circo", "fdp", "nop"]):
            raise ValueError(
                f"gviz_prog must be one of [{', '.join(available_layout)}], got {gviz_prog}."
            )
        try:
            if len(self.rxn.flux(0, self.ival(), self.pval())) > 0:
                use_flux = True
            else:
                use_flux = False
        except AttributeError:
            use_flux = False

        if use_flux is False:
            left_re = r"(?<=dydt\[V.)(.+?)(?=\])"
            right_re = r"(?<=y\[V.)(.+?)\]"
            edges = self._extract_equation(self.path, "ode.py", left_re, right_re)
        elif use_flux is True:
            left_re_flux = r"(?<=v\[)(.+?)(?=\])"
            right_re_flux = r"(?<=y\[V.)(.+?)\]"
            left_re_ode = r"(?<=dydt\[V.)(.+?)(?=\])"
            right_re_ode = r"(?<=v\[)(.+?)\]"
            fluxes = self._extract_equation(
                self.path, "reaction_network.py", left_re_flux, right_re_flux
            )
            odes = self._extract_equation(self.path, "ode.py", left_re_ode, right_re_ode)
            edges = defaultdict(lambda: [])
            for species, velocities in odes.items():
                for velocity in velocities:
                    for participant in fluxes[velocity]:
                        if participant != species:
                            edges[species].append(participant)
        num_participants = []
        for participants in edges.values():
            num_participants += participants
        num_participants += edges.keys()
        assert set(self.species) == set(
            num_participants
        ), f"Not all species are extracted. Expected {len(self.species)}, got {len(set(num_participants))}"
        graph = pgv.AGraph(directed=True)
        for species, partners in edges.items():
            graph.add_node(species)
            for partner in partners:
                graph.add_edge(partner, species)
        self.graph = graph
        graph.layout(prog=gviz_prog, args=gviz_args)
        graph.draw(os.path.join(self.path, file_name))
