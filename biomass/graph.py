import os
import re
import warnings
from collections import defaultdict
from types import ModuleType
from typing import List, Literal, Optional


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

        self.graph = None

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

    def to_graph(self) -> None:
        """Constructs a directed graph of the model.
        Using the pygraphviz library a directed graph of the model is constructed by parsing the equations from
        ode.py/reaction_network.py. Equations will be split at the equal sign and an edge is added between the species on the
        lefthand side to all species on the right hand side. Self references will not be considered.

        Raises
        ------
        AssertionError
            If more species are part of the model, than are detected from ode.py/reaction_network.py or vice versa.

        Warns
        -----
        UserWarning
            If species equations are detected outside of the ODE section.
        """
        try:
            import pygraphviz as pgv
        except ImportError:
            print("pygraphviz is required to run this function.")

        use_flux = False
        try:
            if len(self.rxn.flux(0, self.ival(), self.pval())) > 0:
                use_flux = True
            else:
                use_flux = False
        except AttributeError:
            use_flux = False

        if not use_flux:
            left_re = r"(?<=dydt\[V.)(.+?)(?=\])"
            right_re = r"(?<=y\[V.)(.+?)\]"
            edges = self._extract_equation(self.path, "ode.py", left_re, right_re)
        else:
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

    def static_plot(
        self,
        save_dir: str = "",
        file_name: str = "model_graph.png",
        gviz_args: str = "",
        gviz_prog: Literal["neato", "dot", "twopi", "circo", "fdp", "nop"] = "dot",
    ) -> None:
        """Saves a static image of the network.

        Static image is created using pygraphviz.

        Parameters
        ----------
        save_dir : string
            Name of the directory in which the image will be stored. If empty will be the path of the model folder.
        file_name : string
            Name as which the image of the graph will be stored.
        gviz_args : string, optional, default=""
            Used to specify command line options for gviz, see https://graphviz.org/pdf/dot.1.pdf for available options.
        gviz_prog : {"neato", "dot", "twopi", "circo", "fdp", "nop"}, default="dot"
            Layout engine with which the graph will be arranged. For details see https://graphviz.org/docs/layouts/ .

        Raises
        ------
        ValueError
            If something is passed as the gviz_prog that is not a viable layout program.

        Examples
        --------
        >>> model.static_plot("path/to/", "graph.png")
        Creates graph with dot layout and default options.
        >>> model.static_plot("path/to/", "graph.pdf", gviz_prog="-Nshape=box -Nstyle=filled -Nfillcolor="#ffe4c4" -Edir=none")
        Creates graph with dot layout in pdf file format. Nodes will be rectangular and colored bisque, edges will have no arrows indicating direction.

        """
        if self.graph is None:
            self.to_graph()
        if gviz_prog not in (available_layout := ["neato", "dot", "twopi", "circo", "fdp", "nop"]):
            raise ValueError(
                f"gviz_prog must be one of [{', '.join(available_layout)}], got {gviz_prog}."
            )
        if not save_dir:
            save_dir = self.path
        self.graph.layout(prog=gviz_prog, args=gviz_args)
        self.graph.draw(os.path.join(save_dir, file_name))

    def dynamic_plot(
        self,
        save_dir: str = "",
        file_name: str = "network.html",
        show: bool = True,
        annotate_nodes: bool = True,
        show_controls: bool = False,
        which_controls: Optional[List[str]] = None,
    ) -> None:
        """Saves a dynamic and interactive image of the network graph.
        Graph is read by pyvis. Using pyvis a dynamic and interactive representation of the biological network
        is created in html format.

        Parameters
        ----------
        save_dir : string
            Name of the directory in which the image will be stored. If empty will be the path of the model folder.
        file_name : string
            Name as which the image of the graph will be stored.
        show: bool, default=True
            If true the plot will immediately be displayed in the webbrowser.
        annotate_nodes : bool, default=True
            If true nodes will be scaled according to number of edges and hovering over a node will show interaction partners.
        show_controls : bool, default=False
            If true control buttons will be displayed.
        which_controls : List(str), optional, default=None
            Used to specify which control buttons should be displayed. If empty all buttons will be displayed.

        Examples
        --------
        >>> model.dynamic_plot("path/to/", "graph.html")
        Creates graph and shows interactive graph with default options.
        >>> model.dynamic_plot("path/to/", "graph.html", show=False, show_controls=True, which_controls=["physics", "manipulation", "interaction"])
        Creates interactive graph. Controls for physics, manipulation and interaction will be available.
        """
        try:
            from pyvis.network import Network
        except ImportError:
            print("pyvis is required to run this function.")

        if self.graph is None:
            self.to_graph()
        if os.path.splitext(file_name)[1] != ".html":
            file_name = file_name + ".html"
        if not save_dir:
            save_dir = self.path
        network = Network(directed=True)
        network.add_nodes(self.graph.nodes())
        network.add_edges(self.graph.edges())
        if not isinstance(show_controls, bool):
            raise TypeError(f"show_controls is type {type(show_controls)}, needs to be boolean")
        if show_controls:
            network.show_buttons(filter_=which_controls)
        if annotate_nodes:
            neighbor_map = network.get_adj_list()
            for node in network.nodes:
                node["title"] = " Neighbors:\n"
                node["title"] += "\n".join(neighbor_map[node["id"]])
                node["value"] = len(neighbor_map[node["id"]])
        if show:
            network.show(os.path.join(save_dir, file_name))
        else:
            network.save_graph(os.path.join(save_dir, file_name))
