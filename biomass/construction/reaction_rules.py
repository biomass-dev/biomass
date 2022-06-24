import re
import sys
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from typing import Dict, List, NamedTuple, NoReturn, Optional

import numpy as np

from .thermodynamic_restrictions import ComplexFormation, DuplicateError, ThermodynamicRestrictions


class UnregisteredRule(NamedTuple):
    expected: Optional[str]
    original: Optional[str]


class DetectionError(Exception):
    pass


class ArrowError(ValueError):
    pass


PREPOSITIONS: List[str] = [
    "to",
    "for",
    "from",
    "up",
    "down",
    "in",
    "on",
    "at",
    "off",
    "into",
    "around",
    "among",
    "between",
    "of",
    "over",
    "above",
    "below",
    "under",
    "through",
    "across",
    "along",
    "near",
    "by",
    "beside",
]


@dataclass
class ReactionRules(ThermodynamicRestrictions):
    """Create an executable biochemical model from text.

    .. list-table:: Available reaction rules
        :widths: 25 50 25
        :header-rows: 1

        * - Rule
          - Example sentence
          - Parameters (optional)
        * - :func:`~pasmopy.construction.reaction_rules.dimerize`
          - *A* dimerizes <--> *AA*
          - .. math:: kf, kr
        * - :func:`~pasmopy.construction.reaction_rules.bind`
          - *A* binds *B* <--> *AB*
          - .. math:: kf, kr
        * - :func:`~pasmopy.construction.reaction_rules.dissociate`
          - *AB* dissociates to *A* and *B*
          - .. math:: kf, kr
        * - :func:`~pasmopy.construction.reaction_rules.is_phosphorylated`
          - *uA* is phosphorylated <--> *pA*
          - .. math:: kf, kr
        * - :func:`~pasmopy.construction.reaction_rules.is_dephosphorylated`
          - *pA* is dephosphorylated --> *uA*
          - .. math:: V, K
        * - :func:`~pasmopy.construction.reaction_rules.phosphorylate`
          - *B* phosphorylates *uA* --> *pA*
          - .. math:: V, K
        * - :func:`~pasmopy.construction.reaction_rules.dephosphorylate`
          - *B* dephosphorylates *pA* --> *uA*
          - .. math:: V, K
        * - :func:`~pasmopy.construction.reaction_rules.transcribe`
          - *B* transcribes *a*
          - .. math:: V, K, n, (KF, nF)
        * - :func:`~pasmopy.construction.reaction_rules.synthesize`
          - *B* synthesizes *A*
          - .. math:: kf
        * - :func:`~pasmopy.construction.reaction_rules.is_synthesized`
          - *A* is synthesized
          - .. math:: kf
        * - :func:`~pasmopy.construction.reaction_rules.degrade`
          - *B* degrades *A*
          - .. math:: kf
        * - :func:`~pasmopy.construction.reaction_rules.is_degraded`
          - *A* is degraded
          - .. math:: kf
        * - :func:`~pasmopy.construction.reaction_rules.translocate`
          - *Acyt* translocates from cytoplasm to nucleus (Vcyt, Vnuc) <--> *Anuc*
          - .. math:: kf, kr, (V_{pre}, V_{post})
        * - :func:`~pasmopy.construction.reaction_rules.user_defined`
          - @rxn Reactant --> Product: *define rate equation here*
          - *user-defined*

    From v0.2.2, you can specify directionality in binding-dissociation reaction via different arrows:

    .. code-block:: python

        E + S ⇄ ES | kf=0.003, kr=0.001 | E=100, S=50  # bi-directional
        ES → E + P | kf=0.002  # unidirectional

    Attributes
    ----------
    input_txt : str
        Model description file (*.txt), e.g.,
        `Kholodenko1999.txt <https://github.com/pasmopy/pasmopy/blob/master/tests/text_files/Kholodenko1999.txt>`_.
    similarity_threshold : float
        If all match_scores are below this value, expected_word will not be returned.
    parameters : list of strings
        ``x`` : model parameters.
    species : list of strings
        ``y`` : model species.
    reactions : list of strings
        ``v`` : flux vector.
    differential_equations : list of strings
        ``dydt`` : right-hand side of the differential equation.
    obs_desc : list of List[str]
        Description of observables.
    param_info : list of strings
        Information about parameter values.
    init_info : list of strings
        Information about initial values.
    param_constraints : list of strings
        Information about parameter constraints.
    param_excluded : list of strings
        List of parameters excluded from search params because of parameter constraints.
    fixed_species : list of strings
        List of species which should be held fixed (never consumed) during simulation.
    sim_tspan : list of strings ['t0', 'tf']
        Interval of integration.
    sim_conditions : list of List[str]
        Simulation conditions with stimulation.
    sim_unperturbed : str
        Untreated conditions to get steady state.
    rule_words : dict
        Words to identify reaction rules.
    nothing : List[str]
        Available symbol for degradation/creation to/from nothing.
    fwd_arrows : List[str]
        Available arrows for unidirectional reactions.
    double_arrows : List[str]
        Available arrows for bi-directional reactions.

    """

    input_txt: str
    similarity_threshold: float

    parameters: List[str] = field(
        default_factory=list,
        init=False,
    )
    species: List[str] = field(
        default_factory=list,
        init=False,
    )
    reactions: List[str] = field(
        default_factory=list,
        init=False,
    )
    differential_equations: List[str] = field(
        default_factory=list,
        init=False,
    )
    obs_desc: List[List[str]] = field(
        default_factory=list,
        init=False,
    )
    param_info: List[str] = field(
        default_factory=list,
        init=False,
    )
    init_info: List[str] = field(
        default_factory=list,
        init=False,
    )
    param_constraints: List[str] = field(
        default_factory=list,
        init=False,
    )
    param_excluded: List[str] = field(
        default_factory=list,
        init=False,
    )
    fixed_species: List[str] = field(
        default_factory=list,
        init=False,
    )
    # Information about simulation
    sim_tspan: List[str] = field(
        default_factory=list,
        init=False,
    )
    sim_conditions: List[List[str]] = field(
        default_factory=list,
        init=False,
    )
    sim_unperturbed: str = field(
        default_factory=str,
        init=False,
    )
    # Words to identify reaction rules
    rule_words: Dict[str, List[str]] = field(
        default_factory=lambda: dict(
            _bind_and_dissociate=[
                " +",
            ],
            dimerize=[
                " dimerizes",
                " homodimerizes",
                " forms a dimer",
                " forms dimers",
            ],
            bind=[
                " binds",
                " forms complexes with",
            ],
            dissociate=[
                " is dissociated into",
                " dissociates to",
            ],
            is_phosphorylated=[
                " is phosphorylated",
            ],
            is_dephosphorylated=[
                " is dephosphorylated",
            ],
            phosphorylate=[
                " phosphorylates",
            ],
            dephosphorylate=[
                " dephosphorylates",
            ],
            transcribe=[
                " transcribe",
                " transcribes",
            ],
            synthesize=[
                " synthesizes",
                " promotes synthesis of",
                " increases",
                " is translated into",
            ],
            is_synthesized=[
                " is synthesized",
            ],
            degrade=[
                " degrades",
                " promotes degradation of",
                " decreases",
            ],
            is_degraded=[
                " is degraded",
            ],
            translocate=[
                " translocates",
                " is translocated",
            ],
        ),
        init=False,
    )
    nothing: List[str] = field(
        default_factory=lambda: ["∅", "0"],
        init=False,
    )
    fwd_arrows: List[str] = field(
        default_factory=lambda: [
            " → ",
            " ↣ ",
            " ↦ ",
            " ⇾ ",
            " ⟶ ",
            " ⟼ ",
            " ⥟ ",
            " ⥟ ",
            " ⇀ ",
            " ⇁ ",
            " ⇒ ",
            " ⟾ ",
        ]
        + [" {}> ".format("-" * i) for i in range(3)],
        init=False,
    )
    double_arrows: List[str] = field(
        default_factory=lambda: [
            " ↔ ",
            " ⟷ ",
            " ⇄ ",
            " ⇆ ",
            " ⇌ ",
            " ⇋ ",
            " ⇔ ",
            " ⟺ ",
        ]
        + [" <{}> ".format("-" * i) for i in range(3)],
        init=False,
    )

    def __post_init__(self) -> None:
        if not 0.0 < self.similarity_threshold < 1.0:
            raise ValueError("similarity_threshold must lie within (0, 1).")

    @staticmethod
    def _isfloat(string: str) -> bool:
        """
        Checking if a string can be converted to float.
        """
        try:
            float(string)
            return True
        except ValueError:
            return False

    @staticmethod
    def _remove_prefix(text: str, prefix: str) -> str:
        """
        Remove prefix from a text.
        """
        if text.startswith(prefix):
            return text[len(prefix) :]
        assert False

    def _available_arrows(self) -> List[str]:
        """
        Return all available arrow types.
        """
        return self.fwd_arrows + self.double_arrows

    def _set_params(self, line_num: Optional[int], func_name: Optional[str], *args: str) -> None:
        """
        Set model parameters.
        """
        for p_name in args:
            p_name_used = p_name + (
                f"{line_num:d}"
                if not p_name[-1].isdecimal()
                and isinstance(line_num, int)
                and func_name != "user_defined"
                else ""
            )
            if p_name_used not in self.parameters:
                self.parameters.append(p_name_used)

    def _set_species(self, *args: str) -> None:
        """
        Set model species.
        """
        for s_name in args:
            if s_name not in self.species:
                self.species.append(s_name)

    def _raise_detection_error(self, line_num: int, line: str) -> NoReturn:
        """
        Raise `DetectionError` when a keyword is invalid.
        """
        unregistered_rule = self._get_partial_similarity(line)
        raise DetectionError(
            f"Unregistered words in line{line_num:d}: {line}"
            + (
                f"\nMaybe: '{unregistered_rule.expected.lstrip()}'."
                if unregistered_rule.expected is not None
                else ""
            )
        )

    def _process_pval_section(self, func_name: str, line_num: int, line: str, *args: str) -> None:

        param_values = line.split("|")[1].strip().split(",")
        if all("=" in pval for pval in param_values):
            for pval in param_values:
                base_param = pval.split("=")[0].strip(" ")
                if base_param.startswith("const "):
                    # Parameter names with 'const' will be added to param_excluded.
                    base_param = base_param.split("const ")[-1]
                    fixed = True
                else:
                    fixed = False
                if base_param in args:
                    if self._isfloat(pval.split("=")[1].strip(" ")):
                        self.param_info.append(
                            "x[C."
                            + base_param
                            + (f"{line_num:d}]" if func_name != "user_defined" else "]")
                            + " = "
                            + pval.split("=")[1].strip(" ")
                        )
                        # If a parameter value is initialized to 0.0 or fixed,
                        # then add it to param_excluded.
                        if float(pval.split("=")[1].strip(" ")) == 0.0 or fixed:
                            self.param_excluded.append(
                                base_param
                                + (f"{line_num:d}" if func_name != "user_defined" else "")
                            )
                    else:
                        raise ValueError(
                            f"line{line_num:d}: Parameter value must be int or float."
                        )
                else:
                    raise ValueError(
                        f"line{line_num:d}: '{pval.split('=')[0].strip(' ')}'\n"
                        f"Available parameters are: {', '.join(args)}."
                    )
        elif param_values[0].strip(" ").isdecimal():
            # Parameter constraints
            for param_name in args:
                # base_pname = self._get_base_pname(param_name)
                if (
                    f"{self._get_base_pname(param_name)}{int(param_values[0]):d}"
                ) not in self.parameters:
                    raise ValueError(
                        f"Line {line_num:d} and {int(param_values[0]):d} : "
                        "Different reaction rules in parameter constraints."
                    )
                else:
                    if f"x[C.{param_name}" + (
                        f"{line_num:d}]" if func_name != "user_defined" else "]"
                    ) != (
                        f"x[C.{self._get_base_pname(param_name)}" + f"{int(param_values[0]):d}]"
                    ):
                        self.param_excluded.append(
                            f"{param_name}"
                            + (f"{line_num:d}" if func_name != "user_defined" else "")
                        )
                        self.param_info.append(
                            f"x[C.{param_name}"
                            + (f"{line_num:d}]" if func_name != "user_defined" else "]")
                            + " = "
                            + f"x[C.{self._get_base_pname(param_name)}"
                            + f"{int(param_values[0]):d}]"
                        )
                        self.param_constraints.append(
                            f"x[C.{param_name}"
                            + (f"{line_num:d}]" if func_name != "user_defined" else "]")
                            + " = "
                            + f"x[C.{self._get_base_pname(param_name)}"
                            + f"{int(param_values[0]):d}]"
                        )
        else:
            raise ValueError(
                f"line{line_num:d}: {line}\nInvalid expression in the input parameter."
            )

    def _process_ival_section(self, line_num: int, line: str) -> None:

        initial_values = line.split("|")[2].strip().split(",")
        for ival in initial_values:
            if ival.startswith("fixed "):
                ival = ival.split("fixed ")[-1]
                self.fixed_species.append(ival.split("=")[0].strip(" "))
            if ival.split("=")[0].strip(" ") in line.split("|")[0]:
                if self._isfloat(ival.split("=")[1].strip(" ")):
                    self.init_info.append(
                        "y0[V."
                        + ival.split("=")[0].strip(" ")
                        + "] = "
                        + ival.split("=")[1].strip(" ")
                    )
                else:
                    raise ValueError(f"line{line_num:d}: Initial value must be int or float.")
            else:
                raise NameError(
                    f"line{line_num:d}: " f"Name'{ival.split('=')[0].strip(' ')}' is not defined."
                )

    @staticmethod
    def _get_base_pname(param_name: str) -> str:
        while param_name[-1].isdecimal():
            param_name = param_name[:-1]
        return param_name

    def _preprocessing(self, func_name: str, line_num: int, line: str, *args: str) -> List[str]:
        """
        Extract the information about parameter and/or initial values
        if '|' in the line and find a keyword to identify reaction rules.

        Parameters
        ----------
        func_name : str
            Name of the rule function.

        line_num : int
            Line number.

        line : str
            Each line of the input text.

        Returns
        -------
        description : list of strings

        """
        self._set_params(line_num, func_name, *args)
        if "|" in line:
            if line.split("|")[1].strip():
                self._process_pval_section(func_name, line_num, line, *args)
            if line.count("|") > 1 and line.split("|")[2].strip():
                self._process_ival_section(line_num, line)
            line = line.split("|")[0]
        if func_name != "user_defined":
            hit_words: List[str] = []
            for word in self.rule_words[func_name]:
                # Choose longer word
                if word in line:
                    hit_words.append(word)
            description = line.strip().split(max(hit_words, key=len))
            if description[1] and not description[1].startswith(" "):
                self._raise_detection_error(line_num, line)
        else:
            description = line.strip().split(":")
        return description

    @staticmethod
    def _word2scores(word: str, sentence: str) -> List[float]:
        """
        Calculate similarity scores between word and sentence.

        Parameters
        ----------
        word : str
            User-defined word.
        sentence : str
            Textual unit consisting of two or more words.

        returns
        -------
        scores : list
            List containing similarity scores.

        """
        scores = [
            SequenceMatcher(None, word, sentence[i : i + len(word)]).ratio()
            for i in range(len(sentence) - len(word) + 1)
        ]
        return scores

    def _get_partial_similarity(self, line: str) -> UnregisteredRule:
        """
        Suggest similar rule word when user-defined word is not registered
        in rule_words.

        Parameters
        ----------
        line : str
            Each line of the input text.

        Returns
        -------
        unregistered_rule : UnregisteredRule
            Rule word with the highest similarity score.

        """
        match_words = []
        match_scores = []
        str_subset = []
        for rules in self.rule_words.values():
            for word in rules:
                ratio = self._word2scores(word, line)
                if ratio:
                    match_words.append(word)
                    match_scores.append(max(ratio))
                    str_subset.append(line[np.argmax(ratio) : np.argmax(ratio) + len(word)])
        expected_word = (
            None
            if all([score < self.similarity_threshold for score in match_scores])
            else match_words[np.argmax(match_scores)]
        )
        original_word = (
            None if expected_word is None else str_subset[match_words.index(expected_word)]
        )
        unregistered_rule = UnregisteredRule(expected_word, original_word)

        return unregistered_rule

    @staticmethod
    def _remove_prepositions(sentence: str) -> str:
        """
        Remove preposition from text not to use it for identifying reaction rules.
        """
        for preposition in PREPOSITIONS:
            if sentence.endswith(f" {preposition}"):
                return sentence[: -len(preposition) - 1]
        return sentence

    def _get_arrow_error_message(self, line_num: int) -> str:
        message = (
            f"line{line_num}: Use one of ({', '.join(self.fwd_arrows)}) for unidirectional "
            f"reaction or ({', '.join(self.double_arrows)}) for bi-directional reaction"
        )
        return message

    def _bind_and_dissociate(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'A + B --> AB'  # bind, unidirectional
        >>> 'AB --> A + B'  # dissociate, unidirectional
        >>> 'A + B <--> AB' # bind and dissociate, bidirectional
        """
        for arrow in self._available_arrows():
            if arrow in line:
                params_used = ["kf"] if arrow in self.fwd_arrows else ["kf", "kr"]
                break
        else:
            raise ArrowError(self._get_arrow_error_message(line_num) + ".")
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, *params_used
        )
        is_binding: bool
        is_unidirectional: bool
        for arrow in self._available_arrows():
            if arrow in description[1]:
                is_binding = True
                is_unidirectional = True if arrow in self.fwd_arrows else False
                component1 = description[0].strip(" ")
                component2 = description[1].split(arrow)[0].strip(" ")
                complex = description[1].split(arrow)[1].strip(" ")
                break
            elif arrow in description[0]:
                is_binding = False
                is_unidirectional = True if arrow in self.fwd_arrows else False
                component1 = description[0].split(arrow)[1].strip(" ")
                component2 = description[1].strip(arrow)
                complex = description[0].split(arrow)[0].strip(" ")
                break
        else:
            raise ArrowError(self._get_arrow_error_message(line_num) + ".")
        if component1 == complex or component2 == complex:
            raise ValueError(f"line{line_num:d}: {complex} <- Use a different name.")
        else:
            self._set_species(component1, component2, complex)
            self.complex_formations.append(
                ComplexFormation(line_num, set([component1, component2]), complex, is_binding)
            )
            self.reactions.append(
                f"v[{line_num:d}] = "
                f"x[C.kf{line_num:d}] * y[V.{component1}] * y[V.{component2}]"
                + (f" - x[C.kr{line_num:d}] * y[V.{complex}]" if not is_unidirectional else "")
                if is_binding
                else f"v[{line_num:d}] = "
                f"x[C.kf{line_num:d}] * y[V.{complex}]"
                + (
                    f" - x[C.kr{line_num:d}] * y[V.{component1}] * y[V.{component2}]"
                    if not is_unidirectional
                    else ""
                )
            )
            counter_component1, counter_component2, counter_complex = (0, 0, 0)
            for i, eq in enumerate(self.differential_equations):
                if f"dydt[V.{component1}]" in eq:
                    counter_component1 += 1
                    self.differential_equations[i] = (
                        eq + (" - " if is_binding else " + ") + f"v[{line_num:d}]"
                    )
                elif f"dydt[V.{component2}]" in eq:
                    counter_component2 += 1
                    self.differential_equations[i] = (
                        eq + (" - " if is_binding else " + ") + f"v[{line_num:d}]"
                    )
                elif f"dydt[V.{complex}]" in eq:
                    counter_complex += 1
                    self.differential_equations[i] = (
                        eq + (" + " if is_binding else " - ") + f"v[{line_num:d}]"
                    )
            if counter_component1 == 0:
                self.differential_equations.append(
                    f"dydt[V.{component1}] = - v[{line_num:d}]"
                    if is_binding
                    else f"dydt[V.{component1}] = + v[{line_num:d}]"
                )
            if counter_component2 == 0:
                self.differential_equations.append(
                    f"dydt[V.{component2}] = - v[{line_num:d}]"
                    if is_binding
                    else f"dydt[V.{component2}] = + v[{line_num:d}]"
                )
            if counter_complex == 0:
                self.differential_equations.append(
                    f"dydt[V.{complex}] = + v[{line_num:d}]"
                    if is_binding
                    else f"dydt[V.{complex}] = - v[{line_num:d}]"
                )

    def dimerize(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'A dimerizes <--> AA'
        >>> 'A homodimerizes <--> AA'
        >>> 'A forms a dimer <--> AA'
        >>> 'A forms dimers <--> AA'

        Notes
        -----
        * Parameters
            .. math:: kf, kr

        * Rate equation
            .. math:: v = kf * [A] * [A] - kr * [AA]

        * Differential equation
            .. math::

                d[A]]/dt = - 2 * v

                d[AA]/dt = + v

        """
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, "kf", "kr"
        )
        is_unidirectional: bool
        monomer = description[0].strip(" ")
        for arrow in self._available_arrows():
            if arrow in description[1]:
                is_unidirectional = True if arrow in self.fwd_arrows else False
                dimer = description[1].split(arrow)[1].strip(" ")
                break
        else:
            raise ArrowError(
                self._get_arrow_error_message(line_num) + " to specify the name of the dimer."
            )
        if monomer == dimer:
            raise ValueError(f"{dimer} <- Use a different name.")
        self._set_species(monomer, dimer)
        self.complex_formations.append(ComplexFormation(line_num, set(monomer), dimer, True))
        self.reactions.append(
            f"v[{line_num:d}] = "
            f"x[C.kf{line_num:d}] * y[V.{monomer}] * y[V.{monomer}]"
            + (f" - x[C.kr{line_num:d}] * y[V.{dimer}]" if not is_unidirectional else "")
        )
        counter_monomer, counter_dimer = (0, 0)
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{monomer}]" in eq:
                counter_monomer += 1
                self.differential_equations[i] = eq + f" - 2 * v[{line_num:d}]"
            elif f"dydt[V.{dimer}]" in eq:
                counter_dimer += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_monomer == 0:
            self.differential_equations.append(f"dydt[V.{monomer}] = - v[{line_num:d}]")
        if counter_dimer == 0:
            self.differential_equations.append(f"dydt[V.{dimer}] = + v[{line_num:d}]")

    def bind(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'A binds B <--> AB'
        >>> 'A forms complexes with B <--> AB'

        Notes
        -----
        * Parameters
            .. math:: kf, kr

        * Rate equation
            .. math:: v = kf * [A] * [B] - kr * [AB]

        * Differential equation
            .. math::

                d[A]/dt = - v

                d[B]/dt = - v

                d[AB]/dt = + v

        """
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, "kf", "kr"
        )
        is_unidirectional: bool
        component1 = description[0].strip(" ")
        for arrow in self._available_arrows():
            if arrow in description[1]:
                is_unidirectional = True if arrow in self.fwd_arrows else False
                component2 = description[1].split(arrow)[0].strip(" ")
                complex = description[1].split(arrow)[1].strip(" ")
                break
        else:
            raise ArrowError(
                self._get_arrow_error_message(line_num)
                + " to specify the name of the protein complex."
            )
        if component1 == complex or component2 == complex:
            raise ValueError(f"line{line_num:d}: {complex} <- Use a different name.")
        elif component1 == component2:
            self.dimerize(line_num, line)
        else:
            self._set_species(component1, component2, complex)
            self.complex_formations.append(
                ComplexFormation(line_num, set([component1, component2]), complex, True)
            )
            self.reactions.append(
                f"v[{line_num:d}] = "
                f"x[C.kf{line_num:d}] * y[V.{component1}] * y[V.{component2}]"
                + (f" - x[C.kr{line_num:d}] * y[V.{complex}]" if not is_unidirectional else "")
            )
            counter_component1, counter_component2, counter_complex = (0, 0, 0)
            for i, eq in enumerate(self.differential_equations):
                if f"dydt[V.{component1}]" in eq:
                    counter_component1 += 1
                    self.differential_equations[i] = eq + f" - v[{line_num:d}]"
                elif f"dydt[V.{component2}]" in eq:
                    counter_component2 += 1
                    self.differential_equations[i] = eq + f" - v[{line_num:d}]"
                elif f"dydt[V.{complex}]" in eq:
                    counter_complex += 1
                    self.differential_equations[i] = eq + f" + v[{line_num:d}]"
            if counter_component1 == 0:
                self.differential_equations.append(f"dydt[V.{component1}] = - v[{line_num:d}]")
            if counter_component2 == 0:
                self.differential_equations.append(f"dydt[V.{component2}] = - v[{line_num:d}]")
            if counter_complex == 0:
                self.differential_equations.append(f"dydt[V.{complex}] = + v[{line_num:d}]")

    def dissociate(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'AB dissociates to A and B'
        >>> 'AB is dissociated into A and B'

        Notes
        -----
        * Parameters
            .. math:: kf, kr

        * Rate equation
            .. math:: v = kf * [AB] - kr * [A] * [B]

        * Differential equation
            .. math::

                d[A]/dt = + v

                d[B]/dt = + v

                d[AB]/dt = - v

        """
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, "kf", "kr"
        )
        complex = description[0].strip(" ")
        if " and " not in description[1]:
            raise ValueError(
                f"Use 'and' in line{line_num:d}:\ne.g., AB is dissociated into A and B"
            )
        else:
            component1 = description[1].split(" and ")[0].strip(" ")
            component2 = description[1].split(" and ")[1].strip(" ")
        self._set_species(complex, component1, component2)
        self.complex_formations.append(
            ComplexFormation(line_num, set([component1, component2]), complex, False)
        )
        self.reactions.append(
            f"v[{line_num:d}] = "
            f"x[C.kf{line_num:d}] * y[V.{complex}]"
            f" - x[C.kr{line_num:d}] * y[V.{component1}] * y[V.{component2}]"
        )
        counter_complex, counter_component1, counter_component2 = (0, 0, 0)
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{complex}]" in eq:
                counter_complex += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
            elif f"dydt[V.{component1}]" in eq:
                counter_component1 += 1
                self.differential_equations[i] = (
                    eq + f" + v[{line_num:d}]"
                    if component1 != component2
                    else eq + f" + 2 * v[{line_num:d}]"
                )
            elif f"dydt[V.{component2}]" in eq:
                counter_component2 += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_complex == 0:
            self.differential_equations.append(f"dydt[V.{complex}] = - v[{line_num:d}]")
        if counter_component1 == 0:
            self.differential_equations.append(f"dydt[V.{component1}] = + v[{line_num:d}]")
        if counter_component2 == 0:
            self.differential_equations.append(f"dydt[V.{component2}] = + v[{line_num:d}]")

    def is_phosphorylated(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'uA is phosphorylated <--> pA'

        Notes
        -----
        * Parameters
            .. math:: kf, kr

        * Rate equation
            .. math:: v = kf * [uA] - kr * [pA]

        * Differential equation
            .. math::

                d[uA]/dt = - v

                d[pA]/dt = + v

        """
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, "kf", "kr"
        )
        is_unidirectional: bool
        unphosphorylated_form = description[0].strip(" ")
        for arrow in self._available_arrows():
            if arrow in description[1]:
                is_unidirectional = True if arrow in self.fwd_arrows else False
                phosphorylated_form = description[1].split(arrow)[1].strip(" ")
                break
        else:
            raise ArrowError(
                self._get_arrow_error_message(line_num)
                + " to specify the name of the phosphorylated protein."
            )
        self._set_species(unphosphorylated_form, phosphorylated_form)

        self.reactions.append(
            f"v[{line_num:d}] = "
            f"x[C.kf{line_num:d}] * y[V.{unphosphorylated_form}]"
            + (
                f" - x[C.kr{line_num:d}] * y[V.{phosphorylated_form}]"
                if not is_unidirectional
                else ""
            )
        )
        counter_unphosphorylated_form, counter_phosphorylated_form = (0, 0)
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{unphosphorylated_form}]" in eq:
                counter_unphosphorylated_form += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
            elif "dydt[V.{phosphorylated_form}]" in eq:
                counter_phosphorylated_form += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_unphosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{unphosphorylated_form}] = - v[{line_num:d}]"
            )
        if counter_phosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{phosphorylated_form}] = + v[{line_num:d}]"
            )

    def is_dephosphorylated(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'pA is dephosphorylated --> uA'

        Notes
        -----
        * Parameters
            .. math:: V, K

        * Rate equation
            .. math:: v = V * [pA] / (K + [pA])

        * Differential equation
            .. math::

                d[uA]/dt = + v

                d[pA]/dt = - v

        """
        description = self._preprocessing(sys._getframe().f_code.co_name, line_num, line, "V", "K")
        phosphorylated_form = description[0].strip(" ")
        for arrow in self.fwd_arrows:
            if arrow in description[1]:
                unphosphorylated_form = description[1].split(arrow)[1].strip(" ")
                break
        else:
            raise ArrowError(
                f"line{line_num:d}: Use one of ({', '.join(self.fwd_arrows)}) "
                "to specify the name of the dephosphorylated protein."
            )
        self._set_species(phosphorylated_form, unphosphorylated_form)

        self.reactions.append(
            f"v[{line_num:d}] = "
            f"x[C.V{line_num:d}] * y[V.{phosphorylated_form}] / "
            f"(x[C.K{line_num:d}] + y[V.{phosphorylated_form}])"
        )
        counter_unphosphorylated_form, counter_phosphorylated_form = (0, 0)
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{unphosphorylated_form}]" in eq:
                counter_unphosphorylated_form += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
            elif f"dydt[V.{phosphorylated_form}]" in eq:
                counter_phosphorylated_form += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
        if counter_unphosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{unphosphorylated_form}] = + v[{line_num:d}]"
            )
        if counter_phosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{phosphorylated_form}] = - v[{line_num:d}]"
            )

    def phosphorylate(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'B phosphorylates uA --> pA'

        Notes
        -----
        * Parameters
            .. math:: V, K

        * Rate equation
            .. math:: v = V * [B] * [uA] / (K + [uA])

        * Differential equation
            .. math::

                d[uA]/dt = - v

                d[pA]/dt = + v

        """
        description = self._preprocessing(sys._getframe().f_code.co_name, line_num, line, "V", "K")
        kinase = description[0].strip(" ")
        for arrow in self.fwd_arrows:
            if arrow in description[1]:
                unphosphorylated_form = description[1].split(arrow)[0].strip(" ")
                phosphorylated_form = description[1].split(arrow)[1].strip(" ")
                break
        else:
            raise ValueError(
                f"line{line_num:d}: "
                f"Use one of {', '.join(self.fwd_arrows)} to specify "
                "the name of the phosphorylated (or activated) protein."
            )
        if unphosphorylated_form == phosphorylated_form:
            raise ValueError(f"line{line_num:d}: {phosphorylated_form} <- Use a different name.")
        self._set_species(kinase, unphosphorylated_form, phosphorylated_form)

        self.reactions.append(
            f"v[{line_num:d}] = "
            f"x[C.V{line_num:d}] * y[V.{kinase}] * y[V.{unphosphorylated_form}] / "
            f"(x[C.K{line_num:d}] + y[V.{unphosphorylated_form}])"
        )
        counter_unphosphorylated_form, counter_phosphorylated_form = (0, 0)
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{unphosphorylated_form}]" in eq:
                counter_unphosphorylated_form += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
            elif f"dydt[V.{phosphorylated_form}]" in eq:
                counter_phosphorylated_form += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_unphosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{unphosphorylated_form}] = - v[{line_num:d}]"
            )
        if counter_phosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{phosphorylated_form}] = + v[{line_num:d}]"
            )

    def dephosphorylate(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'B dephosphorylates pA --> uA'

        Notes
        -----
        * Parameters
            .. math:: V, K

        * Rate equation
            .. math:: v = V * [B] * [pA] / (K + [pA])

        * Differential equation
            .. math::

                d[uA]/dt = + v

                d[pA]/dt = - v

        """
        description = self._preprocessing(sys._getframe().f_code.co_name, line_num, line, "V", "K")
        phosphatase = description[0].strip(" ")
        for arrow in self.fwd_arrows:
            if arrow in description[1]:
                phosphorylated_form = description[1].split(arrow)[0].strip(" ")
                unphosphorylated_form = description[1].split(arrow)[1].strip(" ")
                break
        else:
            raise ArrowError(
                f"line{line_num:d}: "
                f"Use one of {', '.join(self.fwd_arrows)} to specify "
                "the name of the dephosphorylated (or deactivated) protein."
            )
        if phosphorylated_form == unphosphorylated_form:
            raise ValueError(f"line{line_num:d}: {unphosphorylated_form} <- Use a different name.")
        self._set_species(phosphatase, phosphorylated_form, unphosphorylated_form)

        self.reactions.append(
            f"v[{line_num:d}] = "
            f"x[C.V{line_num:d}] * y[V.{phosphatase}] * y[V.{phosphorylated_form}] / "
            f"(x[C.K{line_num:d}] + y[V.{phosphorylated_form}])"
        )
        counter_phosphorylated_form, counter_unphosphorylated_form = (0, 0)
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{phosphorylated_form}]" in eq:
                counter_phosphorylated_form += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
            elif f"dydt[V.{unphosphorylated_form}]" in eq:
                counter_unphosphorylated_form += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_phosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{phosphorylated_form}] = - v[{line_num:d}]"
            )
        if counter_unphosphorylated_form == 0:
            self.differential_equations.append(
                f"dydt[V.{unphosphorylated_form}] = + v[{line_num:d}]"
            )

    def transcribe(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'B transcribes a'
        >>> 'B1 & B2 transcribe a'  # (AND-gate)
        >>> 'B transcribes a, repressed by C'  # (Negative regulation)

        Notes
        -----
        * Parameters
            .. math:: V, K, n, (KF, nF)

        * Rate equation
            .. math::

                v = V * [B] ^ {n} / (K ^ {n} + [B] ^ {n})

                v = V * ([B1] * [B2]) ^ {n} / (K ^ {n} + ([B1] * [B2]) ^ {n})

                v = V * [B] ^ {n} / (K ^ {n} + [B] ^ {n} + ([C] / KF) ^ {nF})

        * Differential equation
            .. math:: d[a]/dt = + v

        """
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, "V", "K", "n", "KF", "nF"
        )
        repressor: Optional[str] = None
        ratio = self._word2scores(", repressed by", description[1])
        if not ratio or max(ratio) < 1.0:
            self.parameters.remove(f"KF{line_num:d}")
            self.parameters.remove(f"nF{line_num:d}")
            mRNA = description[1].strip()
            if " " in mRNA:
                # Fix typo in line{line_num:d}
                raise ValueError(
                    f"line{line_num:d}: "
                    "Add ', repressed by XXX' to describe negative regulation from XXX."
                )
        else:
            # Add negative regulation from repressor
            mRNA = description[1].split(", repressed by")[0].strip()
            repressor = description[1].split(", repressed by")[1].strip()
        if " & " not in description[0]:
            TF = description[0].strip(" ")
            self._set_species(mRNA, TF)
            if repressor is not None:
                self._set_species(repressor)
            self.reactions.append(
                f"v[{line_num:d}] = "
                f"x[C.V{line_num:d}] * y[V.{TF}] ** x[C.n{line_num:d}] / "
                f"(x[C.K{line_num:d}] ** x[C.n{line_num:d}] + "
                f"y[V.{TF}] ** x[C.n{line_num:d}]"
                + (
                    ")"
                    if repressor is None
                    else f" + (y[V.{repressor}] / x[C.KF{line_num:d}]) ** x[C.nF{line_num:d}])"
                )
            )
        else:
            # AND-gate
            TFs = [TF.strip(" ") for TF in description[0].split(" & ")]
            self._set_species(mRNA, *TFs)
            if repressor is not None:
                self._set_species(repressor)
            self.reactions.append(
                f"v[{line_num:d}] = "
                f"x[C.V{line_num:d}] * ({'y[V.' + '] * y[V.'.join(TFs) + ']'}) ** x[C.n{line_num:d}] / "
                f"(x[C.K{line_num:d}] ** x[C.n{line_num:d}] + "
                f"({'y[V.' + '] * y[V.'.join(TFs) + ']'}) ** x[C.n{line_num:d}]"
                + (
                    ")"
                    if repressor is None
                    else f" + (y[V.{repressor}] / x[C.KF{line_num:d}]) ** x[C.nF{line_num:d}])"
                )
            )
        counter_mRNA = 0
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{mRNA}]" in eq:
                counter_mRNA += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_mRNA == 0:
            self.differential_equations.append(f"dydt[V.{mRNA}] = + v[{line_num:d}]")

    def synthesize(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'B synthesizes A'

        Notes
        -----
        * Parameters
            .. math:: kf

        * Rate equation
            .. math:: v = kf * [B]

        * Differential equation
            .. math:: d[A]/dt = + v

        """
        description = self._preprocessing(sys._getframe().f_code.co_name, line_num, line, "kf")
        catalyst = description[0].strip(" ")
        product = description[1].strip(" ")
        self._set_species(catalyst, product)
        self.reactions.append(f"v[{line_num:d}] = x[C.kf{line_num:d}] * y[V.{catalyst}]")
        counter_product = 0
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{product}]" in eq:
                counter_product += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_product == 0:
            self.differential_equations.append(f"dydt[V.{product}] = + v[{line_num:d}]")

    def is_synthesized(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'A is synthesized'

        Notes
        -----
        * Parameters
            .. math:: kf

        * Rate equation
            .. math:: v = kf

        * Differential equation
            .. math:: d[A]/dt = + v

        """
        description = self._preprocessing(sys._getframe().f_code.co_name, line_num, line, "kf")
        chemical_species = description[0].strip(" ")
        self._set_species(chemical_species)
        self.reactions.append(f"v[{line_num:d}] = x[C.kf{line_num:d}]")
        counter_chemical_species = 0
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{chemical_species}]" in eq:
                counter_chemical_species += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_chemical_species == 0:
            self.differential_equations.append(f"dydt[V.{chemical_species}] = + v[{line_num:d}]")

    def degrade(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'B degrades A'

        Notes
        -----
        * Parameters
            .. math:: kf

        * Rate equation
            .. math:: v = kf * [B] * [A]

        * Differential equation
            .. math:: d[A]/dt = - v

        """
        description = self._preprocessing(sys._getframe().f_code.co_name, line_num, line, "kf")
        protease = description[0].strip(" ")
        protein = description[1].strip(" ")
        self._set_species(protease, protein)
        self.reactions.append(
            f"v[{line_num:d}] = x[C.kf{line_num:d}] * y[V.{protease}] * y[V.{protein}]"
        )
        counter_protein = 0
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{protein}]" in eq:
                counter_protein += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
        if counter_protein == 0:
            self.differential_equations.append(f"dydt[V.{protein}] = - v[{line_num:d}]")

    def is_degraded(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'A is degraded'

        Notes
        -----
        * Parameters
            .. math:: kf

        * Rate equation
            .. math:: v = kf * [A]

        * Differential equation
            .. math:: d[A]/dt = - v

        """
        description = self._preprocessing(sys._getframe().f_code.co_name, line_num, line, "kf")
        chemical_species = description[0].strip(" ")
        self._set_species(chemical_species)
        self.reactions.append(f"v[{line_num:d}] = x[C.kf{line_num:d}] * y[V.{chemical_species}]")
        counter_chemical_species = 0
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{chemical_species}]" in eq:
                counter_chemical_species += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
        if counter_chemical_species == 0:
            self.differential_equations.append(f"dydt[V.{chemical_species}] = - v[{line_num:d}]")

    def translocate(self, line_num: int, line: str) -> None:
        r"""
        Examples
        --------
        >>> 'A_at_cyt translocates from cytoplasm to nucleus (V_cyt, V_nuc) <--> A_at_nuc'
        >>> 'A_at_cyt is translocated from cytoplasm to nucleus (V_cyt, V_nuc) <--> A_at_nuc'

        Notes
        -----
        * Parameters
            .. math:: kf, kr, (V_{pre}, V_{post})

        * Rate equation
            .. math:: v = kf * [A\_at\_pre] - kr * (V_{post} / V_{pre}) * [A\_at\_post]

        * Differential equation
            .. math::

                d[A\_at\_pre]/dt = - v

                d[A\_at\_post]/dt = + v * (V_{pre} / V_{post})

        """
        for arrow in self._available_arrows():
            if arrow in line:
                is_unidirectional = True if arrow in self.fwd_arrows else False
                params_used = ["kf"] if is_unidirectional else ["kf", "kr"]
                break
        else:
            raise ArrowError(self._get_arrow_error_message(line_num))
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, *params_used
        )
        pre_translocation = description[0].strip(" ")
        for arrow in self._available_arrows():
            if arrow in description[1]:
                post_translocation = description[1].split(arrow)[1].strip(" ")
                break
        else:
            assert False
        if pre_translocation == post_translocation:
            raise ValueError(f"line{line_num:d}: {post_translocation} <- Use a different name.")
        # Information about compartment volumes
        if "(" in description[1] and ")" in description[1]:
            [pre_volume, post_volume] = description[1].split("(")[-1].split(")")[0].split(",")
            if not self._isfloat(pre_volume.strip(" ")) or not self._isfloat(
                post_volume.strip(" ")
            ):
                raise ValueError("pre_volume and post_volume must be float or int.")
        else:
            [pre_volume, post_volume] = ["1", "1"]
        self._set_species(pre_translocation, post_translocation)
        self.reactions.append(
            f"v[{line_num:d}] = x[C.kf{line_num:d}] * y[V.{pre_translocation}]"
            + (
                f" - x[C.kr{line_num:d}] * y[V.{post_translocation}]"
                if not is_unidirectional
                else ""
            )
        )
        if float(pre_volume.strip(" ")) != float(post_volume.strip(" ")):
            self.reactions[
                -1
            ] = f"v[{line_num:d}] = " f"x[C.kf{line_num:d}] * y[V.{pre_translocation}]" + (
                f" - x[C.kr{line_num:d}] * "
                f"({post_volume.strip()} / {pre_volume.strip()}) * "
                f"y[V.{post_translocation}]"
                if not is_unidirectional
                else ""
            )
        counter_pre_translocation, counter_post_translocation = (0, 0)
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{pre_translocation}]" in eq:
                counter_pre_translocation += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
            elif f"dydt[V.{post_translocation}]" in eq:
                counter_post_translocation += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
                if float(pre_volume.strip(" ")) != float(post_volume.strip(" ")):
                    self.differential_equations[
                        i
                    ] += f" * ({pre_volume.strip()} / {post_volume.strip()})"
        if counter_pre_translocation == 0:
            self.differential_equations.append(f"dydt[V.{pre_translocation}] = - v[{line_num:d}]")
        if counter_post_translocation == 0:
            self.differential_equations.append(f"dydt[V.{post_translocation}] = + v[{line_num:d}]")
            if float(pre_volume.strip(" ")) != float(post_volume.strip(" ")):
                self.differential_equations[
                    -1
                ] += f" * ({pre_volume.strip()} / {post_volume.strip()})"

    def user_defined(self, line_num: int, line: str) -> None:
        """
        Examples
        --------
        >>> 'Reactant --> Product: define rate equation here'

        Notes
        -----
        * Use p[xxx] and u[xxx] for describing parameters and species, respectively.

        * Use '0' or '∅' for degradation/creation to/from nothing.

        * Differential equation
            .. math::

                d[Reactant]/dt = - v

                d[Product]/dt = + v
        """
        all_params = re.findall(r"p\[(.*?)\]", line)
        description = self._preprocessing(
            sys._getframe().f_code.co_name, line_num, line, *all_params
        )
        balance = description[0].strip()
        rate_equation = description[1].strip()
        for arrow in self.fwd_arrows:
            if arrow in balance:
                reactant, product = balance.split(arrow)
                self._set_species(reactant.strip(), product.strip())
                break
        else:
            raise ArrowError(f"line{line_num:d}: Use one of {', '.join(self.fwd_arrows)}.")
        rate_equation = (
            rate_equation.replace("p[", "x[C.").replace("u[", "y[V.").replace("^", "**")
        )
        self.reactions.append(f"v[{line_num:d}] = " + rate_equation.strip())
        counter_reactant = 0
        counter_product = 0
        for i, eq in enumerate(self.differential_equations):
            if f"dydt[V.{reactant}]" in eq and reactant not in self.nothing:
                counter_reactant += 1
                self.differential_equations[i] = eq + f" - v[{line_num:d}]"
            elif f"dydt[V.{product}]" in eq and product not in self.nothing:
                counter_product += 1
                self.differential_equations[i] = eq + f" + v[{line_num:d}]"
        if counter_reactant == 0 and reactant not in self.nothing:
            self.differential_equations.append(f"dydt[V.{reactant}] = - v[{line_num:d}]")
        if counter_product == 0 and product not in self.nothing:
            self.differential_equations.append(f"dydt[V.{product}] = + v[{line_num:d}]")

    def _extract_event(self, line_num: int, line: str):
        # About biochemical event
        if line.startswith("@rxn "):
            line = self._remove_prefix(line, "@rxn ")
            if line.count(":") != 1:
                raise SyntaxError(f"line{line_num:d}: Missing colon")
            else:
                self.user_defined(line_num, line)
        # About observables
        elif line.startswith("@obs "):
            line = self._remove_prefix(line, "@obs ")
            if line.count(":") != 1:
                raise SyntaxError(
                    f"line{line_num:d}: Missing colon\n"
                    "Should be `@obs <observable name>: <expression>`."
                )
            else:
                self.obs_desc.append(line.split(":"))
        # About simulation info.
        elif line.startswith("@sim "):
            line = self._remove_prefix(line, "@sim ")
            if line.count(":") != 1:
                raise SyntaxError(f"line{line_num:d}: Missing colon")
            else:
                if line.startswith("tspan"):
                    t_info = line.split(":")[-1].strip()
                    if "[" in t_info and "]" in t_info:
                        [t0, tf] = t_info.split("[")[-1].split("]")[0].split(",")
                        if t0.strip(" ").isdecimal() and tf.strip(" ").isdecimal():
                            self.sim_tspan.append(t0)
                            self.sim_tspan.append(tf)
                        else:
                            raise TypeError("@sim tspan: [t0, tf] must be a list of integers.")
                    else:
                        raise ValueError(
                            "`tspan` must be a two element vector [t0, tf] "
                            "specifying the initial and final times."
                        )
                elif line.startswith("unperturbed"):
                    self.sim_unperturbed += line.split(":")[-1].strip()
                elif line.startswith("condition "):
                    self.sim_conditions.append(self._remove_prefix(line, "condition ").split(":"))
                else:
                    raise ValueError(
                        f"(line{line_num:d}) Available options are: "
                        "'@sim tspan:', '@sim unperturbed:', or '@sim condition XXX:'."
                    )
        # Additional species
        elif line.startswith("@add "):
            line = self._remove_prefix(line, "@add ")
            if line.startswith("species "):
                line = self._remove_prefix(line, "species ")
                new_species = line.strip()
                if new_species not in self.species:
                    self._set_species(new_species)
                else:
                    raise NameError(f"{new_species} is already defined.")
            elif line.startswith("param "):
                line = self._remove_prefix(line, "param ")
                new_param = line.strip()
                if new_param not in self.parameters:
                    self._set_params(None, None, new_param)
                    self.param_excluded.append(new_param)
                else:
                    raise NameError(f"{new_param} is already defined.")
            else:
                raise ValueError(f"(line{line_num:d}) Must be either @add param or @add species.")
        else:
            raise ValueError("Available symbols are: @rxn, @add, @obs, @sim.")

    def create_ode(self) -> None:
        """
        Find a keyword in each line to identify the reaction rule and
        construct an ODE model.

        """
        with open(self.input_txt, encoding="utf-8") as f:
            lines = f.readlines()
        for line_num, line in enumerate(lines, start=1):
            # Remove double spaces
            while True:
                if "  " not in line:
                    break
                else:
                    line = line.replace("  ", " ")
            # Comment out
            line = line.split("#")[0].rstrip(" ")
            if not line.strip():
                # Skip blank lines
                continue
            elif lines.count(line) > 1:
                # Find duplicate lines
                raise DuplicateError(
                    f"Reaction '{line}' is duplicated in lines "
                    + ", ".join([str(i + 1) for i, rxn in enumerate(lines) if rxn == line])
                )
            elif line.startswith("@"):
                self._extract_event(line_num, line)
            # Detect reaction rule
            else:
                for reaction_rule, words in self.rule_words.items():
                    if any([self._remove_prepositions(word) in line for word in words]):
                        exec("self." + reaction_rule + "(line_num, line)")
                        break
                else:
                    self._raise_detection_error(line_num, line)
