from collections import Counter
from dataclasses import dataclass, field
from typing import List, NamedTuple, Tuple


class DuplicateError(Exception):
    pass


class ComplexFormation(NamedTuple):
    rxn_no: int
    components: set
    complex: str
    is_binding: bool


@dataclass
class ThermodynamicRestrictions(object):
    """
    Thermodynamic restrictions along cyclic pathways in a kinetic scheme.

    Notes
    -----
    If a kinetic scheme includes "true" cycles, in which the initial and final states are
    identical, the equilibrium constants of the reactions along any cycle satisfy so-called
    "detailed balance" relationships. These detailed balance relations require the product
    of the equilibrium constants along a cycle to be equal to 1, since at equilibrium the net flux
    through any cycle vanishes.
    """

    _rxn_indices: dict = field(
        default_factory=dict,
        init=False,
    )
    _tree: dict = field(
        default_factory=dict,
        init=False,
    )
    restrictions: list = field(
        default_factory=list,
        init=False,
    )
    complex_formations: List[ComplexFormation] = field(
        default_factory=list,
        init=False,
    )

    def _create_tree(
        self,
        tree: dict,
        parent_pattern: ComplexFormation,
    ) -> dict:
        """
        Create tree.
        """
        for component in parent_pattern.components:
            for another_pattern in self.complex_formations:
                if (
                    component == another_pattern.complex
                    and parent_pattern.rxn_no != another_pattern.rxn_no
                ):
                    for k in another_pattern.components:
                        if component not in tree.keys():
                            tree.setdefault(component, {})
                        tree[component].setdefault(k, {f"{another_pattern.rxn_no}": {}})
                        tree[component][k] = self._create_tree(tree[component][k], another_pattern)
        return tree

    def _add_to_rxn_indices1(self, complex_name: str, tree: dict) -> None:
        for component in tree.keys():
            if not tree[component]:
                self._rxn_indices[complex_name].append(component)
            else:
                self._add_to_rxn_indices1(complex_name, tree[component])

    def _add_to_rxn_indices2(self, complex_name: str, monomer: str, pair: str) -> None:
        for another_complex in self._tree[complex_name].keys():
            if (
                monomer in self._tree[complex_name][another_complex].keys()
                and len(self._tree[complex_name][another_complex][monomer]) > 1
            ):
                for name in self._tree[complex_name][another_complex][monomer].keys():
                    if name in self._tree[complex_name][pair].keys():
                        self._append_reaction_number(
                            complex_name,
                            monomer,
                            pair,
                            name,
                            another_complex,
                        )

    def _add_to_rxn_indices3(self, complex_name: str, monomer: str, pair: str) -> None:
        for another_complex in self._tree[complex_name].keys():
            if (
                monomer in self._tree[complex_name][another_complex].keys()
                and len(self._tree[complex_name][another_complex][monomer]) == 1
            ):
                for name in self._tree[complex_name][another_complex].keys():
                    if not name.isdecimal() and name != monomer:
                        for components in self._tree[complex_name][another_complex].keys():
                            if name == components:
                                for k in self._tree[complex_name][pair][name].keys():
                                    if not k.isdecimal() and k in self._tree[complex_name].keys():
                                        self._append_reaction_number(
                                            complex_name,
                                            monomer,
                                            pair,
                                            name,
                                            another_complex,
                                        )

    def _append_reaction_number(
        self,
        complex_name: str,
        monomer: str,
        pair: str,
        name: str,
        another_complex: str,
    ) -> None:
        for val in self._tree[complex_name][another_complex].keys():
            if val.isdecimal():
                self._rxn_indices[monomer].append(val)
        for val in self._tree[complex_name][another_complex][monomer].keys():
            if val.isdecimal():
                self._rxn_indices[monomer].append(val)
        for val in self._tree[complex_name][pair][name].keys():
            if val.isdecimal():
                self._rxn_indices[monomer].append(val)

    def _get_complex_patterns(self) -> List[Tuple[ComplexFormation, ComplexFormation]]:
        complex_patterns = []
        for i, pattern_a in enumerate(self.complex_formations):
            for j in range(i + 1, len(self.complex_formations)):
                pattern_b = self.complex_formations[j]
                if (
                    pattern_a.components == pattern_b.components
                    and pattern_a.complex == pattern_b.complex
                ):
                    raise DuplicateError(
                        "Duplicate binding-dissociation events are detected "
                        f"in lines {pattern_a.rxn_no:d} and {pattern_b.rxn_no:d}."
                    )
                elif (
                    pattern_a.is_binding != pattern_b.is_binding
                    and pattern_a.complex == pattern_b.complex
                ):
                    complex_patterns.append((pattern_a, pattern_b))
        return complex_patterns

    def _impose_restrictions(self, complex_name: str) -> None:
        count = Counter(self._rxn_indices[complex_name])
        if all(n == 2 for n in count.values()):
            self.restrictions.append(list(set(self._rxn_indices[complex_name])))
        else:
            _monomers = {}
            for name in self._tree[complex_name].keys():
                if len(self._tree[complex_name][name].values()) == 1:
                    _monomers[name] = list(self._tree[complex_name][name].keys())[0]
            for monomer, ridx in _monomers.items():
                self._rxn_indices[monomer] = [ridx]
                pair = None
                for another_complex in self._tree[complex_name].keys():
                    if (
                        another_complex != monomer
                        and ridx in self._tree[complex_name][another_complex]
                    ):
                        pair = another_complex
                if pair is not None:
                    self._add_to_rxn_indices2(complex_name, monomer, pair)
                    _reactions = list(set(self._rxn_indices[monomer]))
                    if len(_reactions) > 2 and _reactions not in self.restrictions:
                        self.restrictions.append(_reactions)
            # initialize monomers
            for monomer, ridx in _monomers.items():
                self._rxn_indices[monomer] = [ridx]
                pair = None
                for another_complex in self._tree[complex_name].keys():
                    if (
                        another_complex != monomer
                        and ridx in self._tree[complex_name][another_complex]
                    ):
                        pair = another_complex
                if pair is not None:
                    self._add_to_rxn_indices3(complex_name, monomer, pair)
                    _reactions = list(set(self._rxn_indices[monomer]))
                    if len(_reactions) > 2 and _reactions not in self.restrictions:
                        self.restrictions.append(_reactions)

    def find_cyclic_reaction_routes(self) -> None:
        """
        Find cyclic pathways in a reaction network.
        """
        complex_patterns = self._get_complex_patterns()
        _tree = {}
        for patterns in complex_patterns:
            for pattern in patterns:
                _tree.setdefault(pattern.complex, {})
                for component in pattern.components:
                    _tree[pattern.complex].setdefault(component, {f"{pattern.rxn_no}": {}})
                _tree[pattern.complex] = self._create_tree(_tree[pattern.complex], pattern)
        self._tree = _tree
        for complex_name in self._tree.keys():
            self._rxn_indices.setdefault(complex_name, [])
            self._add_to_rxn_indices1(complex_name, self._tree[complex_name])
            self._impose_restrictions(complex_name)
