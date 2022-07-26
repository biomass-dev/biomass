from sbmlutils.io import write_sbml, read_sbml
import sbmlutils.factory as sbml
from biomass import Text2Model
import tempfile, os
from pathlib import Path
import re
from typing import NamedTuple


class ExportModel(object):
    def __init__(
        self,
        model: Text2Model,
    ) -> None:
        """
        Transfer biomass model to SBML format.

        Attributes
        -----------
            model (Text2Model):
                Model that will be transferred to SBML.

        Examples
        --------
        >>> model = Text2Model('model.txt')
        >>> biomass2sbml = ExportModel(model)
        >>> biomass2sbml.model2sbml()
        >>> biomass2sbml.safe_sbml('home', 'model.xml')
        """
        self.model = model
        self.sbmlmodel = None

    def _get_species(
        self,
        species: str,
    ) -> sbml.Species:
        """Helper to generate species object

        Parameters
        -----------
            species (str): Name of the species

        Returns:
            sbml.Species: sbmlutils species object for model creation
        """
        match = "y0[V." + species + "]"
        y0 = 0
        for init in self.model.init_info:
            if match in init:
                y0 = float(init.split("=")[1])
                break
        sbmlspecies = sbml.Species(
            sid=species,
            name=species,
            initialConcentration=y0,
            compartment="C",
        )
        return sbmlspecies

    def _get_param(
        self,
        param: str,
        _recurse: bool = False,
    ) -> sbml.Parameter:
        """Helper function to generate parameter object

        Args:
            param (str): Name of parameter
            _recurse (bool, optional): Used internally to recursively search for parameter references. Defaults to False.

        Returns:
            sbml.Parameter: sbmlutils Parameter object
        """
        match = "x[C." + param + "]"
        paramval = 1

        for parameter in self.model.param_info:
            if match in parameter.split("=")[0]:
                try:
                    paramval = float(parameter.split("=")[1])
                except ValueError:
                    reg = r"(?<=x\[C\.)(.+)(?=\])"
                    righthand = re.search(reg, parameter.split("=")[1]).group()
                    paramval = self._get_param(righthand, _recurse=True)
                break
        if _recurse is True:
            return paramval
        else:
            sbmlparam = sbml.Parameter(
                sid=param,
                name=param,
                value=paramval,
            )
            return sbmlparam

    def _get_reaction(
        self,
        num: int,
        kinetics: NamedTuple,
    ) -> sbml.Reaction:
        """Helper function to generate reaction object

        Parameters
        ----------
            num (int): Nr of the reaction, used for naming
            kinetics (NamedTuple): Participants in reaction and rate law

        Returns:
            sbml.Reaction: sbmlutils Reaction object
        """
        lefthand = " + ".join(kinetics.reactants)
        righthand = " + ".join(kinetics.products)
        equation = lefthand + " => " + righthand
        if kinetics.modifiers:
            modifiers = "[" + ",".join(kinetics.modifiers) + "]"
            equation = equation + " " + modifiers
        rate = kinetics.rate.replace("**", "^")
        reaction = sbml.Reaction(
            sid="R" + str(num), name="R" + str(num), equation=equation, formula=(rate, None)
        )
        return reaction

    def model2sbml(self) -> None:
        """Uses helper functions to get model species, parameters and reactions to generate SBML model"""
        self.sbmlmodel = sbml.Model(
            "menten_model",
            name="menten_model",
            units=sbml.Units,
            model_units=sbml.ModelUnits(
                time=sbml.Units.second,
                substance=sbml.Units.mole,
                extent=sbml.Units.mole,
                volume=sbml.Units.litre,
            ),
            compartments=[sbml.Compartment(sid="C", value=1.0)],
            species=[self._get_species(species) for species in self.model.species],
            parameters=[self._get_param(param) for param in self.model.parameters],
            reactions=[
                self._get_reaction(i, kinetics) for i, kinetics in enumerate(self.model.kinetics)
            ],
        )

    def save_sbml(
        self,
        savepath: str,
        filename: str,
    ) -> None:
        """Generates SBML file and validates it.

        Args:
            savepath (str): Directory in which it will be stored.
            filename (str): Output filename
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            results = sbml.create_model(
                models=self.sbmlmodel,
                output_dir=Path(tmp_path),
                tmp=False,
                units_consistency=False,
                sbml_level=3,
                sbml_version=1,
            )
            doc = read_sbml(source=results.sbml_path, validate=False)
        write_sbml(doc, Path(os.path.join(savepath, filename)))
