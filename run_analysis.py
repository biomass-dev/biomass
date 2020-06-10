from biomass.models.Nakakuki_Cell_2010 import *

if __name__ == '__main__':
    from biomass.analysis.reaction import ReactionSensitivity
    reaction = ReactionSensitivity(
        reaction_system=set_model,
        obs=observables,
        sim=NumericalSimulation(),
        sp=SearchParam(),
        rxn=ReactionNetwork()
    )
    reaction.analyze(metric='integral', style='barplot')

    """
    from biomass.analysis.nonzero_init import NonZeroInitSensitivity
    nonzero_init = NonZeroInitSensitivity(
        species = V.species,
        ival = initial_values,
        obs = observables,
        sim = NumericalSimulation(),
        sp = SearchParam()
    )
    nonzero_init.analyze(metric='integral', style='barplot')
    """