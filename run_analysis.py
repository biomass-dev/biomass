import warnings

from biomass.models.Nakakuki_Cell_2010 import *
from biomass.models.Nakakuki_Cell_2010 import __path__ as MODEL_PATH

if __name__ == '__main__':
    from biomass.analysis.reaction import ReactionSensitivity
    reaction = ReactionSensitivity(
        model_path=MODEL_PATH[0],
        reaction_system=set_model,
        obs=observables,
        sim=NumericalSimulation(),
        sp=SearchParam(),
        rxn=ReactionNetwork()
    )
    warnings.filterwarnings('ignore')
    reaction.analyze(metric='integral', style='barplot')

    """
    from biomass.analysis.nonzero_init import NonZeroInitSensitivity
    nonzero_init = NonZeroInitSensitivity(
        model_path=MODEL_PATH[0],
        species=V.NAMES,
        ival=initial_values,
        obs=observables,
        sim=NumericalSimulation(),
        sp=SearchParam()
    )
    warnings.filterwarnings('ignore')
    nonzero_init.analyze(metric='integral', style='barplot')
    """