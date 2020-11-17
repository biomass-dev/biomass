import os
import shutil
from distutils.dir_util import copy_tree
import warnings
from biomass import run_simulation, run_analysis

import pytest

copy_tree('tests/out', 'biomass/models/Nakakuki_Cell_2010/out')
from biomass.models import Nakakuki_Cell_2010

MODEL_PATH = Nakakuki_Cell_2010.__path__[0]

if os.path.isdir(MODEL_PATH + '/figure'):
    shutil.rmtree(MODEL_PATH + '/figure')
    

def test_run_simulation():
    run_simulation(Nakakuki_Cell_2010, viz_type='average', stdev=True)
    assert os.path.isfile(
        MODEL_PATH + '/simulation_data/simulations_all.npy'
    )

'''
def test_sensitivity_analysis():
    run_analysis(
        Nakakuki_Cell_2010, target='initial_condition', metric='integral'
    )
    assert os.path.isfile(
        MODEL_PATH
        + '/sensitivity_coefficients/initial_condition/integral/sc.npy'
    )
'''