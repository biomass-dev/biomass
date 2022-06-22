import os, shutil
import pytest
from biomass import create_model
def test_graph():
    model_lib_list = [name for name in os.listdir('biomass/models') if '__' not in name]
    if not os.path.isdir(os.path.join('biomass/models', '__graphs')):
        os.mkdir(os.path.join('biomass/models', '__graphs'))
    for model_name in model_lib_list:
        if model_name in ['prolif_quies', 'g1s_transition']:
            model = create_model('.'.join(('biomass', 'models', model_name)))
            with pytest.warns(UserWarning):
                model.to_graph(os.path.join('biomass/models', '__graphs', model_name + 'test.png'))
                assert os.stat(os.path.join('biomass/models', '__graphs', model_name + 'test.png')).st_size > 1024*10
            continue
        model = create_model('.'.join(('biomass', 'models', model_name)))
        model.to_graph(os.path.join('biomass/models', '__graphs', model_name + 'test.png'))
        assert os.stat(os.path.join('biomass/models', '__graphs', model_name + 'test.png')).st_size > 1024*10
        
def test_cleanup():
    if os.path.isdir(os.path.join('biomass/models', '__graphs')):
        shutil.rmtree(os.path.join('biomass/models', '__graphs'))