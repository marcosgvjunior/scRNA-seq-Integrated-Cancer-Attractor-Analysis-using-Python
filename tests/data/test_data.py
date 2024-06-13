import pytest
import sys
import os

current_dir = os.path.abspath(os.path.dirname(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.append(project_root)

from src.data.data_module import Data

import pandas as pd
import numpy as np

@pytest.fixture
def load_data_fixture():
    data = Data().count_matrix
    yield data
    del data

# @pytest.mark.xfail
def test_data_type(load_data_fixture):
    msg = "Incorrect data type!"
    assert type(load_data_fixture) == np.ndarray, msg

   






