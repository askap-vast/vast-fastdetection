import pytest
import numpy as np
import os

PATH = os.path.dirname(__file__)
FILE_0 = PATH + "/fixtures/chimap.npy"
FILE_1 = PATH + "/../chimap.npy"

@pytest.fixture
def original_chimap():
    o_chimap = np.load(FILE_0)
    return o_chimap

def test_chimap(original_chimap) -> None:
    new_chimap = np.load(FILE_1)
    assert np.array_equal(original_chimap, new_chimap)

