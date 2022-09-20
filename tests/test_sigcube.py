import pytest
import os
import numpy as np

PATH = os.path.dirname(__file__)

FILE_0 = PATH + "/fixtures/sigcube.npy"
FILE_1 = PATH  + "/../new_sigcube.npy"

@pytest.fixture
def original_sigcube():
    ori_arr = np.load(FILE_0)
    return ori_arr

def test_sigcube(original_sigcube) -> None:
    new_sigcube = np.load(FILE_1)
    assert np.array_equal(original_sigcube, new_sigcube)