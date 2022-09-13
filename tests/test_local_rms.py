import pytest
import numpy as np
import os

PATH = os.path.dirname(__file__)

# data to compare with
FILE_1 = PATH + "/fixtures/rmscube.npy"

FILE_2 = PATH + "/../new_rmscube.npy"

FILE_3 = PATH + "/../rms_0.npy"

@pytest.fixture
def original_rmscube():
    rmscube = np.load(FILE_1)
    return rmscube

@pytest.fixture
def original_rms_0():
    rmscube = np.load(FILE_1)
    return rmscube[0,:,:]

def test_cube_local_rms(original_rmscube) -> None:
    new_rmscube = np.load(FILE_2)
    assert np.array_equal(original_rmscube, new_rmscube)

def test_get_local_rms(original_rms_0) -> None:
    rms_0 = np.load(FILE_3)
    assert np.array_equal(original_rms_0, rms_0)
