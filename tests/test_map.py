import pytest
import numpy as np
import os

PATH = os.path.dirname(__file__)
FILE_0 = PATH + "/fixtures/chimap.npy"
FILE_1 = PATH + "/../chimap.npy"

FILE_2 = PATH + "/fixtures/peakmap.npy"
FILE_3 = PATH + "/../peakmap.npy"

FILE_4 = PATH + "/fixtures/stdmap.npy"
FILE_5 = PATH + "/../stdmap.npy"

@pytest.fixture
def original_chimap():
    o_chimap = np.load(FILE_0)
    return o_chimap

def test_chimap(original_chimap) -> None:
    # update from numpy.power to numpy.square, the result is not exactly the same, but very close to
    new_chimap = np.load(FILE_1)
    assert np.allclose(original_chimap, new_chimap)

@pytest.fixture
def original_peakmap():
    o_peakmap = np.load(FILE_2)
    return o_peakmap

def test_peakmap(original_peakmap) -> None:
    new_peakmap = np.load(FILE_3)
    assert np.array_equal(original_peakmap, new_peakmap)

@pytest.fixture
def original_stdmap():
    o_stdmap = np.load(FILE_4)
    return o_stdmap

def test_stdmap(original_stdmap) -> None:
    # dask.array.nanstd produces slightly different result to numpy.nanstd
    new_stdmap = np.load(FILE_5)
    assert np.allclose(original_stdmap, new_stdmap)
