import os
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal

# result file to compare with
PATH = os.path.dirname(__file__)
RES_FILE_1a = PATH + "/fixtures/beam14/output_chisquare_cand.csv"
RES_FILE_2a = PATH + "/fixtures/beam14/output_peak_cand.csv"

RES_FILE_4a = PATH + "/fixtures/beam14/output_final.csv"
RES_FILE_5a = PATH + "/fixtures/beam14/output_lightcurve_local_rms.csv"
RES_FILE_6a = PATH + "/fixtures/beam14/output_lightcurve_peak_flux.csv"

# generated data
RES_FILE_1b = PATH + "/../output/output_beam14_chisquare_cand.csv"
RES_FILE_2b = PATH + "/../output/output_beam14_peak_cand.csv"

RES_FILE_4b = PATH + "/../output/output_beam14_final.csv"
RES_FILE_5b = PATH + "/../output/output_beam14_lightcurve_local_rms.csv"
RES_FILE_6b = PATH + "/../output/output_beam14_lightcurve_peak_flux.csv"


@pytest.fixture
def result_chi():
    df = pd.read_csv(RES_FILE_1a)
    return df 

@pytest.fixture
def result_peak():
    df = pd.read_csv(RES_FILE_2a)
    return df 

@pytest.fixture
def result_final():
    df = pd.read_csv(RES_FILE_4a)
    return df 

@pytest.fixture
def result_local_rms():
    df = pd.read_csv(RES_FILE_5a)
    return df 

@pytest.fixture
def result_flux():
    df = pd.read_csv(RES_FILE_6a)
    return df 

def test_regression_chi(result_chi) -> None:
    df_new = pd.read_csv(RES_FILE_1b)
    assert_frame_equal(result_chi, df_new)

def test_regression_peak(result_peak) -> None:
    df_new = pd.read_csv(RES_FILE_2b)
    assert_frame_equal(result_peak, df_new)

def test_regression_final(result_final) -> None:
    df_new = pd.read_csv(RES_FILE_4b)
    assert_frame_equal(result_final, df_new)

def test_regression_rms(result_local_rms) -> None:
    df_new = pd.read_csv(RES_FILE_5b)
    assert_frame_equal(result_local_rms, df_new)

def test_regression_flux(result_flux) -> None:
    df_new = pd.read_csv(RES_FILE_6b)
    assert_frame_equal(result_flux, df_new)






