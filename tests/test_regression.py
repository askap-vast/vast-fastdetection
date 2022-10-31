import os
import pytest
import pandas as pd
from pandas.testing import assert_frame_equal

# result file to compare with
PATH = os.path.dirname(__file__)
RES_FILE_1 = PATH + "/fixtures/detection_results/output_cand.csv"
RES_FILE_3 = PATH + "/fixtures/"

# generated data
RES_FILE_2 = PATH + "/../output/output_beam14_cand.csv"


@pytest.fixture
def result_df():
    df = pd.read_csv(RES_FILE_1)
    return df 

def test_regression_run_cube(result_df) -> None:
    df_new = pd.read_csv(RES_FILE_2)
    assert_frame_equal(result_df, df_new)

