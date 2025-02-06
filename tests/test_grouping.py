import pytest
import pandas as pd
import numpy as np
from MSCI.Grouping_MS1.Grouping_mw_irt import make_data_compatible, within_ppm, within_tolerance, find_combinations_kdtree, process_peptide_combinations

def test_make_data_compatible():
    df = pd.DataFrame({
        'Name': ['Peptide1', 'Peptide2'],
        'MW': [500.1, 600.2],
        'iRT': [20.5, 25.7]
    })
    df.index = [10, 20]
    expected = [(10, 500.1, 20.5), (20, 600.2, 25.7)]
    result = make_data_compatible(df)
    assert result == expected, f"Expected {expected}, got {result}"

@pytest.mark.parametrize("pair, ppm_tolerance1, ppm_tolerance2, expected", [
    (((0, 500.0, 20.0), (1, 500.5, 20.2)), 1000, 0.5, True),
    (((0, 500.0, 20.0), (1, 510.0, 21.0)), 100, 1.0, False)
])
def test_within_ppm(pair, ppm_tolerance1, ppm_tolerance2, expected):
    assert within_ppm(pair, ppm_tolerance1, ppm_tolerance2) == expected

@pytest.mark.parametrize("pair, tolerance1, tolerance2, expected", [
    (((0, 500.0, 20.0), (1, 500.5, 20.2)), 1.0, 1.0, True),
    (((0, 500.0, 20.0), (1, 502.0, 22.0)), 1.0, 1.0, False)
])
def test_within_tolerance(pair, tolerance1, tolerance2, expected):
    assert within_tolerance(pair, tolerance1, tolerance2) == expected

def test_find_combinations_kdtree():
    data = [(0, 500.0, 20.0), (1, 500.5, 20.2), (2, 510.0, 21.0)]
    result = find_combinations_kdtree(data, 1.0, 1.0, use_ppm=False)
    assert len(result) == 1, "Expected one valid combination, but found none."

def test_process_peptide_combinations():
    df = pd.DataFrame({
        'Name': ['Peptide1', 'Peptide2', 'Peptide3'],
        'MW': [500.1, 500.5, 510.0],
        'iRT': [20.5, 20.8, 21.0]
    })
    df.index = [10, 20, 30]
    result_df = process_peptide_combinations(df, 1.0, 1.0, use_ppm=False)
    assert not result_df.empty, "Expected a non-empty DataFrame."
    assert 'index1' in result_df.columns, "Missing column 'index1'."
    assert 'peptide 1' in result_df.columns, "Missing column 'peptide 1'."