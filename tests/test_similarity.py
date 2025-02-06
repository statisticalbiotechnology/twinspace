import pandas as pd
import numpy as np
from MSCI.Similarity.spectral_angle_similarity import ndotproduct, nspectraangle, joinPeaks, process_spectra_pairs, process_spectra_pairs_cosine
from matchms import Spectrum

def test_ndotproduct():
    x = pd.DataFrame({'mz': [100.0, 200.0, 300.0], 'intensities': [0.1, 0.5, 0.9]})
    y = pd.DataFrame({'mz': [100.0, 200.0, 300.0], 'intensities': [0.2, 0.4, 0.8]})
    result = ndotproduct(x, y)
    assert isinstance(result, float), "ndotproduct should return a float"
    assert result > 0, "ndotproduct result should be positive"

def test_nspectraangle():
    x = pd.DataFrame({'mz': [100.0, 200.0, 300.0], 'intensities': [0.1, 0.5, 0.9]})
    y = pd.DataFrame({'mz': [100.0, 200.0, 300.0], 'intensities': [0.2, 0.4, 0.8]})
    result = nspectraangle(x, y)
    assert isinstance(result, float), "nspectraangle should return a float"
    assert 0 <= result <= 1, "nspectraangle should return a value between 0 and 1"

def test_joinPeaks():
    x = pd.DataFrame({'mz': [100.0, 200.0, 300.0], 'intensities': [0.1, 0.5, 0.9]})
    y = pd.DataFrame({'mz': [100.0, 201.0, 299.0], 'intensities': [0.2, 0.4, 0.8]})
    matcher = joinPeaks(tolerance=1)
    x_matched, y_matched = matcher.match(x, y)
    assert len(x_matched) == len(y_matched), "Matched spectra should have the same length"

def test_process_spectra_pairs():
    spectra = [Spectrum(mz=np.array([100.0, 200.0, 300.0], dtype=float), intensities=np.array([0.1, 0.5, 0.9], dtype=float)),
               Spectrum(mz=np.array([100.0, 201.0, 299.0], dtype=float), intensities=np.array([0.2, 0.4, 0.8], dtype=float))]
    mz_irt_df = pd.DataFrame({'Name': ['Peptide1', 'Peptide2'], 'MW': [500, 600], 'iRT': [20.5, 25.7]})
    chunk = [(0, 1)]
    result_df = process_spectra_pairs(chunk, spectra, mz_irt_df, tolerance=1)
    assert not result_df.empty, "Expected a non-empty DataFrame"
    assert 'similarity_score' in result_df.columns, "Missing similarity_score column"

def test_process_spectra_pairs_cosine():
    spectra = [Spectrum(mz=np.array([100.0, 200.0, 300.0], dtype=float), intensities=np.array([0.1, 0.5, 0.9], dtype=float)),
               Spectrum(mz=np.array([100.0, 201.0, 299.0], dtype=float), intensities=np.array([0.2, 0.4, 0.8], dtype=float))]
    mz_irt_df = pd.DataFrame({'Name': ['Peptide1', 'Peptide2'], 'MW': [500, 600], 'iRT': [20.5, 25.7]})
    chunk = [(0, 1)]
    result_df = process_spectra_pairs_cosine(chunk, spectra, mz_irt_df, tolerance=1)
    assert not result_df.empty, "Expected a non-empty DataFrame"
    assert 'similarity_score' in result_df.columns, "Missing similarity_score column"
