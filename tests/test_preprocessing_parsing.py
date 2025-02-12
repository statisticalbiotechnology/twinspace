import pytest
import os
from MSCI.Preprocessing.Parsing import read_msp_file, read_mgf_file, read_mzml_file

# Set the correct path to the Example_data directory
data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "Example_data")

@pytest.fixture
def mgf_test_file():
    return os.path.join(data_dir, "example.mgf")

@pytest.fixture
def mzml_test_file():
    return os.path.join(data_dir, "example.mzML")

def test_read_mgf_file(mgf_test_file):
    assert os.path.exists(mgf_test_file), f"File not found: {mgf_test_file}"  # Debugging
    spectra = read_mgf_file(mgf_test_file)
    assert isinstance(spectra, list), "MGF file should return a list"
    assert len(spectra) > 0, "No spectra found in MGF file"
    assert "mz_values" in spectra[0], "Missing mz_values in parsed spectrum"
    assert "intensities" in spectra[0], "Missing intensities in parsed spectrum"

def test_read_mzml_file(mzml_test_file):
    assert os.path.exists(mzml_test_file), f"File not found: {mzml_test_file}"  # Debugging
    spectra = read_mzml_file(mzml_test_file)
    assert isinstance(spectra, list), "MZML file should return a list"
    assert len(spectra) > 0, "No spectra found in MZML file"
