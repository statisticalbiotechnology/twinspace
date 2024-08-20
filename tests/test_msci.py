import pytest
import sys
import os
sys.path.append('/home/zahra/Downloads/MSCI')
from MSCI.Preprocessing.Koina import PeptideProcessor
from MSCI.Preprocessing.read_msp_file import read_msp_file
from MSCI.Grouping_MS1.Grouping_mw_irt import process_peptide_combinations
from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs
from matchms.importing import load_from_msp
import pandas as pd
# Ensure pandas has the version attribute
if not hasattr(pd, 'version'):
    class Version:
        def __init__(self, version):
            self.version = version
    pd.version = Version(pd.__version__)
@pytest.fixture
def peptide_processor():
    """Fixture to create a PeptideProcessor instance."""
    processor = PeptideProcessor(
        input_file="/home/zahra/Downloads/random_tryptic_peptides.txt",
        collision_energy=30,
        charge=2,
        model_intensity="Prosit_2020_intensity_HCD",
        model_irt="Prosit_2019_irt"
    )
    return processor

def test_peptide_processor(peptide_processor):
    """Test the PeptideProcessor's processing."""
    output_file = 'output.msp'
    peptide_processor.process(output_file)
    
    assert os.path.exists(output_file), "The output file should be created"

    # Additional assertions to verify file content
    with open(output_file, 'r') as file:
        content = file.read()
        assert 'BEGIN IONS' in content, "Output file should contain expected MSP format content"


    # Additional assertions could be added here to check the content of 'result'

def test_read_msp_file():
    """Test reading the .msp file."""
    mz_irt_df = read_msp_file('output.msp')
    assert not mz_irt_df.empty, "The DataFrame should not be empty"
    assert 'Name' in mz_irt_df.columns, "The 'Name' column should be present in the DataFrame"
    assert 'MW' in mz_irt_df.columns, "The 'MW' column should be present in the DataFrame"
    assert 'iRT' in mz_irt_df.columns, "The 'iRT' column should be present in the DataFrame"

@pytest.mark.parametrize("mz_tolerance, irt_tolerance", [(1, 10), (0.5, 5), (2, 20)])
def test_process_peptide_combinations_param(mz_tolerance, irt_tolerance):
    mz_irt_df = read_msp_file('output.msp')
    Groups_df = process_peptide_combinations(mz_irt_df, mz_tolerance, irt_tolerance, use_ppm=False)
    assert not Groups_df.empty, f"The Groups DataFrame should not be empty for mz_tolerance={mz_tolerance} and irt_tolerance={irt_tolerance}"

def test_process_spectra_pairs():
    """Test processing and comparing spectra pairs."""
    spectra = list(load_from_msp('output.msp'))
    mz_irt_df = read_msp_file('output.msp')
    mz_tolerance = 1
    irt_tolerance = 10
    Groups_df = process_peptide_combinations(mz_irt_df, mz_tolerance, irt_tolerance, use_ppm=False)    
    Groups_df.columns = Groups_df.columns.str.strip()
    mz_irt_df = read_msp_file('output.msp')
    index_array = Groups_df[['index1','index2']].values.astype(int)
    result = process_spectra_pairs(index_array, spectra, mz_irt_df, tolerance=0, ppm=10)
    assert result is not None, "The function should return a result"
    assert 'similarity_score' in result.columns, "The result should contain 'similarity_score' column"