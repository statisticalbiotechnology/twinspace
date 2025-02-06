import pytest
import pandas as pd
import io
from unittest import mock
import streamlit as st
from MSCI.gui.peptide_checker import (
    parse_fasta,
    extract_peptide_and_charge,
    find_peptide_in_proteins,
    find_colliding_peptides,
    peptide_twins_checker
)

# ----------------------
# Test parse_fasta
# ----------------------
def test_parse_fasta():
    fasta_content = ">Protein1\nMKWVTFISLLFLFSSAYSR\n>Protein2\nAGTCTGACGTAGCTAGTCG\n"
    fasta_file = io.StringIO(fasta_content)
    protein_sequences = parse_fasta(fasta_file)
    assert isinstance(protein_sequences, dict), "Parsed FASTA should return a dictionary."
    assert "Protein1" in protein_sequences, "Protein1 is missing in parsed sequences."
    assert protein_sequences["Protein1"] == "MKWVTFISLLFLFSSAYSR", "Protein1 sequence mismatch."

# ----------------------
# Test extract_peptide_and_charge
# ----------------------
def test_extract_peptide_and_charge():
    assert extract_peptide_and_charge("SDPYGIIR/2") == ("SDPYGIIR", 2), "Failed to extract peptide and charge."
    assert extract_peptide_and_charge("MKLVWR") == ("MKLVWR", None), "Failed to handle missing charge."

# ----------------------
# Test find_peptide_in_proteins
# ----------------------
def test_find_peptide_in_proteins():
    protein_sequences = {"Protein1": "MKWVTFISLLFLFSSAYSR", "Protein2": "AGTCTGACGTAGCTAGTCG"}
    matches = find_peptide_in_proteins("FISLLF", protein_sequences)
    assert "Protein1" in matches, "Peptide should be found in Protein1."
    assert "Protein2" not in matches, "Peptide should not be found in Protein2."

# ----------------------
# Test find_colliding_peptides
# ----------------------
def test_find_colliding_peptides():
    data = """x_peptide,y_peptide
    SDPYGIIR/2,AGTCTGAC/3
    MKWVTFISLL/2,SDPYGIIR/2
    """
    df = pd.read_csv(io.StringIO(data))
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)  # Strip extra spaces
    colliding = find_colliding_peptides(df, "SDPYGIIR", 2)
    assert ("MKWVTFISLL", 2) in colliding, "Expected MKWVTFISLL/2 to collide with SDPYGIIR/2."
    assert ("AGTCTGAC", 3) not in colliding, "Unexpected collision detected."

# ----------------------
# Mock Streamlit GUI Components for Testing peptide_twins_checker
# ----------------------
@mock.patch("streamlit.text_input", return_value="SDPYGIIR")
@mock.patch("streamlit.number_input", return_value=2)
@mock.patch("streamlit.selectbox", return_value="Reference Human Canonical proteome")
@mock.patch("streamlit.multiselect", return_value=[28])
@mock.patch("streamlit.file_uploader", return_value=None)  # Simulate no file uploaded
@mock.patch("streamlit.button", return_value=True)  # Simulate button click
def test_peptide_twins_checker(
    mock_text_input,
    mock_number_input,
    mock_selectbox,
    mock_multiselect,
    mock_file_uploader,
    mock_button
):
    """Test that the Streamlit UI loads correctly and responds to input."""
    
    with mock.patch("streamlit.write") as mock_write, mock.patch("streamlit.warning") as mock_warning:
        peptide_twins_checker()
        
        # Verify UI elements were called
        mock_text_input.assert_called_with("Enter Peptide:", key="peptide", value="SDPYGIIR")
        mock_number_input.assert_called_with("Enter Charge:", min_value=1, step=1, value=2, key="charge")
        mock_selectbox.assert_called_with("Select Universe:", options=list(DATASETS.keys()), key="Universe")
        mock_multiselect.assert_called()
        mock_button.assert_called()

        # Check if warnings are triggered when no file is uploaded
        mock_warning.assert_called_with("Please enter a peptide, charge, and upload a FASTA file to check.")
