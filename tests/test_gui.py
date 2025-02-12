import pytest
import pandas as pd
import io
import logging
from unittest import mock
import streamlit 
from MSCI.gui.peptide_checker import (
    extract_peptide_and_charge,
    find_colliding_peptides,
    peptide_twins_checker
)

# Configure logging for debugging
logging.basicConfig(level=logging.DEBUG)


# ----------------------
# Test extract_peptide_and_charge
# ----------------------
def test_extract_peptide_and_charge():
    assert extract_peptide_and_charge("SDPYGIIR/2") == ("SDPYGIIR", 2), "Failed to extract peptide and charge."
    assert extract_peptide_and_charge("MKLVWR") == ("MKLVWR", None), "Failed to handle missing charge."


# ----------------------
# Test find_colliding_peptides
# ----------------------
def test_find_colliding_peptides():
    """Test the collision detection function with a sample dataset."""
    
    # Mocked input dataset
    data = """x_peptide,y_peptide
    SDPYGIIR/2,SDPYALVR/2
    MKWVTFISLL/2,SDPYGILR/2
    """
    
    # Read dataset into DataFrame
    df = pd.read_csv(io.StringIO(data))
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)  # Updated: Use .map() instead of deprecated .applymap()

    # Run function under test
    colliding = find_colliding_peptides(df, "SDPYGIIR", 2, "Reference Human Canonical proteome")

    # Debugging Output
    print("Colliding peptides:", colliding)  # Print to console
    logging.debug(f"Colliding peptides: {colliding}")  # Log output

    # Assertions to validate function behavior
    assert ("SDPYALVR", 2, None) in colliding, "Expected SDPYALVR/2 to collide with SDPYGIIR/2."
    assert ("AGTCTGAC", 3, None) not in colliding, "Unexpected collision detected."
