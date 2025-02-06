import pytest
import pandas as pd
from MSCI.Preprocessing.Koina import PeptideProcessor

@pytest.fixture
def peptide_processor():
    return PeptideProcessor("test.csv", collision_energy=30, charge=2, model_intensity="Prosit_2020", model_irt="iRT")

@pytest.mark.parametrize("peptide, expected_mass", [
    ("ACDE", 436.12639967),  # Updated expected mass including N-terminus and C-terminus
    ("M[UNIMOD:35]K", 293.14092767),
    ("K[UNIMOD:737]R", 531.36957067),
])
def test_calculate_peptide_mass(peptide_processor, peptide, expected_mass):
    result = peptide_processor.calculate_peptide_mass(peptide)
    assert pytest.approx(result, rel=1e-5) == expected_mass, f"Expected {expected_mass}, got {result}"


def test_get_request_body(peptide_processor):
    peptides = ["ACDE", "KLMN"]
    body = peptide_processor.get_request_body(peptides, "Prosit_2020")
    assert "inputs" in body, "Request body is missing 'inputs' key"
    assert len(body["inputs"][0]["data"]) == len(peptides), "Mismatch in peptide sequences"


def test_format_msp(peptide_processor):
    formatted_msp = peptide_processor.format_msp(
        peptide="ACDE",
        charge=2,
        collision_energy=30,
        mz_values=[100.1, 200.2, 300.3],
        intensities=[0.1, 0.5, 0.9],
        irt=45.6
    )
    assert "Name: ACDE/2" in formatted_msp, "MSP format missing expected name header"
    assert "MW:" in formatted_msp, "Missing MW field"
    assert "Collision_energy: 30.00" in formatted_msp, "Incorrect collision energy format"