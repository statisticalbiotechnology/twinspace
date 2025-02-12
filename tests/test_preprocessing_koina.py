import pytest
import pandas as pd
import requests
from MSCI.Preprocessing.Koina import PeptideProcessor, MASSES  

@pytest.fixture
def peptide_processor():
    return PeptideProcessor(
        input_file="test_input.csv",
        collision_energy=30,
        charge=2,
        model_intensity="Prosit_2020",
        model_irt="iRT_Predictor"
    )

def test_calculate_peptide_mass(peptide_processor):
    peptide = "ACDE"
    expected_mass = (71.037114 + 103.009185 + 115.026943 + 129.042593
                     + MASSES["N_TERMINUS"]  # Fix: Use global MASSES
                     + MASSES["C_TERMINUS"])
    
    calculated_mass = peptide_processor.calculate_peptide_mass(peptide)
    
    assert calculated_mass is not None, "calculate_peptide_mass returned None"
    assert pytest.approx(calculated_mass, rel=1e-6) == expected_mass

def test_get_request_body(peptide_processor):
    peptides = ["ACDE", "WXYZ"]
    model_type = "Prosit_2020"
    request_body = peptide_processor.get_request_body(peptides, model_type)
    
    assert isinstance(request_body, dict)
    assert "inputs" in request_body
    assert len(request_body["inputs"]) > 0
    assert request_body["inputs"][0]["name"] == "peptide_sequences"

def test_format_msp(peptide_processor):
    peptide = "ACDE"
    charge = 2
    collision_energy = 30
    mz_values = [100.1, 200.2, 300.3]
    intensities = [1000.0, 500.0, 100.0]
    irt = 12.34
    
    msp_entry = peptide_processor.format_msp(peptide, charge, collision_energy, mz_values, intensities, irt)
    
    assert "Name:" in msp_entry
    assert "MW:" in msp_entry
    assert "Collision_energy:" in msp_entry
    assert "iRT:" in msp_entry
    assert "Num peaks:" in msp_entry
    assert "100.10" in msp_entry
    assert "200.20" in msp_entry
    assert "300.30" in msp_entry

def test_get_predictions(monkeypatch, peptide_processor):
    def mock_post(url, json):
        class MockResponse:
            status_code = 200
            def json(self):
                return {
                    "outputs": [
                        {},  # Placeholder for response structure
                        {"data": [100.1, 200.2, 300.3]},  
                        {"data": [1000.0, 500.0, 100.0]}                      ]
                }
        return MockResponse()

    monkeypatch.setattr(requests, "post", mock_post)  # Fix: Correctly mock requests.post
    
    peptides = ["ACDE"]
    df = peptide_processor.get_predictions(peptides)
    
    assert df is not None, "get_predictions returned None"
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "mz_values" in df.columns
    assert "intensities" in df.columns
    assert len(df["mz_values"][0]) == 3
    assert len(df["intensities"][0]) == 3

def test_get_irt_predictions(monkeypatch, peptide_processor):
    def mock_post(url, json):
        class MockResponse:
            status_code = 200
            def json(self):
                return {"outputs": [{"data": [12.34]}]}  # Ensure output format matches expectations
        return MockResponse()

    monkeypatch.setattr(requests, "post", mock_post)
    
    peptides = ["ACDE"]
    irt_values = peptide_processor.get_irt_predictions(peptides)
    
    assert isinstance(irt_values, list)
    assert len(irt_values) == 1
    assert pytest.approx(irt_values[0], rel=1e-6) == 12.34
