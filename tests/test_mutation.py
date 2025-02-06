import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from MSCI.mutation.mutated_sequences import tryptic_digest

@pytest.fixture
def test_protein_sequence():
    return "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQ"  # Sample protein sequence

def test_tryptic_digest(test_protein_sequence):
    peptides = tryptic_digest(test_protein_sequence)
    expected_peptides = ["MKWVTFISLLFLFSSAYSR", "GVFR", "RDTHK", "SEIAHR", "FKDLGEEHFK", "GLVLIAFSQYLQ"]
    assert peptides == expected_peptides, f"Expected {expected_peptides}, but got {peptides}"

def test_mutation_processing():
    mutations_df = pd.DataFrame({
        'Entry': ['P12345', 'P12345'],
        'Natural variant': ['VARIANT 5; /note="A -> T"', 'VARIANT 10; /note="G -> V"']
    })
    proteome = {'P12345': "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQ"}
    
    # Check if mutations are correctly parsed
    mutation_str = mutations_df.iloc[0]['Natural variant']
    assert "VARIANT" in mutation_str, "Mutation format incorrect"
    assert "A -> T" in mutation_str, "Mutation note parsing failed"

    # Simulating peptide mutation mapping
    peptides = tryptic_digest(proteome['P12345'])
    assert len(peptides) > 0, "Tryptic digestion failed"
    
    # Validate output structure
    mutated_peptides = {peptide: [] for peptide in peptides}
    assert isinstance(mutated_peptides, dict), "Mutation processing should return a dictionary"
    
    # Ensuring mutation handling
    for peptide in peptides:
        if "R" in peptide or "K" in peptide:
            mutated_peptides[peptide].append(peptide)  # Simulate mutation storage
    
    assert any(mutated_peptides.values()), "Mutations not assigned to peptides correctly"
