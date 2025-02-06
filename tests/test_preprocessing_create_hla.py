import pytest
from MSCI.Preprocessing.create_hla import generate_variable_length_peptides, extract_peptides_from_fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

def test_generate_variable_length_peptides():
    protein_sequence = "ABCDEFGHIJKL"
    peptides = generate_variable_length_peptides(protein_sequence, 8, 11)
    
    expected_peptides = [
        "ABCDEFGH", "BCDEFGHI", "CDEFGHIJ", "DEFGHIJK", "EFGHIJKL", 
        "ABCDEFGHI", "BCDEFGHIJ", "CDEFGHIJK", "DEFGHIJKL", 
        "ABCDEFGHIJ", "BCDEFGHIJK", "CDEFGHIJKL", 
        "ABCDEFGHIJK", "BCDEFGHIJKL"
    ]
    
    assert set(peptides) == set(expected_peptides), "Peptides generation failed"


def test_extract_peptides_from_fasta(tmp_path):
    fasta_path = tmp_path / "test.fasta"
    
    # Create a mock FASTA file
    records = [
        SeqRecord(Seq("ABCDEFGHIJKL"), id="test1", description=""),
        SeqRecord(Seq("MNOPQRSTUVWXYZ"), id="test2", description="")
    ]
    
    with open(fasta_path, "w") as f:
        SeqIO.write(records, f, "fasta")
    
    peptides = extract_peptides_from_fasta(str(fasta_path), 8, 11)
    
    expected_peptides = set(generate_variable_length_peptides("ABCDEFGHIJKL", 8, 11) +
                            generate_variable_length_peptides("MNOPQRSTUVWXYZ", 8, 11))
    
    assert set(peptides) == expected_peptides, "Extracted peptides do not match expected values"
