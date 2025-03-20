import MSCI
from MSCI.Preprocessing.Koina import PeptideProcessor
from MSCI.Grouping_MS1.Grouping_mw_irt import process_peptide_combinations
from MSCI.Preprocessing.read_msp_file import read_msp_file
from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs
from MSCI.data.digest import parse_fasta_and_digest, tryptic_digest, peptides_to_csv
from matchms.importing import load_from_msp
import random
import numpy as np
import pandas as pd
import requests


# Parameters
fasta_file = "https://raw.githubusercontent.com/proteomicsunitcrg/MSCI/main/tutorial/sp_human_2023_04.fasta"
collision_energy = 30
charge = 2
model_intensity = "Prosit_2020_intensity_HCD"
model_irt = "Prosit_2019_irt"


# Function to read and digest a fasta file
def read_fasta_file(fasta_file):
    # Read fasta file
    fasta_file = open(fasta_file, "r")
    fasta_file = fasta_file.read()
    # Digest fasta file
    peptides = parse_fasta_and_digest(fasta_file, enzyme="trypsin", missed_cleavages=0)
    return peptides


## Main script
peptides = read_fasta_file(fasta_file)
peptides_to_csv(peptides, "peptides.csv")

# filter peptides for length
for key in result.keys():
    result[key] = list(set(filter(lambda pep: 6 <= len(pep) <= 60, result[key])))

# Predict spectra
processor = PeptideProcessor(
    input_file="not_needed",
    collision_energy=collision_energy,
    charge=charge,
    model_intensity=model_intensity,
    model_irt=model_irt,
)

batch_size = 1000
total_batches = len(peptides) // batch_size + (
    1 if len(peptides) % batch_size != 0 else 0
)
predictions = []
for i, start in enumerate(range(0, len(peptides), batch_size)):
    batch_peptides = peptides[start : start + batch_size]
    irt_values = processor.get_irt_predictions(batch_peptides)
    df = processor.get_predictions(batch_peptides)

    if df is not None and irt_values is not None:
        df["irt"] = irt_values
        predictions.append(df)
    print(f"Progress: {i + 1}/{total_batches}")
