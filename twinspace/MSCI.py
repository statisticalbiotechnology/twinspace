import numpy as np
import pandas as pd
import random

from MSCI.Preprocessing.Koina import PeptideProcessor
from matchms.importing import load_from_msp
from MSCI.Preprocessing.read_msp_file import read_msp_file
from MSCI.Grouping_MS1.Grouping_mw_irt import process_peptide_combinations
from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs, process_spectra_pairs_parallel


# Parse fasta file
import requests
# URL of the FASTA file
url = "https://raw.githubusercontent.com/proteomicsunitcrg/MSCI/main/tutorial/sp_human_2023_04.fasta"

# Download the file
response = requests.get(url)
if response.status_code == 200:
    fasta_file = "sp_human_2023_04.fasta"  # Save as a local file
    with open(fasta_file, "w") as f:
        f.write(response.text)
else:
    raise Exception(f"Failed to download file: {response.status_code}")

# Now pass the local file path to your function
result = parse_fasta_and_digest(fasta_file, digest_type="trypsin")

# Save output
peptides_to_csv(result, "random_tryptic_peptides.txt")

# Predict with Koina
processor = PeptideProcessor(
    input_file="random_tryptic_peptides.txt",
    collision_energy=30,
    charge=2,
    model_intensity="Prosit_2020_intensity_HCD",
    model_irt="Prosit_2019_irt"
)

processor.process('random_tryptic_peptides.msp')

# Load dataset
File= 'random_tryptic_peptides.msp'
spectra = list(load_from_msp(File))

# Group within MS1 tolerance
mz_tolerance = 1
irt_tolerance = 5


mz_irt_df = read_msp_file(File)
Groups_df = process_peptide_combinations(mz_irt_df, mz_tolerance, irt_tolerance, use_ppm=False)

Groups_df

# Calculate the similarity within fragment tolerance
Groups_df.columns = Groups_df.columns.str.strip()
index_array = Groups_df[['index1','index2']].values.astype(int)
#result = process_spectra_pairs(index_array, spectra,  mz_irt_df, tolerance =0, ppm=10)
result = process_spectra_pairs_parallel(index_array, spectra,  mz_irt_df, tolerance =0, ppm=10)
result.to_csv("output.csv", index=False)
result

# Plot spectra of interest using matchms
import matplotlib.pyplot as plt
print(mz_irt_df.iloc[19])
print(mz_irt_df.iloc[36])
spectra[19].plot_against(spectra[36])
plt.savefig('spectra_comparison.png')
