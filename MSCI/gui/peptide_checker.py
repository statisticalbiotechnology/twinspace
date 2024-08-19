# peptide_checker.py

import streamlit as st
import pandas as pd
from typing import Dict, List, Tuple
from Bio import SeqIO
import io

DATASETS = {
    "Human": {
        25: "Z:/zelhamraoui/MSCA_Package/results/Tryptic/01_08_2024/NSA_HRHR_NCE25.csv",
        28: "Z:/zelhamraoui/MSCA_Package/results/Tryptic/01_08_2024/NSA_HRHR_NCE28.csv",
        30: "Z:/zelhamraoui/MSCA_Package/results/Tryptic/01_08_2024/NSA_HRHR_NCE30.csv",
        32: "Z:/zelhamraoui/MSCA_Package/results/Tryptic/01_08_2024/NSA_HRHR_NCE32.csv",
        35: "Z:/zelhamraoui/MSCA_Package/results/Tryptic/01_08_2024/NSA_HRHR_NCE35.csv"
    },
    "Immunopeptidome": {
        28: "Z:/zelhamraoui/MSCA_Package/results/Tryptic/01_08_2024/NSA_HRHR_NCE28.csv",
    },
    "Mutated human proteome": {
        28: "Z:/zelhamraoui/MSCA_Package/results/Tryptic/01_08_2024/NSA_HRHR_NCE28.csv",
    }
}

def parse_fasta(fasta_file) -> Dict[str, str]:
    """Parse a FASTA file and return a dictionary of protein sequences."""
    protein_sequences = {}
    fasta_text = fasta_file.read().decode('utf-8')
    fasta_stream = io.StringIO(fasta_text)

    for record in SeqIO.parse(fasta_stream, "fasta"):
        protein_sequences[record.id] = str(record.seq)

    return protein_sequences

def extract_peptide_and_charge(peptide_str: str) -> Tuple[str, int]:
    """Extract peptide sequence and charge from a string."""
    try:
        peptide, charge = peptide_str.rsplit('/', 1)
        return peptide, int(charge)
    except ValueError:
        return peptide_str, None

def find_peptide_in_proteins(peptide: str, protein_sequences: Dict[str, str]) -> List[str]:
    """Find proteins containing a given peptide sequence."""
    return [protein_id for protein_id, sequence in protein_sequences.items() if peptide in sequence]

def find_colliding_peptides(df: pd.DataFrame, peptide: str, charge: int) -> set:
    """Find peptides that collide with the given peptide."""
    if df is None:
        st.error("The provided DataFrame is None. Please ensure data is loaded correctly.")
        return set()

    matches = df[((df['x_peptide'].apply(lambda x: extract_peptide_and_charge(x) == (peptide, charge))) |
                  (df['y_peptide'].apply(lambda x: extract_peptide_and_charge(x) == (peptide, charge)))) &
                 (df['angle'] > 0.7)]

    colliding_peptides = {
        (row['x_peptide'].rsplit('/', 1)[0], int(row['x_peptide'].rsplit('/', 1)[1]))
        for _, row in matches.iterrows() if row['x_peptide'].rsplit('/', 1)[0] != peptide
    }.union({
        (row['y_peptide'].rsplit('/', 1)[0], int(row['y_peptide'].rsplit('/', 1)[1]))
        for _, row in matches.iterrows() if row['y_peptide'].rsplit('/', 1)[0] != peptide
    })

    return colliding_peptides

def peptide_twins_checker():
    """Render the Peptide Twins Checker page."""
    st.header("Peptide Twins Checker")

    organism = st.selectbox("Select Organism:", options=list(DATASETS.keys()), key='organism')
    energies = st.multiselect("Select Collision Energies:", options=list(DATASETS[organism].keys()), key='energies')

    peptide = st.text_input("Enter Peptide:", key='peptide', value="SDPYGIIR")
    charge = st.number_input("Enter Charge:", min_value=1, step=1, value=2, key='charge')

    fasta_file = st.file_uploader("Upload FASTA file with protein sequences", type=["fasta"])

    if st.button("Check twins"):
        if peptide and fasta_file:
            with st.spinner("Parsing FASTA file..."):
                protein_sequences = parse_fasta(fasta_file)

            colliding_info = {}

            for energy in energies:
                df_path = DATASETS[organism][energy]
                with st.spinner(f"Loading data for NCE {energy}..."):
                    df = pd.read_csv(df_path)

                if df.empty:
                    st.warning(f"No data found in the file for NCE {energy}. Skipping.")
                    continue

                with st.spinner(f"Checking for twins in NCE {energy}..."):
                    colliding_peptides = find_colliding_peptides(df, peptide, charge)
                    if colliding_peptides:
                        colliding_info[energy] = colliding_peptides

            if colliding_info:
                st.info(f"Peptides that are similar to {peptide} in charge {charge}:")
                for energy, peptides in colliding_info.items():
                    st.write(f"**NCE {energy}:**")
                    for colliding_peptide, colliding_charge in peptides:
                        st.write(f"- {colliding_peptide} in charge {colliding_charge}")
                        matching_proteins = find_peptide_in_proteins(colliding_peptide, protein_sequences)
                        if matching_proteins:
                            st.write(f"  Found in proteins: {', '.join(matching_proteins)}")
                        else:
                            st.write("  No proteins found for this peptide.")
            else:
                st.success("No twin peptides detected in the selected energies.")
        else:
            st.warning("Please enter a peptide, charge, and upload a FASTA file to check.")

