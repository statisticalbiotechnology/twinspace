import streamlit as st
import pandas as pd
from typing import Dict, List, Tuple
from Bio import SeqIO
import io
import requests

DATASETS = {
    "Reference Human Canonical proteome": {
        25: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/NSA_HRHR_NCE25.csv",
        28: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/NSA_HRHR_NCE28.csv",
        30: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/NSA_HRHR_NCE30.csv",
        32: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/NSA_HRHR_NCE32.csv",
        35: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/NSA_HRHR_NCE35.csv"
    },
    "Reference Human Canonical proteome with natural variants": {
        28: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/mutation_peptides_HRHR_NSA_score.csv",
    },

    "Immunopeptidome": {
        28: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/HLA_peptides_HRHR_NSA_score.csv",
    },
    "Human Oral microbiome": {
        28: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/Oral_microbiom.csv",
    }
}

FASTA_URL = "https://raw.githubusercontent.com/proteomicsunitcrg/MSCI/refs/heads/main/tutorial/sp_human_2023_04.fasta"

def parse_fasta(fasta_content: str) -> Dict[str, str]:
    """Parse a FASTA content and return a dictionary of protein sequences."""
    protein_sequences = {}
    fasta_stream = io.StringIO(fasta_content)

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
    matches = df[((df['x_peptide'].apply(lambda x: extract_peptide_and_charge(x) == (peptide, charge))) |
                  (df['y_peptide'].apply(lambda x: extract_peptide_and_charge(x) == (peptide, charge))))]

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
    st.markdown("""
    In the Peptide Twin Checker tool users can enter a peptide of interest whether the provided peptide is indistinguishable from any other peptide within a pre-calculated search space (e.g., human canonical proteome, human immunopeptidome)    """)
    
    organism = st.selectbox("Select Universe:", options=list(DATASETS.keys()), key='Universe')
    energies = st.multiselect("Select Collision Energies:", options=list(DATASETS[organism].keys()), key='energies')

    peptide = st.text_input("Enter Peptide:", key='peptide', value="SDPYGIIR")

    # Automatically convert to uppercase if the user enters lowercase letters
    if peptide != peptide.upper():
        st.warning(f"Peptide sequence converted to uppercase: `{peptide.upper()}`")
        peptide = peptide.upper()

    charge = st.number_input("Enter Charge:", min_value=1, step=1, value=2, key='charge')
    
    fasta_option = st.radio("Upload FASTA for peptide-protein annotation", ("Upload File", "Use Default (Human Proteome)"))
    fasta_file = None
    fasta_content = ""
    st.markdown("""The uploaded FASTA file is used to identify the proteins that contain the colliding peptides found in the selected dataset     """)   
    if fasta_option == "Upload File":
        fasta_file = st.file_uploader("Upload FASTA file with protein sequences", type=["fasta"])
    else:
        with st.spinner("Fetching default FASTA file..."):
            response = requests.get(FASTA_URL)
            if response.status_code == 200:
                fasta_content = response.text
            else:
                st.error("Failed to load default FASTA file.")
                return

    if st.button("Check twins"):
        if peptide and (fasta_file or fasta_content):
            with st.spinner("Parsing FASTA file..."):
                if fasta_file:
                    fasta_content = fasta_file.read().decode('utf-8')
                protein_sequences = parse_fasta(fasta_content)

            colliding_info = {}

            for energy in energies:
                df_path = DATASETS[organism][energy]
                with st.spinner(f"Loading data for NCE {energy}..."):
                    try:
                        df = pd.read_csv(df_path, delimiter=',')
                    except Exception as e:
                        st.error(f"An error occurred while loading the CSV file: {e}")
                        continue

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
            st.warning("Please enter a peptide, charge, and ensure a FASTA file is available to check.")

if __name__ == "__main__":
    peptide_twins_checker()
