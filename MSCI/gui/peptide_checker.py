import streamlit as st
import pandas as pd
from typing import Dict, List, Tuple
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
        30: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/mutation_peptides_HRHR_NSA_score.csv",
    },
    "Immunopeptidome": {
        30: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/HLA_peptides_HRHR_NSA_score.csv",
    },
    "Human Oral microbiome": {
        30: "https://raw.githubusercontent.com/zahrael97/MSCI/master/Database/Oral_microbiom.csv",
    }
}

def extract_peptide_and_charge(peptide_str: str) -> Tuple[str, int]:
    try:
        peptide, charge = peptide_str.rsplit('/', 1)
        return peptide, int(charge)
    except ValueError:
        return peptide_str, None

def find_colliding_peptides(df: pd.DataFrame, peptide: str, charge: int, organism: str) -> set:
    angle_available = organism in ["Reference Human Canonical proteome with natural variants", "Immunopeptidome"]

    matches = df[((df['x_peptide'].apply(lambda x: extract_peptide_and_charge(x) == (peptide, charge))) |
                  (df['y_peptide'].apply(lambda x: extract_peptide_and_charge(x) == (peptide, charge))))]

    colliding_peptides = {
        (row['x_peptide'].rsplit('/', 1)[0], int(row['x_peptide'].rsplit('/', 1)[1]), row['angle'] if angle_available else None)
        for _, row in matches.iterrows() if row['x_peptide'].rsplit('/', 1)[0] != peptide
    }.union({
        (row['y_peptide'].rsplit('/', 1)[0], int(row['y_peptide'].rsplit('/', 1)[1]), row['angle'] if angle_available else None)
        for _, row in matches.iterrows() if row['y_peptide'].rsplit('/', 1)[0] != peptide
    })

    return colliding_peptides

def peptide_twins_checker():
    st.header("Peptide Twins Checker")
    st.markdown("""
    This tool checks whether a provided peptide is indistinguishable from any other peptide
    within a pre-calculated search space (e.g., human canonical proteome, human immunopeptidome).
    """)
    
    organism = st.selectbox("Select Universe:", options=list(DATASETS.keys()), key='Universe')
    energies = st.multiselect("Select Collision Energies:", options=list(DATASETS[organism].keys()), key='energies')

    peptide = st.text_input("Enter Peptide:", key='peptide', value="SDPYGIIR")
    
    if peptide != peptide.upper():
        st.warning(f"Peptide sequence converted to uppercase: {peptide.upper()}")
        peptide = peptide.upper()

    charge = st.number_input("Enter Charge:", min_value=1, step=1, value=2, key='charge')

    if st.button("Check twins"):
        if peptide:
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
                    colliding_peptides = find_colliding_peptides(df, peptide, charge, organism)
                    if colliding_peptides:
                        colliding_info[energy] = colliding_peptides

            if colliding_info:
                st.info(f"Peptides that are similar to {peptide} in charge {charge}:")
                for energy, peptides in colliding_info.items():
                    st.write(f"**NCE {energy}:**")
                    for colliding_peptide, colliding_charge, angle in peptides:
                        angle_text = f" (Spectral Angle: {angle:.4f})" if angle is not None else ""
                        st.write(f"- {colliding_peptide} in charge {colliding_charge}{angle_text}")
            else:
                st.success("No twin peptides detected in the selected energies.")
        else:
            st.warning("Please enter a peptide and charge to check.")

if __name__ == "__main__":
    peptide_twins_checker()
