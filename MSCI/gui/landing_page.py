# landing_page.py

import streamlit as st
from pathlib import Path
import sys
sys.path.append('Z:/zelhamraoui/MSCA_Package/MSCI_package/MSCI')

from MSCI.gui.utils import load_image

def landing_page():
    """Render the landing page of the application."""
    st.title("Welcome to MSCI")

    # URL of the image on GitHub
    workflow_url = "https://github.com/proteomicsunitcrg/MSCI/raw/main/docs/INTRODUCTION.png"

    st.write("""
    Peptide identification by mass spectrometry relies on the interpretation of fragmentation spectra based on the m/z pattern, relative intensities, and retention time (RT). Given a proteome, we wondered how many peptides generate very similar fragmentation spectra with current MS methods. MSCI is a Python package built to assess the information content of peptide fragmentation spectra. We aimed to calculate an information-content index for all peptides in a given proteome, which would enable us to design data acquisition and data analysis strategies that generate and prioritize the most informative fragment ions to be queried for peptide quantification.
    """)

    # Load image from URL
    workflow_image = load_image(workflow_url)

    if workflow_image:
        st.markdown(f"""
        <p align="center">
            <img src="data:image/png;base64,{workflow_image}" alt="workflow illustration" width="1200">
        </p>
        """, unsafe_allow_html=True)
    else:
        st.error(f"Workflow image not found at {workflow_url}")

    st.subheader("Installation")
    st.write("""
    **Prerequisites:**

    - Python 3.8 - 3.11
    - Matchms
    """)


    st.subheader("Implementation and Example")
    st.write("""
    **Open the Notebook**: Click on the following [link to Google Colab](https://colab.research.google.com/drive/1ny97RNgvnpD7ZrHW8TTRXWCAQvIcavkk?usp=sharing)
    """)

    st.subheader("Contribution")
    st.write("""
    If you would like to contribute to this project, feel free to fork the repository on GitHub and submit a pull request.
    """)
