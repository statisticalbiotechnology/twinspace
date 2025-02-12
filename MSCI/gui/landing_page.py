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
    The MSCI Python library was developed to address the challenges of peptide identification in mass spectrometry-based proteomics, particularly regarding the issue of indistinguishable peptides that exhibit similar analytical values and fragmentation patterns. MSCI provides a comprehensive toolset to streamline the workflow from data import to spectral analysis, enabling researchers to effectively evaluate fragmentation similarity scores, identify indistinguishable peptide pairs, and design data acquisition and analysis strategies that prioritize the most informative fragment ions for accurate peptide quantification.
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
    **Open the Notebook**: Click on the following [link to Google Colab](https://colab.research.google.com/drive/1ny97RNgvnpD7ZrHW8TTRXWCAQvIcavkk)
    """)

    st.subheader("Contribution")
    st.write("""
    If you would like to contribute to this project, feel free to fork the repository on GitHub and submit a pull request.
    """)

    st.subheader("Please Cite")
    st.write("""
    If you use the MSCI package for your research, please cite the following work:
    
    MSCI: an open-source Python package for information content assessment of peptide fragmentation spectra, Zahra Elhamraoui, Eva Borràs, Mathias Wilhelm, Eduard Sabidó
    
    """)

if __name__ == "__main__":
    landing_page()
