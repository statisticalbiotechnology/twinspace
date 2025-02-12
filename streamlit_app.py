import streamlit as st
from pathlib import Path
import base64
import requests
from MSCI.gui.utils import load_image_from_url, add_custom_css
from MSCI.gui.landing_page import landing_page
from MSCI.gui.peptide_analysis import peptide_twins_analysis, plot_spectra
from MSCI.gui.peptide_checker import peptide_twins_checker

def main():
    st.set_page_config(layout="wide")
    add_custom_css()

    # URL of the image on GitHub
    logo_url = "https://github.com/proteomicsunitcrg/MSCI/raw/main/docs/MSCI_logor.png"
    logo_image = load_image_from_url(logo_url)

    with st.sidebar:
        if logo_image:
            st.markdown(f"""
            <p align="center">
                <img src="data:image/png;base64,{logo_image}" alt="logo" width="300" height="300">
            </p>
            """, unsafe_allow_html=True)

        # Set the default option to an empty string
        option = st.radio("Choose an option", ("HOME", "Peptide Twins Analysis", "Peptide Twins Checker"), index=0)

    # Default to landing_page when no option is selected
    if option == "HOME":
        landing_page()
    elif option == "Peptide Twins Analysis":
        peptide_twins_analysis()
        plot_spectra()
    elif option == "Peptide Twins Checker":
        peptide_twins_checker()

if __name__ == "__main__":
    main()
