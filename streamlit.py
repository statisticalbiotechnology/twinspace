import streamlit as st
from pathlib import Path
import MSCI
from MSCI.gui.utils import load_image, add_custom_css
from MSCI.gui.landing_page import landing_page
from MSCI.gui.peptide_analysis import peptide_twins_analysis, plot_spectra
from MSCI.gui.peptide_checker import peptide_twins_checker

def main():
    st.set_page_config(layout="wide")
    add_custom_css()

    logo_path = Path("images/MSCI_logor.png")
    logo_image = load_image(str(logo_path))

    with st.sidebar:
        if logo_image:
            st.markdown(f"""
            <p align="center">
                <img src="data:image/png;base64,{logo_image}" alt="logo" width="300" height="300">
            </p>
            """, unsafe_allow_html=True)
        st.header("MSCI")
        option = st.radio("Choose an option", ("MSCI", "Peptide Twins Analysis", "Peptide Twins Checker"))

    if option == "MSCI":
        landing_page()
    elif option == "Peptide Twins Analysis":
        peptide_twins_analysis()
        plot_spectra()
    elif option == "Peptide Twins Checker":
        peptide_twins_checker()

if __name__ == "__main__":
    main()
