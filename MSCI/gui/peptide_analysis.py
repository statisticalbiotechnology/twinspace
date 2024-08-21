# peptide_analysis.py
import streamlit as st
import tempfile
import pandas as pd
from MSCI.Preprocessing.Koina import PeptideProcessor
from .utils import load_image
from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs
from MSCI.Preprocessing.read_msp_file import read_msp_file
from matchms.importing import load_from_msp
from matchms.exporting import save_as_msp
from matchms.filtering import reduce_to_number_of_peaks, select_by_relative_intensity
import os
import matplotlib.pyplot as plt
import sys
from matchms import Spectrum
from MSCI.Preprocessing.Koina import PeptideProcessor
from MSCI.Grouping_MS1.Grouping_mw_irt import process_peptide_combinations
from MSCI.Preprocessing.read_msp_file import read_msp_file
from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs
# Constants
INTENSITY_MODELS = [
    "Prosit_2020_intensity_HCD", "ms2pip_HCD2021", "ms2pip_timsTOF2023", "ms2pip_iTRAQphospho",
    "ms2pip_Immuno_HCD", "ms2pip_TTOF5600", "ms2pip_timsTOF2024", "ms2pip_CID_TMT",
    "Prosit_2019_intensity", "Prosit_2023_intensity_timsTOF",
    "Prosit_2020_intensity_CID", "Prosit_2023_intensity_XL_CMS2", "Prosit_2020_intensity_TMT",
    "Prosit_2023_intensity_XL_CMS3"
]
IRT_MODELS = ["Prosit_2019_irt", "Deeplc_hela_hf", "AlphaPeptDeep_rt_generic", "Prosit_2020_irt_TMT"]


import os
def keep_top_n_peaks(spectrum: Spectrum, n: int) -> Spectrum:
    """Keep only the top n most intense peaks in the spectrum."""
    return reduce_to_number_of_peaks(spectrum, n_required=n, n_max=n)


def keep_above_intensity(spectrum: Spectrum, threshold: float) -> Spectrum:
    """Keep only peaks above a certain relative intensity threshold."""
    spectrum = select_by_relative_intensity(spectrum, intensity_from=threshold, intensity_to=1.0)
    return spectrum

def filter_spectra_by_intensity(input_file_path: str, output_file_path: str, intensity_threshold: float):
    """Filter the spectra by keeping only peaks above a certain intensity threshold."""
    try:
        # Load spectra from the MSP file
        spectra = list(load_from_msp(input_file_path))

        # Loop through all spectra and apply keep_above_intensity function
        processed_spectra = []
        for spectrum in spectra:
            processed_spectrum = keep_above_intensity(spectrum, intensity_threshold)
            processed_spectra.append(processed_spectrum)

        # Ensure the output file is either non-existent or empty before saving
        if os.path.exists(output_file_path):
            os.remove(output_file_path)

        # Save the processed spectra as an MSP file
        save_as_msp(processed_spectra, output_file_path)

    except Exception as e:
        st.error(f"An error occurred while filtering spectra: {e}")
        return None

def filter_spectra_by_top_peaks(input_file_path: str, output_file_path: str, n_peaks: int):
    """Filter the spectra by keeping only the top n peaks and format them correctly."""
    try:
        # Load spectra from the MSP file
        spectra = list(load_from_msp(input_file_path))

        # Loop through all spectra and apply keep_top_n_peaks function
        processed_spectra = []
        for spectrum in spectra:
            processed_spectrum = keep_top_n_peaks(spectrum, n_peaks)
            processed_spectra.append(processed_spectrum)

        # Ensure the output file is either non-existent or empty before saving
        if os.path.exists(output_file_path):
            os.remove(output_file_path)

        # Save the processed spectra as an MSP file
        save_as_msp(processed_spectra, output_file_path)

    except Exception as e:
        st.error(f"An error occurred while filtering spectra: {e}")
        return None

def perform_analysis(mz_tolerance: float, irt_tolerance: float, use_ppm: bool):
    """Perform the peptide twins analysis with detailed debugging."""
    if not st.session_state.spectra_cache or st.session_state.mz_irt_df_cache.empty:
        st.error("Spectra data is missing or invalid. Please ensure the MSP file is correctly loaded.")
        return

    with st.spinner("Processing peptide combinations..."):
        try:
            Groups_df = process_peptide_combinations(
                st.session_state.mz_irt_df_cache, mz_tolerance, irt_tolerance, use_ppm=use_ppm
            )

            st.write(f"Grouped data shape: {Groups_df.shape}")
#            st.write(Groups_df.head())  # Log the first few rows for debugging

            if not {'index1', 'index2'}.issubset(Groups_df.columns):
                st.error("Unexpected DataFrame structure. 'index1' and 'index2' columns not found.")
                return

            index_array = Groups_df[['index1', 'index2']].values.astype(int)
#            st.write(f"Number of combinations: {len(index_array)}")

            similarity_progress = st.progress(0)
            total_combinations = len(index_array)
            results = []

            with st.spinner("Calculating spectra similarities..."):
                for i, (idx1, idx2) in enumerate(index_array):
                    result = process_spectra_pairs(
                        [(idx1, idx2)], st.session_state.spectra_cache, st.session_state.mz_irt_df_cache, tolerance=mz_tolerance, ppm=use_ppm
                    )

                    results.append(result)
                    similarity_progress.progress((i + 1) / total_combinations)

                st.session_state.analysis_results = pd.concat(results, ignore_index=True)
#                st.subheader("Spectra Similarity Results:")
#                st.dataframe(st.session_state.analysis_results)

        except Exception as e:
            st.error(f"An error occurred during analysis: {str(e)}")


import tempfile

def peptide_twins_analysis():
    """Render the Peptide Twins Analysis page."""
    st.session_state.setdefault('spectra_cache', None)
    st.session_state.setdefault('mz_irt_df_cache', None)
    st.session_state.setdefault('temp_file_path', None)

    st.header("Peptide Twins Analysis")

    uploaded_file = st.file_uploader("Upload your peptide file", type=["txt"])
    model_intensity = st.selectbox("Select Intensity Model", INTENSITY_MODELS)
    model_irt = st.selectbox("Select iRT Model", IRT_MODELS)
    collision_energy = st.number_input("Set Collision Energy", min_value=1, step=1, value=30)
    charge = st.number_input("Set Charge", min_value=1, step=1, value=2)

    filter_option = st.selectbox("Filter Option", ["None", "Top N Peaks", "Intensity Threshold"])
    n_peaks = st.number_input("Number of Peaks to Keep", min_value=1, step=1, value=6, key="n_peaks")
    intensity_threshold = st.number_input("Intensity Threshold", min_value=0.0, max_value=1.0, step=0.01, value=0.1, key="intensity_threshold")

    if uploaded_file:
        st.session_state.peptide_data = uploaded_file.read().decode("utf-8")

        if not st.session_state.temp_file_path:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode='w') as temp_file:
                temp_file.write(st.session_state.peptide_data)
                st.session_state.temp_file_path = temp_file.name

        st.subheader("Prediction Results")

        processor = PeptideProcessor(
            input_file=st.session_state.temp_file_path,
            collision_energy=collision_energy,
            charge=charge,
            model_intensity=model_intensity,
            model_irt=model_irt
        )

        try:
            with st.spinner("Running prediction..."):
                with tempfile.NamedTemporaryFile(delete=False, suffix=".msp") as temp_output:
                    output_msp_path = temp_output.name
                    processor.process(output_msp_path)
                st.success("Prediction Completed Successfully")

            st.subheader("Spectra Analysis")

            mz_tolerance = st.number_input("Set Mass Tolerance", min_value=0.0, step=0.1, value=0.1)
            irt_tolerance = st.number_input("Set iRT Tolerance", min_value=0.0, step=0.1, value=0.1)
            use_ppm = st.checkbox("Use PPM for mass tolerance", value=False)

            similarity_method = st.selectbox(
                "Select Similarity Method",
                ("Spectral Angle", "Greedy Cosine Similarity")
            )

            st.write(f"Selected similarity method: {similarity_method}")

            if st.button("Start Analysis"):
                st.session_state.similarity_method = similarity_method

                spectra_file = output_msp_path
                
                if filter_option == "Top N Peaks":
                    try:
                        with tempfile.NamedTemporaryFile(delete=False, suffix=".msp") as temp_filtered:
                            filtered_file = temp_filtered.name
                            filter_spectra_by_top_peaks(
                                spectra_file,
                                filtered_file,
                                st.session_state.n_peaks
                            )
                            spectra_file = filtered_file
                    except Exception as e:
                        st.error(f"An error occurred while filtering spectra: {e}")
                        return

                elif filter_option == "Intensity Threshold":
                    try:
                        with tempfile.NamedTemporaryFile(delete=False, suffix=".msp") as temp_filtered:
                            filtered_file = temp_filtered.name
                            filter_spectra_by_intensity(
                                spectra_file,
                                filtered_file,
                                st.session_state.intensity_threshold
                            )
                            spectra_file = filtered_file
                    except Exception as e:
                        st.error(f"An error occurred while filtering spectra: {e}")
                        return

                try:
                    st.session_state.spectra_cache = list(load_from_msp(spectra_file))
                    st.session_state.mz_irt_df_cache = read_msp_file(spectra_file)
                    st.write(f"Loaded {len(st.session_state.spectra_cache)} spectra from the MSP file.")
                except Exception as e:
                    st.error(f"An error occurred while loading spectra: {e}")
                    return

                perform_analysis(mz_tolerance, irt_tolerance, use_ppm)

        except Exception as e:
            st.error(f"An error occurred during prediction: {e}")
    else:
        st.info("Please upload a peptide text file to proceed.")
