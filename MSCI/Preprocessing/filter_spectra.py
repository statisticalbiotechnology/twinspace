import matchms.filtering as filtering
import pickle

def keep_top_n_peaks(spectrum, n: int):
    """Keep only the top n most intense peaks."""
    return filtering.reduce_to_number_of_peaks(spectrum, n_required=n, n_max=n)

def filter_spectra_by_top_peaks(input_file_path: str, output_file_path: str, n_peaks: int):
    """Filter spectra to keep only the top n peaks."""
    with open(input_file_path, "rb") as open_file:
        spectra = pickle.load(open_file)
    
    processed_spectra = [
        keep_top_n_peaks(spectrum, n_peaks) for spectrum in spectra
    ]
    
    with open(output_file_path, "wb") as output_file:
        pickle.dump(processed_spectra, output_file)

    return processed_spectra
