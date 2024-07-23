import numpy as np
import pandas as pd 

"""
ndotproduct(x, y, m=0, n=0.5, na_rm=True): Computes the normalized dot product between two spectra x and y. The function uses a weighting mechanism defined by _weightxy to weight the m/z and intensity values of the spectra.
nspectraangle(x, y, m=0, n=0.5, na_rm=True): Calculates the normalized spectral angle between two spectra. The angle is derived from the normalized dot product.
_weightxy(x, y, m=0, n=0.5): A helper function that returns weighted values based on the provided m/z (x) and intensity (y) values.

The joinPeaks class is designed for matching peaks between two spectra based on a specified tolerance and parts per million (ppm).
Upon initialization, the class sets the tolerance and ppm values and initializes an empty dictionary mz_index to store m/z values and their corresponding indices.
The match method merges two spectra, sorts them by m/z values, and assigns indices to peaks. Peaks within the specified tolerance are matched and assigned the same index. The method then filters out duplicated indices with NaN intensities and updates m/z values where intensity is missing. The final matched spectra are returned as two separate data frames.
Usage:
# For functions
dot_product_value = ndotproduct(x, y)
angle_value = nspectraangle(x, y)
joiner = joinPeaks(tolerance=value1, ppm=value2)
matched_x, matched_y = joiner.match(x, y)
"""
def ndotproduct(x, y, m=0, n=0.5, na_rm=True):
    wx = _weightxy(x.iloc[:,0], x.iloc[:,1], m, n)
    wy = _weightxy(y.iloc[:,0], y.iloc[:,1], m, n)
    wx2 = wx**2
    wy2 = wy**2
    return (np.sum(wx * wy)**2) / (np.sum(wx2, axis=0) * np.sum(wy2, axis=0))



def nspectraangle(x, y, m=0, n=0.5, na_rm=True):
    return 1 - 2 * np.arccos(ndotproduct(x, y, m, n, na_rm)) / np.pi

def _weightxy(x, y, m=0, n=0.5):
    return x**m * y**n

import pandas as pd
import numpy as np

class joinPeaks:
    def __init__(self, tolerance=0, ppm=0):
        self.tolerance = tolerance
        self.ppm = ppm

    def match(self, x, y):
        # Initialize the result lists
        x_matched = []
        y_matched = []

        # Create a list to keep track of matched y indices
        matched_y_indices = set()

        # Calculate ppm tolerance for each mz value in x and y
        ppm_tolerance_x = x['mz'] * self.ppm * 1e-6

        # Perform the matching
        for i, x_row in x.iterrows():
            mz_x = x_row['mz']
            intensity_x = x_row['intensities']
            match_found = False

            for j, y_row in y.iterrows():
                if j in matched_y_indices:
                    continue

                mz_y = y_row['mz']
                intensity_y = y_row['intensities']
                ppm_tolerance_y = mz_y * self.ppm * 1e-6

                if abs(mz_y - mz_x) <= self.tolerance or abs(mz_y - mz_x) <= ppm_tolerance_x[i] or abs(mz_y - mz_x) <= ppm_tolerance_y:
                    x_matched.append([mz_x, intensity_x])
                    y_matched.append([mz_y, intensity_y])
                    matched_y_indices.add(j)
                    match_found = True
                    break

            if not match_found:
                x_matched.append([mz_x, intensity_x])
                y_matched.append([np.nan, np.nan])

        # Add unmatched rows from y
        for j, y_row in y.iterrows():
            if j not in matched_y_indices:
                mz_y = y_row['mz']
                intensity_y = y_row['intensities']
                x_matched.append([np.nan, np.nan])
                y_matched.append([mz_y, intensity_y])

        # Convert to DataFrames
        x_matched_df = pd.DataFrame(x_matched, columns=['mz', 'intensities'])
        y_matched_df = pd.DataFrame(y_matched, columns=['mz', 'intensities'])

        # Set index to have the same structure
        x_matched_df['index'] = range(len(x_matched_df))
        y_matched_df['index'] = range(len(y_matched_df))

        x_matched_df.set_index('index', inplace=True)
        y_matched_df.set_index('index', inplace=True)

        return x_matched_df, y_matched_df


def process_spectra_pairs(chunk, spectra, mz_irt_df, tolerance=0, ppm=10, m=0, n=0.5):
    """
    Processes pairs of spectra, matches their peaks, and calculates the angle between them.

    Parameters:
    chunk (list): List of index pairs (tuples) to be processed.
    spectra (list): List of spectra objects.
    mz_irt_df (pd.DataFrame): DataFrame containing MW, iRT, and peptide names.
    tolerance (int): Tolerance for peak matching.
    ppm (int): Parts per million for peak matching.
    m (float): Parameter for angle calculation.
    n (float): Parameter for angle calculation.

    Returns:
    pd.DataFrame: DataFrame containing indices of the spectra, calculated angle, and similarity score.
    """
    results = []

    for index_pair in chunk:
        i, j = index_pair

        x = spectra[i]
        y = spectra[j]

        x_df = pd.DataFrame({'mz': x.peaks.mz, 'intensities': x.peaks.intensities})
        y_df = pd.DataFrame({'mz': y.peaks.mz, 'intensities': y.peaks.intensities})

        matcher = joinPeaks(tolerance=tolerance, ppm=ppm)
        x_matched, y_matched = matcher.match(x_df, y_df)

        angle = nspectraangle(x_matched, y_matched, m=m, n=n)
        
        # Extract the relevant information for the given index pair
        results.append({
            'index1': i,
            'index2': j,
            'column1_peptide': mz_irt_df.loc[i, 'Name'],
            'column2_peptide': mz_irt_df.loc[j, 'Name'],
            'MW1': mz_irt_df.loc[i, 'MW'],
            'MW2': mz_irt_df.loc[j, 'MW'],
            'iRT1': mz_irt_df.loc[i, 'iRT'],
            'iRT2': mz_irt_df.loc[j, 'iRT'],
            'similarity_score': angle,  
        })

    return pd.DataFrame(results)
