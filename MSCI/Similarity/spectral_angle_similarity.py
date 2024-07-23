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


def process_spectra_pairs(chunk, spectra, mz_irt_df, unique_result_ppm, tolerance=0, ppm=10, m=0, n=0.5):
    """


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
        
        # Find the corresponding unique result ppm pair
        for (index1, mw1, irt1), (index2, mw2, irt2) in unique_result_ppm:
            if (i == index1 and j == index2) or (i == index2 and j == index1):
                similarity_score = angle # You can replace this with your own similarity score calculation if needed
                results.append({
                    'index1': index1,
                    'index2': index2,
                    'column1_peptide': mz_irt_df.loc[index1, 'Name'],
                    'column2_peptide': mz_irt_df.loc[index2, 'Name'],
                    'MW1': mw1,
                    'MW2': mw2,
                    'iRT1': irt1,
                    'iRT2': irt2,
                    'angle': angle,
                    'similarity_score': similarity_score
                })

    return pd.DataFrame(results)
