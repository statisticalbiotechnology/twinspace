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
        # Merge dataframes
        self.mz_index = {}
        merged = pd.merge(x, y, on='mz', how='outer', suffixes=('_x', '_y'))
        merged = merged.sort_values('mz')

        # Assign indices based on tolerance or ppm matching
        indices = []
        for index, row in merged.iterrows():
            matched = False
            for mz, idx in self.mz_index.items():
                ppm_tolerance = mz * self.ppm * 1e-6
                if abs(row['mz'] - mz) <= max(self.tolerance, ppm_tolerance):
                    indices.append(idx)
                    matched = True
                    break
            if not matched:
                new_index = len(self.mz_index)
                self.mz_index[row['mz']] = new_index
                indices.append(new_index)

        merged['index'] = indices

        # Define a custom filter to deduplicate while keeping up to two unique intensities per index
        def deduplicate(group):
            # If there are not more than two unique values, return the group
            if group['intensities_x'].nunique() <= 1 and group['intensities_y'].nunique() <= 1:
                return group.head(2)  # Keep at most two rows
            else:
                # Try to keep up to one of each unique intensity per column
                deduped = pd.concat([
                    group.drop_duplicates(subset=['intensities_x']).head(1),
                    group.drop_duplicates(subset=['intensities_y']).head(1)
                ]).drop_duplicates()
                return deduped

        # Apply the custom deduplication logic
        filtered = merged.groupby('index').apply(deduplicate).reset_index(drop=True)

        # Separate back into x and y dataframes
        x_filtered = filtered[['mz', 'intensities_x', 'index']].rename(columns={'intensities_x': 'intensities'})
        y_filtered = filtered[['mz', 'intensities_y', 'index']].rename(columns={'intensities_y': 'intensities'})

        # Handle NaN intensities and update m/z values accordingly
        x_filtered['intensities'] = x_filtered['intensities'].where(pd.notnull(x_filtered['intensities']), None)
        y_filtered['intensities'] = y_filtered['intensities'].where(pd.notnull(y_filtered['intensities']), None)
        x_filtered['mz'] = np.where(x_filtered['intensities'].isna(), None, x_filtered['mz'])
        y_filtered['mz'] = np.where(y_filtered['intensities'].isna(), None, y_filtered['mz'])

        return x_filtered, y_filtered


def process_combin(pair, spectra, tolerance, ppm):
    matcher = joinPeaks(tolerance=tolerance, ppm=ppm)
    # Unpack the pair into individual indices
    x_idx, y_idx = pair
    
    # Check if either spectra[x_idx] or spectra[y_idx] is None
    if spectra[x_idx] is None or spectra[y_idx] is None:
        # Handle the case where either spectra[x_idx] or spectra[y_idx] is None
        return np.nan
    
    # Extract the data for the two spectra
    x_data = {"mz": spectra[x_idx].peaks.mz, "intensities": spectra[x_idx].peaks.intensities}
    y_data = {"mz": spectra[y_idx].peaks.mz, "intensities": spectra[y_idx].peaks.intensities}
    # Convert the data to DataFrames
    x = pd.DataFrame(x_data)
    y = pd.DataFrame(y_data)
    
    # Match the peaks in the two spectra
    x, y = matcher.match(x, y)
    print(x)
    # Compute the angle between the two spectra
    angle = nspectraangle(x, y, m=0, n=0.5, na_rm=True)    
    return angle




from functools import partial
def process_combin_wrapper(combin, spectra, tolerance, ppm):
    return process_combin(combin, spectra, tolerance ,ppm)
def saveList(myList,filename):
    # the filename should mention the extension 'npy'
    np.save(filename,myList)
    print("Saved successfully!")
    

def loadList(filename):
    # the filename should mention the extension 'npy'
    tempNumpyArray=np.load(filename, allow_pickle=True)
    return tempNumpyArray.tolist()
