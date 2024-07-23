import pandas as pd
from scipy.spatial import cKDTree
import numpy as np

def make_data_compatible(index_df):
    data_tuples = [(index, row['MW'], row['iRT']) for index, row in index_df.iterrows()]
    return data_tuples

def within_ppm(pair, ppm_tolerance1, ppm_tolerance2):
    return (
        abs(pair[0][1] - pair[1][1]) <= (pair[0][1] * ppm_tolerance1) / 1e6 and
        abs(pair[0][2] - pair[1][2]) <= ppm_tolerance2
    )

def within_tolerance(pair, tolerance1, tolerance2):
    return (
        abs(pair[0][1] - pair[1][1]) <= tolerance1 and
        abs(pair[0][2] - pair[1][2]) <= tolerance2
    )

def find_combinations_kdtree(data, tolerance1, tolerance2, use_ppm=True):
    valid_combinations = []
    
    # Create numpy array for k-d tree
    data_array = np.array([(mw, irt) for index, mw, irt in data])
    
    # Build k-d tree
    tree = cKDTree(data_array)
    
    if use_ppm:
        # Pre-compute the square of tolerances to avoid repetitive computation
        ppm_tolerance1_sq = (tolerance1 / 1e6)**2
        ppm_tolerance2_sq = tolerance2**2
    
    # Iterate through the data to find valid combinations
    for i in range(len(data)):
        current_point = data_array[i]
        
        if use_ppm:
            # Calculate a combined radius tolerance for the query
            radius = np.sqrt((ppm_tolerance1_sq * current_point[0]**2) + ppm_tolerance2_sq)
        else:
            radius = np.sqrt(tolerance1**2 + tolerance2**2)
        
        indices = tree.query_ball_point(current_point, radius)
        
        for j in indices:
            if i < j:  # Avoid self-pairing and duplicate pairs
                pair = (data[i], data[j])
                if use_ppm:
                    if within_ppm(pair, tolerance1, tolerance2):
                        valid_combinations.append(pair)
                else:
                    if within_tolerance(pair, tolerance1, tolerance2):
                        valid_combinations.append(pair)

    return valid_combinations

def process_peptide_combinations(mz_irt_df, tolerance1, tolerance2, use_ppm=True):
    compatible_data = make_data_compatible(mz_irt_df)
    result_tolerance = find_combinations_kdtree(compatible_data, tolerance1, tolerance2, use_ppm)
    unique_result_tolerance = list({tuple(sorted(pair)) for pair in result_tolerance})
    
    # Create a DataFrame for the results
    results = []
    for (index1, mw1, irt1), (index2, mw2, irt2) in unique_result_tolerance:
        results.append({
            'index1': index1,
            'index2': index2,
            'peptide 1': mz_irt_df.loc[index1, 'Name'],
            'peptide 2': mz_irt_df.loc[index2, 'Name'],
            'm/z  1': mz_irt_df.loc[i, 'MW'],
            'm/z 2': mz_irt_df.loc[j, 'MW'],
            'iRT 1': mz_irt_df.loc[i, 'iRT'],
            'iRT 2': mz_irt_df.loc[j, 'iRT'],
        })

    results_df = pd.DataFrame(results)
    return results_df
