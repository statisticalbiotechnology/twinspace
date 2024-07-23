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

def find_combinations_kdtree(data, ppm_tolerance1, ppm_tolerance2):
    valid_combinations = []
    
    # Create numpy array for k-d tree
    data_array = np.array([(mw, irt) for index, mw, irt in data])
    
    # Build k-d tree
    tree = cKDTree(data_array)
    
    # Pre-compute the square of tolerances to avoid repetitive computation
    ppm_tolerance1_sq = (ppm_tolerance1 / 1e6)**2
    ppm_tolerance2_sq = ppm_tolerance2**2
    
    # Iterate through the data to find valid combinations
    for i in range(len(data)):
        current_point = data_array[i]
        # Calculate a combined radius tolerance for the query
        radius = np.sqrt((ppm_tolerance1_sq * current_point[0]**2) + ppm_tolerance2_sq)
        indices = tree.query_ball_point(current_point, radius)
        
        for j in indices:
            if i < j:  # Avoid self-pairing and duplicate pairs
                pair = (data[i], data[j])
                if within_ppm(pair, ppm_tolerance1, ppm_tolerance2):
                    valid_combinations.append(pair)

    return valid_combinations

def process_peptide_combinations(mz_irt_df, ppm_tolerance1, ppm_tolerance2):
    compatible_data = make_data_compatible(mz_irt_df)
    result_ppm = find_combinations_kdtree(compatible_data, ppm_tolerance1, ppm_tolerance2)
    unique_result_ppm = list({tuple(sorted(pair)) for pair in result_ppm})
    
    # Create a DataFrame for the results
    results = []
    for (index1, mw1, irt1), (index2, mw2, irt2) in unique_result_ppm:
        results.append({
            'index1': index1,
            'index2': index2,
            'column1_peptide': mz_irt_df.loc[index1, 'Name'],
            'column2_peptide': mz_irt_df.loc[index2, 'Name'],
            'MW1': mw1,
            'MW2': mw2,
            'iRT1': irt1,
            'iRT2': irt2
        })

    results_df = pd.DataFrame(results)
    return results_df
