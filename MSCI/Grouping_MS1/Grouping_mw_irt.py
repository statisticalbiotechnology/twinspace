import pandas as pd
import numpy as np
from itertools import combinations

def make_data_compatible(index_df):
    data_tuples = [(index, row['MW'], row['iRT']) for index, row in index_df.iterrows()]
    return data_tuples

compatible_data = make_data_compatible(mz_irt_df)
def within_ppm(pair, ppm_tolerance1, tolerance2):
    return (
        abs(pair[0][1] - pair[1][1]) <= (pair[0][1] * ppm_tolerance1) / 1e6 and
        abs(pair[0][2] - pair[1][2]) <= tolerance2
    )

def find_combinations_optimized(data, ppm_tolerance1, ppm_tolerance2):
    valid_combinations = []

    # Iterate through the data to find valid combinations
    for i in range(len(data)):
        print(i)
        for j in range(len(data)):
            # Check if the pair satisfies the tolerance conditions
            if within_ppm((data[i], data[j]), ppm_tolerance1, ppm_tolerance2):
                valid_combinations.append((data[i], data[j]))

    return valid_combinations

ppm_tolerance1 = 5
ppm_tolerance2 = 5

result_ppm = find_combinations_optimized(compatible_data, ppm_tolerance1, ppm_tolerance2)

# Calculate the elapsed time

unique_result_ppm = list(set(result_ppm))



