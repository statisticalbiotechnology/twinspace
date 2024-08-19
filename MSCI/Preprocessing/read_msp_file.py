import pandas as pd

def read_msp_file(filename):
    """
    This function is designed to parse and extract specific spectral data from either of two MSP formatted files. 
    It processes the file line-by-line, capturing relevant fields such as the compound name, nominal mass, IRT, molecular weight (MW),
    and collision energy. The function then standardizes the data into columns 'Name', 'MW', and 'iRT', which is returned as the final output.

    Usage:
    df = read_msp_file(filename)
    """
    spectra = []
    current_spectrum = {}
    
    with open(filename) as f:
        for line in f:
            line = line.strip()
            # Detect format type based on the line content
            if line.startswith('Name: '):
                if current_spectrum:
                    spectra.append(current_spectrum)
                current_spectrum = {'Name': line.split('Name: ')[1]}
            elif line.startswith('MW: '):
                current_spectrum['MW'] = float(line.split('MW: ')[1])
            elif line.startswith('iRT: '):
                current_spectrum['iRT'] = float(line.split('iRT: ')[1])
            elif line.startswith('COLLISION_ENERGY: '):
                if current_spectrum:
                    spectra.append(current_spectrum)
                # Assuming COMPOUND_NAME comes after COLLISION_ENERGY in the file
                current_spectrum = {}
            elif line.startswith('IRT: '):
                current_spectrum['iRT'] = float(line.split('IRT: ')[1])
            elif line.startswith('COMPOUND_NAME: '):
                current_spectrum['Name'] = line.split('COMPOUND_NAME: ')[1] if 'COMPOUND_NAME: ' in line else None
            elif line.startswith('NOMINAL_MASS: '):
                current_spectrum['MW'] = float(line.split('NOMINAL_MASS: ')[1]) if 'NOMINAL_MASS: ' in line else None

    if current_spectrum:
        spectra.append(current_spectrum)

    # Standardize the DataFrame columns to 'Name', 'MW', and 'iRT'
    index_df = pd.DataFrame(spectra, columns=['Name', 'MW', 'iRT'])
    
    return index_df
