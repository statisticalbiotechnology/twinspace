import pandas as pd
def read_msp_file(filename):
    """
    This function is designed to parse and extract specific spectral data from an MSP formatted file. 
    It processes the file line-by-line, capturing the name of the spectrum, molecular weight (MW),and  the iRT value. 
    The function then converts the list of dictionaries into pandas DataFrame 
        with columns 'Name', 'MW', and 'iRT', which is returned as the final output.

    Usage:
    df = read_msp_file(filename)
"""
    spectra = []
    current_spectrum = {}
    
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith('Name: '):
                if current_spectrum:
                    spectra.append(current_spectrum)
                current_spectrum = {'Name': line.split('Name: ')[1]}
            elif line.startswith('MW: '):
                current_spectrum['MW'] = float(line.split('MW: ')[1])
            elif line.startswith('iRT: '):
                current_spectrum['iRT'] = float(line.split('iRT: ')[1])

    if current_spectrum:
        spectra.append(current_spectrum)

    index_df = pd.DataFrame(spectra, columns=['Name', 'MW', 'iRT'])
    return index_df

