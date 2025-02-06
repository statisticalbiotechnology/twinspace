import pandas as pd
from pyopenms import MzMLFile, MSExperiment
from pyteomics import mgf

def read_msp_file(filename):
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

def process_spectrum(spectrum):
    precursors = spectrum.getPrecursors()
    if precursors:
        mw = precursors[0].getMZ()
        rt = spectrum.getRT()
        num_peaks = spectrum.size()
        peaks = [(peak.getMZ(), peak.getIntensity()) for peak in spectrum]
        return {'MW': mw, 'RT': rt, 'Num Peaks': num_peaks, 'Peaks': peaks}
    return None

def read_mgf_file(filename):
    spectra = []
    with mgf.read(filename) as mgf_file:
        for spectrum in mgf_file:
            mz_values = spectrum['m/z array']
            intensities = spectrum['intensity array']
            mw = spectrum['params'].get('pepmass', [None])[0]
            rt = spectrum['params'].get('rtinseconds', None)
            if rt is None:
                print(f"Warning: No retention time found for spectrum with MW {mw}.")
            spectra.append({
                'mz_values': mz_values,
                'intensities': intensities,
                'MW': mw,
                'RT': rt
            })
            print(f"Spectrum with {len(mz_values)} m/z values and {len(intensities)} intensities. MW: {mw}, RT: {rt}")
    return spectra

def read_mzml_file(filename):
    spectra_data = []
    exp = MSExperiment()
    MzMLFile().load(filename, exp)

    for spectrum in exp:
        processed_spectrum = process_spectrum(spectrum)
        if processed_spectrum:
            spectra_data.append(processed_spectrum)
    
    return spectra_data

def read_ms_file(filename):
    file_extension = filename.split('.')[-1].lower()
    if file_extension == 'msp':
        return read_msp_file(filename)
    elif file_extension == 'mzml':
        return read_mzml_file(filename)
    elif file_extension == 'mgf':
        return read_mgf_file(filename)
    else:
        raise ValueError(f"Unsupported file format: {file_extension}")

# Example usage
file_path = 'Z:/zelhamraoui/MSCA_Package/real_data/example.mgf'
result = read_ms_file(file_path)
if isinstance(result, pd.DataFrame):
    display(result)
else:
    print(result)
