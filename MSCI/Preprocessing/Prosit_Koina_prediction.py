import re
import requests
import pandas as pd

# Modification masses
MODIFICATION = {
    "CAM": 57.021464,  # Carbamidomethylation
    "OX": 15.994915   # Oxidation
}

# Accurate amino acid masses including common modifications
AMINO_ACID = {
    "G": 57.021464,
    "R": 156.101111,
    "V": 99.068414,
    "P": 97.052764,
    "S": 87.032028,
    "U": 150.95363,
    "L": 113.084064,
    "M": 131.040485,
    "Q": 128.058578,
    "N": 114.042927,
    "Y": 163.063329,
    "E": 129.042593,
    "C": 103.009185 + MODIFICATION["CAM"],
    "F": 147.068414,
    "I": 113.084064,
    "A": 71.037114,
    "T": 101.047679,
    "W": 186.079313,
    "H": 137.058912,
    "D": 115.026943,
    "K": 128.094963,
}
AMINO_ACID["M(ox)"] = AMINO_ACID["M"] + MODIFICATION["OX"]

def calculate_molecular_weight_with_modifications(sequence):
    if not sequence:
        return 0  # Return 0 for empty sequences

    # Define a regular expression pattern to match modifications
    modification_pattern = r'\[UNIMOD:(\d+)\]'

    # Find all modification matches in the sequence
    modifications = re.findall(modification_pattern, sequence)

    # Remove modifications from the sequence
    sequence = re.sub(modification_pattern, '', sequence)

    # Calculate the molecular weight of the sequence without modifications
    base_molecular_weight = sum(AMINO_ACID[aa] for aa in sequence if aa in AMINO_ACID)

    # Add the mass of each modification to the molecular weight
    for mod_id in modifications:
        modification_mass = get_modification_mass(mod_id)
        if modification_mass is not None:
            base_molecular_weight += modification_mass

    # Add the mass of water (H2O) to account for the complete peptide
    molecular_weight_with_water = base_molecular_weight + 18.01056

    return molecular_weight_with_water

def get_modification_mass(mod_id):
    # Dictionary mapping UniMod IDs to modification masses
    modification_masses = {
        "4": MODIFICATION["CAM"],  # Carbamidomethylation on Cysteine
        "35": MODIFICATION["OX"]  # Oxidation on Methionine
        # Add more modifications if needed
    }
    return modification_masses.get(mod_id, 0)  # Return 0 if mod_id not found

def get_prosit_predictions(peptides, charges, collision_energies):
    # Construct the request body
    request_body = {
        "id": "batch_request",
        "inputs": [
            {
                "name": "peptide_sequences",
                "shape": [len(peptides), 1],
                "datatype": "BYTES",
                "data": list(peptides)
            },
            {
                "name": "collision_energies",
                "shape": [len(collision_energies), 1],
                "datatype": "FP32",
                "data": collision_energies
            },
            {
                "name": "precursor_charges",
                "shape": [len(charges), 1],
                "datatype": "INT32",
                "data": charges
            }
        ]
    }

    # Define the URL for the model endpoint
    url = "https://koina.wilhelmlab.org/v2/models/Prosit_2019_intensity/infer"

    # Send the POST request
    response = requests.post(url, json=request_body)

    # Check if the request was successful
    if response.status_code == 200:
        predictions = response.json()
        print(predictions['outputs'])

        
        # Extract and debug the intensity values
        try:
            intensities = predictions['outputs'][2]['data']
            mz_values = predictions['outputs'][1]['data'] 
            print(mz_values)
            print(f"Number of peptides: {len(peptides)}")
        except KeyError:
            print("The response does not contain the expected 'outputs' key.")
            print(predictions)
            return None
        
        # Assume each peptide has equal number of intensity values
        num_intensities_per_peptide = len(intensities) // len(peptides)
        print(f"Number of intensity values per peptide: {num_intensities_per_peptide}")

        # Create a DataFrame with the peptide sequences, charges, collision energies, and intensity values
        rows = []
        for i, peptide in enumerate(peptides):
            peptide_intensities = intensities[i * num_intensities_per_peptide : (i + 1) * num_intensities_per_peptide]
            peptide_mz_values = mz_values[i * num_intensities_per_peptide : (i + 1) * num_intensities_per_peptide]
            
            row = {
                'peptide_sequence': peptide,
                'charge': charges[i],
                'collision_energy': collision_energies[i],
                'mz_values': peptide_mz_values,
                'intensities': peptide_intensities
            }
            rows.append(row)
        
        df = pd.DataFrame(rows)
        return df
    else:
        print(f"Request failed with status code {response.status_code}")
        print(response.text)
        return None

def get_irt_predictions(peptides):
    # Construct the request body
    request_body = {
        "id": "test_id",
        "inputs": [
            {
                "name": "peptide_sequences",
                "shape": [len(peptides), 1],
                "datatype": "BYTES",
                "data": peptides
            }
        ]
    }

    # Define the URL for the model endpoint
    url = "https://koina.wilhelmlab.org/v2/models/Prosit_2019_irt/infer"

    # Send the POST request
    response = requests.post(url, json=request_body)

    # Check if the request was successful
    if response.status_code == 200:
        predictions = response.json()
        irt_values = predictions['outputs'][0]['data']
        return irt_values
    else:
        print(f"Request failed with status code {response.status_code}")
        print(response.text)
        return None

def format_msp(peptide, charge, collision_energy, mz_values, intensities, irt):
    mw = calculate_molecular_weight_with_modifications(peptide)
    mz = (mw + (charge)) / charge  # Correct m/z calculation
    header = f"Name: {peptide}/{charge}\n"
    pepmass = f"MW: {mz:.6f}\n"
    collision_energy_line = f"Collision_energy: {collision_energy}\n"
    irt_line = f"iRT: {irt:.6f}\n"

    # Combine mz_values and intensities, filter out -1.00 entries, and then sort by mz_values
    sorted_peaks = sorted(
        (mz, intensity) for mz, intensity in zip(mz_values, intensities) if mz != -1.00 and intensity != -1.000000
    )

    # Update the number of valid peaks
    num_peaks = len(sorted_peaks)
    peaks_header = f"Num peaks: {num_peaks}\n"

    peaks = ""
    # Normalize intensities
    if sorted_peaks:
        max_intensity = max(intensity for mz, intensity in sorted_peaks) if max(intensity for mz, intensity in sorted_peaks) != 0 else 1
        normalized_intensities = [intensity / max_intensity * 100 for mz, intensity in sorted_peaks]
    else:
        max_intensity = 1
        normalized_intensities = []

    for (mz, _), intensity in zip(sorted_peaks, normalized_intensities):
        peaks += f"{mz:.2f}\t{intensity:.6f}\n"

    return f"{header}{pepmass}{collision_energy_line}{irt_line}{peaks_header}{peaks}\n"

def save_to_msp(df, filename, irt_values):
    with open(filename, 'a') as file:
        for index, row in df.iterrows():
            irt = irt_values[index]
            msp_entry = format_msp(row['peptide_sequence'], row['charge'], row['collision_energy'], row['mz_values'], row['intensities'], irt)
            file.write(msp_entry)

# Initialize the MSP file (clear existing content)
with open('peptides.msp', 'w') as file:
    file.write("")

peptides = df_peptides['peptide'].tolist()
batch_size = 1000

# Process peptides with specified charges
def process_peptides_with_charge(peptides, charge, filename):
    for start in range(0, len(peptides), batch_size):
        batch_peptides = peptides[start:start + batch_size]
        charges = [charge] * len(batch_peptides)
        collision_energies = [28] * len(batch_peptides)

        # Get iRT predictions
        irt_values = get_irt_predictions(batch_peptides)

        # Get Prosit predictions
        df = get_prosit_predictions(batch_peptides, charges, collision_energies)

        if df is not None and irt_values is not None:
            # Append the DataFrame to the MSP file with iRT values
            save_to_msp(df, filename, irt_values)
            print(f"Processed batch {start // batch_size + 1} with charge {charge}")

# Process with charge 2
process_peptides_with_charge(peptides, 2, 'peptides.msp')

# Process with charge 3
process_peptides_with_charge(peptides, 3, 'peptides.msp')

print("MSP file saved successfully.")
