import requests
import pandas as pd
import re

# Dummy data for AA_MASSES, MOD_MASSES, MASSES to avoid errors
AA_MASSES = {'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694, 'C': 103.00919,
             'E': 129.04259, 'Q': 128.05858, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
             'L': 113.08406, 'K': 128.09496, 'M': 131.04049, 'F': 147.06841, 'P': 97.05276,
             'S': 87.03203, 'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841}
MOD_MASSES = {'[UNIMOD:4]': 15.99491}  # Example modification mass
MASSES = {'N_TERMINUS': 1.007825, 'C_TERMINUS': 17.00274}

class PeptideProcessor:
    def __init__(self, peptide_file, ce, charge):
        self.peptide_file = peptide_file
        self.ce = ce
        self.charge = charge
        self.peptides = self.read_peptides_from_file(peptide_file)
        self.df_peptides = pd.DataFrame(self.peptides, columns=["Peptide"]).drop_duplicates()

    def read_peptides_from_file(self, file_path):
        with open(file_path, "r") as file:
            lines = file.readlines()
        lines = [line.strip() for line in lines if line.strip()]
        print(f"Peptides read from file: {lines}")
        return lines

    def calculate_peptide_mass(self, peptide_sequence):
        total_mass = MASSES["N_TERMINUS"]
        pattern = re.compile(r'([A-Z])(\[UNIMOD:\d+\])?')
        for match in pattern.finditer(peptide_sequence):
            aa = match.group(1)
            mod = match.group(2)
            if aa in AA_MASSES:
                total_mass += AA_MASSES[aa]
                if mod:
                    mod_mass = MOD_MASSES.get(mod, 0.0)
                    total_mass += mod_mass
            else:
                print(f"Unknown amino acid: {aa}")
                return None
        total_mass += MASSES["C_TERMINUS"]
        return total_mass

    def get_prosit_predictions(self, peptides):
        print("Fetching Prosit predictions...")
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
                    "shape": [len(peptides), 1],
                    "datatype": "FP32",
                    "data": [self.ce] * len(peptides)
                },
                {
                    "name": "precursor_charges",
                    "shape": [len(peptides), 1],
                    "datatype": "INT32",
                    "data": [self.charge] * len(peptides)
                }
            ]
        }
        url = "https://koina.wilhelmlab.org/v2/models/Prosit_2020_intensity_HCD/infer"
        response = requests.post(url, json=request_body)
        if response.status_code == 200:
            print("Received response from Prosit.")
            predictions = response.json()
            try:
                intensities = predictions['outputs'][2]['data']
                mz_values = predictions['outputs'][1]['data']
            except KeyError:
                print("Error: KeyError in the predictions data.")
                return None
            num_intensities_per_peptide = len(intensities) // len(peptides)
            rows = []
            for i, peptide in enumerate(peptides):
                peptide_intensities = intensities[i * num_intensities_per_peptide : (i + 1) * num_intensities_per_peptide]
                peptide_mz_values = mz_values[i * num_intensities_per_peptide : (i + 1) * num_intensities_per_peptide]
                rows.append({
                    'peptide_sequence': peptide,
                    'charge': self.charge,
                    'collision_energy': self.ce,
                    'mz_values': peptide_mz_values,
                    'intensities': peptide_intensities
                })
            print("Prosit predictions processed.")
            return pd.DataFrame(rows)
        else:
            print(f"Error: Failed to get Prosit predictions. Status code: {response.status_code}")
            return None

    def save_to_msp(self, df, file_path):
        print(f"Saving predictions to {file_path}...")
        with open(file_path, 'w') as file:
            for _, row in df.iterrows():
                mw = self.calculate_peptide_mass(row['peptide_sequence'])
                mz = (mw + (self.charge * 1.007276)) / self.charge
                header = f"Name: {row['peptide_sequence']}/{self.charge}\nMW: {mz:.6f}\nCollision_energy: {self.ce}\n"
                peaks = "".join(f"{mz:.2f}\t{intensity:.6f}\n" for mz, intensity in zip(row['mz_values'], row['intensities']) if mz != -1.00 and intensity != -1.000000)
                file.write(f"{header}Num peaks: {len(peaks)}\n{peaks}\n")
        print("File saved successfully.")

    def process(self, output_file='output.msp'):
        print("Starting peptide processing...")
        df = self.get_prosit_predictions(self.peptides)
        if df is not None:
            print("Predictions DataFrame:")
            print(df.head())
            self.save_to_msp(df, output_file)
        else:
            print("No data to save.")
