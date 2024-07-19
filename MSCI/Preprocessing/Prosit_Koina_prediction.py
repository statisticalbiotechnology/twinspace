import requests
import pandas as pd

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
            predictions = response.json()
            try:
                intensities = predictions['outputs'][2]['data']
                mz_values = predictions['outputs'][1]['data']
            except KeyError:
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
            return pd.DataFrame(rows)
        else:
            return None

    def save_to_msp(self, df, file_path):
        with open(file_path, 'w') as file:
            for _, row in df.iterrows():
                mw = self.calculate_peptide_mass(row['peptide_sequence'])
                mz = (mw + (self.charge * 1.007276)) / self.charge
                header = f"Name: {row['peptide_sequence']}/{self.charge}\nMW: {mz:.6f}\nCollision_energy: {self.ce}\n"
                peaks = "".join(f"{mz:.2f}\t{intensity:.6f}\n" for mz, intensity in zip(row['mz_values'], row['intensities']) if mz != -1.00 and intensity != -1.000000)
                file.write(f"{header}Num peaks: {len(peaks)}\n{peaks}\n")

    def process(self):
        df = self.get_prosit_predictions(self.peptides)
        if df is not None:
            self.save_to_msp(df, 'output.msp')
