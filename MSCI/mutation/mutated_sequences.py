import re
import itertools
from Bio import SeqIO
import pandas as pd


class ProteinMutator:
    def __init__(self, proteome_file, mutations_file, output_dir, digestion_method):
        self.proteome_file = proteome_file
        self.mutations_file = mutations_file
        self.output_dir = output_dir
        self.digestion_method = digestion_method
        self.proteome = {}
        self.mutations_df = None

    def load_proteome(self):
        """Load the proteome from a FASTA file."""
        accession_pattern = r'sp\|([A-Z0-9]+)\|'
        for record in SeqIO.parse(self.proteome_file, "fasta"):
            match = re.search(accession_pattern, record.id)
            if match:
                accession = match.group(1)
                self.proteome[accession] = str(record.seq)

    def load_mutations(self):
        """Load the mutations data from a TSV file."""
        self.mutations_df = pd.read_csv(self.mutations_file, sep="\t")

    def process_protein(self, target_protein_accession):
        """Process a single protein by generating peptides and mutating them based on mutations data."""
        print(f"\nProcessing protein accession: {target_protein_accession}")

        if target_protein_accession in self.proteome:
            protein_sequence = self.proteome[target_protein_accession]
            peptides = self.digestion_method(protein_sequence)  # Use the provided digestion method

            # Save the original peptides to a file
            output_filename = f"{self.output_dir}{target_protein_accession}_peptides.txt"
            with open(output_filename, 'w') as output_file:
                for original_peptide in peptides:
                    output_file.write(original_peptide + '\n')

            print(f"Original peptides saved to: {output_filename}")

            protein_mutations = self.mutations_df[self.mutations_df['Entry'] == target_protein_accession]['Natural variant'].tolist()

            if not pd.isna(protein_mutations[0]):
                total_mutations = self.mutations_df[self.mutations_df['Entry'] == target_protein_accession]['Natural variant'].count()
                print(f"Total number of mutations in the input data: {total_mutations}")

                mutations_per_peptide = {peptide: [] for peptide in peptides}

                for mutation_str in protein_mutations:
                    if not pd.isna(mutation_str):
                        mutations = re.findall(r'VARIANT (\d+); /note="([^"]+)"', mutation_str)

                        for variant, note in mutations:
                            position = int(variant)

                            if position > 0:
                                assigned = False
                                for peptide in peptides:
                                    peptide_start = protein_sequence.index(peptide) + 1
                                    peptide_end = peptide_start + len(peptide) - 1
                                    if peptide_start <= position <= peptide_end:
                                        existing_positions = [pos for pos, _ in mutations_per_peptide[peptide]]
                                        if position not in existing_positions:
                                            mutation_position_in_peptide = position - peptide_start
                                            mutations_per_peptide[peptide].append((mutation_position_in_peptide, note))
                                            assigned = True
                                            break

                                if not assigned:
                                    print(f"Warning: Mutation at position {position} could not be assigned to any peptide.")

                processed_mutations_count = sum(len(mutations) for mutations in mutations_per_peptide.values())
                print(f"Number of mutations after processing: {processed_mutations_count}")

                print("Original Mutations:")
                for mutation_str in protein_mutations:
                    if not pd.isna(mutation_str):
                        mutations = re.findall(r'VARIANT (\d+); /note="([^"]+)"', mutation_str)
                        for variant, note in mutations:
                            print(f"Position: {variant}, Note: {note}")

                print("\nProcessed Mutations:")
                for peptide, mutations in mutations_per_peptide.items():
                    for position, note in mutations:
                        print(f"Peptide: {peptide}, Position: {position}, Note: {note}")

                mutated_peptides = {}

                for peptide, mutations in mutations_per_peptide.items():
                    if len(mutations) > 20:
                        mutated_peptides[peptide] = []
                    else:
                        mutated_peptides[peptide] = []

                        for r in range(1, len(mutations) + 1):
                            for combo in itertools.combinations(mutations, r):
                                mutated_peptide = list(peptide)
                                for position, mutation_note in combo:
                                    mutation_position_in_peptide = position
                                    if 0 <= mutation_position_in_peptide < len(mutated_peptide):
                                        mutated_residue = mutation_note.split(' -> ')[-1].split(" ")[0]
                                        mutated_peptide[mutation_position_in_peptide] = mutated_residue
                                    else:
                                        print(f"Warning: Mutation at position {position} is out of range for the peptide {peptide}. Skipping.")

                                mutated_peptide = ''.join(mutated_peptide)

                                if mutated_peptide[-1] == 'K' or mutated_peptide[-1] == 'R':
                                    mutated_peptides[peptide].append(mutated_peptide)
                                else:
                                    next_peptide_index = peptides.index(peptide) + 1
                                    if next_peptide_index < len(peptides):
                                        next_peptide = peptides[next_peptide_index]
                                        combined_peptide = mutated_peptide + next_peptide
                                        mutated_peptides[peptide].append(combined_peptide)

                output_filename = f"{self.output_dir}{target_protein_accession}_mutated_peptides.txt"
                with open(output_filename, 'w') as output_file:
                    for original_peptide, mutations in mutated_peptides.items():
                        for mutated_peptide in mutations:
                            output_file.write(mutated_peptide + '\n')

                print(f"Mutated peptides saved to: {output_filename}")

            else:
                print(f"Processing protein accession {target_protein_accession} with NaN mutations.")
        else:
            print(f"Target protein accession {target_protein_accession} not found in the proteome.")

    def process_all_proteins(self):
        """Process all proteins in the proteome."""
        self.load_proteome()
        self.load_mutations()

        for target_protein_accession in self.proteome.keys():
            self.process_protein(target_protein_accession)


def tryptic_digest(sequence):
    peptides = []
    aa0 = '\0'
    digest = ""

    for aa1 in sequence:
        if aa0 != '\0':
            digest += aa0

        if (aa1 != 'P' and aa0 == 'R') or aa0 == 'K':
            if digest:
                peptides.append(digest)
            digest = ""

        aa0 = aa1

    if aa0 != '\0':
        digest += aa0

    if digest:
        peptides.append(digest)

    return peptides


proteome_file = "Z:/zelhamraoui/MSCA_Package/mutation/uniprotkb_Human_AND_reviewed_true_AND_m_2023_09_12.fasta"
mutations_file = "Z:/zelhamraoui/MSCA_Package/mutation/uniprotkb_Human_AND_reviewed_true_AND_m_2023_09_12.tsv"
output_dir = "Z:/zelhamraoui/MSCA_Package/mutation/Dataset/one_point_mutation/"

mutator = ProteinMutator(proteome_file, mutations_file, output_dir, tryptic_digest)
mutator.process_all_proteins()
