import re
import itertools
from Bio import SeqIO
import pandas as pd

# Load the proteome
accession_pattern = r'sp\|([A-Z0-9]+)\|'
proteome = {}
for record in SeqIO.parse("Z:/zelhamraoui/MSCA_Package/mutation/uniprotkb_Human_AND_reviewed_true_AND_m_2023_09_12.fasta", "fasta"):
    match = re.search(accession_pattern, record.id)
    if match:
        accession = match.group(1)
        proteome[accession] = str(record.seq)

# Load the mutations data
mutations_df = pd.read_csv("Z:/zelhamraoui/MSCA_Package/mutation/uniprotkb_Human_AND_reviewed_true_AND_m_2023_09_12.tsv", sep="\t")

# Output directory
output_dir = "Z:/zelhamraoui/MSCA_Package/mutation/Dataset/one_point_mutation/"

# Iterate through all protein accessions in the proteome
for target_protein_accession in proteome.keys():
    print(f"\nProcessing protein accession: {target_protein_accession}")

    # Check if the target protein accession is in the proteome
    if target_protein_accession in proteome:
        # Extract the protein sequence
        protein_sequence = proteome[target_protein_accession]

        # Digest the protein sequence into peptides
        def tryptic_digest(sequence):
            peptides = []
            aa0 = '\0'  # Initialize the previous amino acid as the null character
            digest = ""

            for aa1 in sequence:
                if aa0 != '\0':
                    digest += aa0

                if (aa1 != 'P' and aa0 == 'R') or aa0 == 'K':
                    if digest:
                        peptides.append(digest)
                    digest = ""

                aa0 = aa1

            # This is the end of the sequence
            if aa0 != '\0':
                digest += aa0

            if digest:
                peptides.append(digest)

            return peptides

        peptides = tryptic_digest(protein_sequence)

        # Save the original peptides to a file
        output_filename = f"{output_dir}{target_protein_accession}_peptides.txt"
        with open(output_filename, 'w') as output_file:
            for original_peptide in peptides:
                output_file.write(original_peptide + '\n')

        print(f"Original peptides saved to: {output_filename}")

        # Extract mutations for the target protein from mutations_df
        protein_mutations = mutations_df[mutations_df['Entry'] == target_protein_accession]['Natural variant'].tolist()

        # Skip if mutations are NaN
        if not pd.isna(protein_mutations[0]):
            # Count the total number of mutations in the mutations_df
            total_mutations = mutations_df[mutations_df['Entry'] == target_protein_accession]['Natural variant'].count()
            print(f"Total number of mutations in the input data: {total_mutations}")

            # Initialize a dictionary to store mutations per peptide
            mutations_per_peptide = {peptide: [] for peptide in peptides}

            # Map mutations to peptides, ensuring no duplication and position > 0
            for mutation_str in protein_mutations:
                if not pd.isna(mutation_str):  # Skip if mutation_str is NaN
                    mutations = re.findall(r'VARIANT (\d+); /note="([^"]+)"', mutation_str)

                    for variant, note in mutations:
                        position = int(variant)

                        if position > 0:  # Ensure position is greater than 0
                            assigned = False  # Flag to track if the mutation is assigned to any peptide
                            for peptide in peptides:
                                peptide_start = protein_sequence.index(peptide) + 1  # Adjust index to start from 1
                                peptide_end = peptide_start + len(peptide) - 1  # Adjust index to start from 1
                                if peptide_start <= position <= peptide_end:
                                    # Check if the mutation position is not in the list already
                                    existing_positions = [pos for pos, _ in mutations_per_peptide[peptide]]
                                    if position not in existing_positions:
                                        # Adjust the mutation position to be relative to the peptide
                                        mutation_position_in_peptide = position - peptide_start
                                        mutations_per_peptide[peptide].append((mutation_position_in_peptide, note))
                                        assigned = True
                                        break  # Break the loop once the mutation is assigned to a peptide

                            if not assigned:
                                print(f"Warning: Mutation at position {position} could not be assigned to any peptide.")

            # Count the number of mutations after processing
            processed_mutations_count = sum(len(mutations) for mutations in mutations_per_peptide.values())
            print(f"Number of mutations after processing: {processed_mutations_count}")

            # Print the original and processed mutations
            print("Original Mutations:")
            for mutation_str in protein_mutations:
                if not pd.isna(mutation_str):  # Skip if mutation_str is NaN
                    mutations = re.findall(r'VARIANT (\d+); /note="([^"]+)"', mutation_str)
                    for variant, note in mutations:
                        print(f"Position: {variant}, Note: {note}")

            print("\nProcessed Mutations:")
            for peptide, mutations in mutations_per_peptide.items():
                for position, note in mutations:
                    print(f"Peptide: {peptide}, Position: {position}, Note: {note}")

            # Initialize a dictionary to store mutated peptides
            mutated_peptides = {}

            # Iterate through each peptide
            for peptide, mutations in mutations_per_peptide.items():
                if len(mutations) > 20:
                    mutated_peptides[peptide] = []
                else:
                    mutated_peptides[peptide] = []

                    # Generate all possible combinations of mutations within the peptide
                    for r in range(1, len(mutations) + 1):
                        for combo in itertools.combinations(mutations, r):
                            mutated_peptide = list(peptide)  # Convert the peptide to a list to modify it
                            for position, mutation_note in combo:
                                # Apply mutations (mutation_position_in_peptide to adjust for 0-based indexing)
                                mutation_position_in_peptide = position
                                if 0 <= mutation_position_in_peptide < len(mutated_peptide):
                                    mutated_residue = mutation_note.split(' -> ')[-1].split(" ")[0]
                                    mutated_peptide[mutation_position_in_peptide] = mutated_residue
                                else:
                                    print(f"Warning: Mutation at position {position} is out of range for the peptide {peptide}. Skipping.")

                            # Join the mutated peptide list back to a string
                            mutated_peptide = ''.join(mutated_peptide)

                            # Check if the mutated peptide is still tryptic
                            if mutated_peptide[-1] == 'K' or mutated_peptide[-1] == 'R':
                                mutated_peptides[peptide].append(mutated_peptide)
                            else:
                                # Combine with the next tryptic peptide
                                next_peptide_index = peptides.index(peptide) + 1
                                if next_peptide_index < len(peptides):
                                    next_peptide = peptides[next_peptide_index]
                                    combined_peptide = mutated_peptide + next_peptide
                                    mutated_peptides[peptide].append(combined_peptide)

            # Save the output to a file
            output_filename = f"{output_dir}{target_protein_accession}_peptides.txt"
            with open(output_filename, 'w') as output_file:
                for original_peptide, mutations in mutated_peptides.items():
                    output_file.write(original_peptide + '\n')

            print(f"Output saved to: {output_filename}")

        else:
            print(f"Processing protein accession {target_protein_accession} with NaN mutations.")

    else:
        print(f"Target protein accession {target_protein_accession} not found in the proteome.")
