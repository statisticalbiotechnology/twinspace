from Bio import SeqIO

class PeptideGenerator:
    def __init__(self, fasta_file, min_length=8, max_length=11):
        """
        Initialize the PeptideGenerator with a FASTA file and a range for peptide lengths.

        :param fasta_file: Path to the FASTA file containing protein sequences.
        :param min_length: Minimum length of peptides to generate (default is 8).
        :param max_length: Maximum length of peptides to generate (default is 11).
        """
        self.fasta_file = fasta_file
        self.min_length = min_length
        self.max_length = max_length
        self.all_peptides = []

    def generate_variable_length_peptides(self, protein_sequence):
        """
        Generate all variable length peptides from a given protein sequence.

        :param protein_sequence: The protein sequence to generate peptides from.
        :return: A list of peptides of varying lengths.
        """
        peptides = []
        for length in range(self.min_length, self.max_length + 1):
            peptides.extend([protein_sequence[i:i+length] for i in range(len(protein_sequence) - length + 1)])
        return peptides

    def parse_fasta_and_generate_peptides(self):
        """
        Parse the FASTA file and generate variable length peptides for each protein sequence.
        Stores all peptides in the all_peptides attribute.
        """
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            protein_sequence = str(record.seq)
            self.all_peptides.extend(self.generate_variable_length_peptides(protein_sequence))

    def save_peptides_to_file(self, output_file):
        """
        Save the generated peptides to a text file, one peptide per line.

        :param output_file: Path to the output file where peptides will be saved.
        """
        with open(output_file, 'w') as f:
            for peptide in self.all_peptides:
                f.write(peptide + '\n')

# Example usage:
fasta_file = "Z:/zelhamraoui/MSCA_Package/IMMUNO/immuno/sp_human_2023_04.fasta"
output_file = "Z:/zelhamraoui/MSCA_Package/IMMUNO/immuno/variable_length_peptides.txt"

# Create an instance of PeptideGenerator with desired peptide length range
peptide_generator = PeptideGenerator(fasta_file, min_length=8, max_length=11)

# Generate peptides
peptide_generator.parse_fasta_and_generate_peptides()

# Save peptides to a file
peptide_generator.save_peptides_to_file(output_file)
