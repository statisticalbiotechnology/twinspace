from Bio import SeqIO

def generate_variable_length_peptides(protein_sequence, min_length=8, max_length=11):
    """Generate all variable length peptides from a protein sequence."""
    peptides = []
    for length in range(min_length, max_length + 1):
        peptides.extend([protein_sequence[i:i+length] for i in range(len(protein_sequence) - length + 1)])
    return peptides

# Path to the human proteome FASTA file
fasta_file = "Z:/zelhamraoui/MSCA_Package/IMMUNO/immuno/sp_human_2023_04.fasta"

# List to store all variable length peptides
all_peptides = []

# Parse the FASTA file and generate variable length peptides
for record in SeqIO.parse(fasta_file, "fasta"):
    protein_sequence = str(record.seq)
    all_peptides.extend(generate_variable_length_peptides(protein_sequence))
