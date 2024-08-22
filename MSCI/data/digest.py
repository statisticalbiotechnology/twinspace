from Bio import SeqIO

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

def chymotrypsin_digest(sequence):
    # Example of a different digestion method
    peptides = []
    aa0 = '\0'
    digest = ""

    for aa1 in sequence:
        if aa0 != '\0':
            digest += aa0

        if aa0 in ['F', 'W', 'Y', 'L']:
            if digest:
                peptides.append(digest)
            digest = ""

        aa0 = aa1

    if aa0 != '\0':
        digest += aa0

    if digest:
        peptides.append(digest)

    return peptides

def parse_fasta_and_digest(fasta_file, digest_type="trypsin"):
    digestion_function = None
    
    # Choose the digestion function based on the digest_type argument
    if digest_type == "trypsin":
        digestion_function = tryptic_digest
    elif digest_type == "chymotrypsin":
        digestion_function = chymotrypsin_digest
    else:
        raise ValueError(f"Unknown digestion type: {digest_type}")

    # Parse the FASTA file
    results = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        peptides = digestion_function(sequence)
        results[record.id] = peptides

    return results
import csv

def peptides_to_csv(result, output_file):
    # Open the output CSV file
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # Write the header
        writer.writerow(["Peptide"])
        
        # Write the peptides
        for peptides in result.values():
            for peptide in peptides:
                writer.writerow([peptide])
