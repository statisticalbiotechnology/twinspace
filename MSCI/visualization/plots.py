
# Create a new DataFrame to store the sequences corresponding to x_id and y_id
sequence_mapping = pd.DataFrame(columns=['x_id', 'y_id', 'x_sequence', 'y_sequence'])

# Function to get the sequence for a given ID
def get_sequence(id):
    if id in df.index:
        return spectra[id].get("compound_name")
    else:
        return None

# Populate the new DataFrame with sequences
sequence_mapping['x_id'] = hrhr_filtered['x_id']
sequence_mapping['y_id'] = hrhr_filtered['y_id']
sequence_mapping['x_sequence'] = sequence_mapping['x_id'].apply(get_sequence)
sequence_mapping['y_sequence'] = sequence_mapping['y_id'].apply(get_sequence)

# Remove '/2' from 'x_sequence' and 'y_sequence'
sequence_mapping['x_sequence'] = sequence_mapping['x_sequence'].str.replace('/2', '')
sequence_mapping['y_sequence'] = sequence_mapping['y_sequence'].str.replace('/2', '')

# Display the new DataFrame
print(sequence_mapping)

# Function to compare sequences while ignoring 'I' and 'L'
def sequences_are_same(seq1, seq2):
    if len(seq1) != len(seq2):
        return False
    
    for char1, char2 in zip(seq1, seq2):
        if char1 != 'I' and char1 != 'L' and char1 != char2:
            return False
    
    return True

# Compare sequences and add a new column to the DataFrame
sequence_mapping['sequences_are_same'] = sequence_mapping.apply(lambda row: sequences_are_same(row['x_sequence'], row['y_sequence']), axis=1)

# Display the DataFrame
print(sequence_mapping)

# Count the number of True values in the 'sequences_are_same' column
true_count = sequence_mapping['sequences_are_same'].sum()

# Display the count
print("Number of sequences that are the same while ignoring 'I' and 'L':", true_count)

# Assuming you have already added the 'sequences_are_same' column to your DataFrame as shown earlier

# Filter the DataFrame to keep only the rows where 'sequences_are_same' is True
true_pairs = sequence_mapping[sequence_mapping['sequences_are_same']]

# Extract unique peptides from 'x_sequence' and 'y_sequence' separately
unique_x_peptides = set(true_pairs['x_sequence'])
unique_y_peptides = set(true_pairs['y_sequence'])

# Combine the unique peptides from both sequences
unique_peptides = unique_x_peptides.union(unique_y_peptides)

# Count the number of unique peptides
num_unique_peptides = len(unique_peptides)

# Display the count of unique peptides
print("Number of unique peptides in true pairs:", num_unique_peptides)
