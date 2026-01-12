import os
import random
import pandas as pd
from Bio import SeqIO
from copy import deepcopy

# Define the directory paths
input_dir = r"C:\Users\iulia\PycharmProjects\pythonProject"
output_dir = r"C:\Users\iulia\PycharmProjects\pythonProject"

# Load your fasta file
fasta_files = [f for f in os.listdir(input_dir) if f.endswith('.fasta') or f.endswith('.fas')]
if not fasta_files:
    raise FileNotFoundError("No FASTA file found in input directory!")

fasta_path = os.path.join(input_dir, fasta_files[0])
record = SeqIO.read(fasta_path, "fasta")
dna_sequence = str(record.seq).upper()

# Define mutation function
def mutate_leucines_from_original(dna_sequence, output_dir, run_number):
    # Break into codons
    codons_original = [dna_sequence[i:i + 3] for i in range(0, len(dna_sequence), 3)]

    # Define sets
    leucine_codons = {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
    good_leucines = {"CTA"}
    bad_leucines = {"TTA", "TTG"}

    # Get indices of bad leucines in the original sequence
    bad_indices = [i for i, codon in enumerate(codons_original) if codon in bad_leucines]

    if not bad_indices:
        print("No bad leucines to mutate.")
        return

    variants = []
    variant_labels = []
    mutation_logs = []

    # Save original
    variants.append(''.join(codons_original))
    variant_labels.append("Original")
    mutation_logs.append("None")

    for num_mutations in range(1, len(bad_indices) + 1):
        # Create a fresh copy from original each time
        codons_copy = deepcopy(codons_original)

        # Randomly sample `num_mutations` unique bad indices
        indices_to_mutate = random.sample(bad_indices, num_mutations)

        mutation_details = []
        for idx in indices_to_mutate:
            old_codon = codons_copy[idx]
            codons_copy[idx] = "CTA"
            mutation_details.append(f"Codon {idx + 1}: {old_codon} → CTA")

        new_sequence = ''.join(codons_copy)
        variants.append(new_sequence)
        variant_labels.append(f"Variant {num_mutations}")
        mutation_logs.append('; '.join(mutation_details))

    # Save to Excel
    excel_path = os.path.join(output_dir, f'env_leucine_variants_run_{run_number}.xlsx')
    df = pd.DataFrame({'Label': variant_labels, 'Sequence': variants, 'Mutation': mutation_logs})
    df.to_excel(excel_path, index=False)
    print(f"Run {run_number}: Variants saved to {excel_path}.")

    # Save to FASTA
    fasta_path = os.path.join(output_dir, f'env_leucine_variants_run_{run_number}.fasta')
    with open(fasta_path, 'w') as f:
        for label, sequence in zip(variant_labels, variants):
            f.write(f">{label}\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + '\n')
    print(f"Run {run_number}: FASTA saved to {fasta_path}.")

# Run the mutation process 10 times
for run in range(1, 11):
    mutate_leucines_from_original(dna_sequence, output_dir, run)


