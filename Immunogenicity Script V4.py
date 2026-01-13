import pandas as pd
from Bio import SeqIO
import os
import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl.drawing.image import Image

def count_sequences(gene, sequences):
    sequence_counts = {seq: 0 for seq in sequences}
    gene_length = len(gene)

    for i in range(gene_length - 3):
        substring = gene[i:i + 4]
        if substring in sequence_counts:
            sequence_counts[substring] += 1

    total_occurrences = sum(sequence_counts.values())
    total_percentage = (total_occurrences / gene_length) * 100
    total_motif_count = sum(sequence_counts.values())

    return gene_length, total_motif_count, total_percentage

def count_cpg_motifs(gene):
    cpg_count = gene.count("CG")
    gene_length = len(gene)
    cpg_percentage = (cpg_count / gene_length) * 100 if gene_length > 0 else 0
    return cpg_count, cpg_percentage

def create_scatter_plot(y_data, title, ylabel, output_path):
    plt.figure(figsize=(10, 6))
    plt.scatter(range(len(y_data)), y_data, color='skyblue')
    plt.title(title)
    plt.xlabel('Gene Index')
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# List of sequences to search for
sequences = ["TTGT", "TTTC", "TGTT", "CTGT", "TATT", "TTTA", "TTGC", "GTTT", "ATTT", "ATGT", "CTTT", "TCTT"]

input_folder = r"C:\Users\iulia\PycharmProjects\pythonProject"
output_folder = r"C:\Users\iulia\PycharmProjects\pythonProject"
os.makedirs(output_folder, exist_ok=True)

# Modify to handle both .fasta and .fas files
fasta_files = [f for f in os.listdir(input_folder) if f.endswith('.fasta') or f.endswith('.fas')]

if not fasta_files:
    print("No FASTA files found in the input folder.")
else:
    # Prepare to append data from all files into a single Excel sheet
    all_data = [] 

    # Output path for the Excel file
    output_path = os.path.join(output_folder, "combined_immunostimulatory_motif_analysis.xlsx")

    # Initialize Excel writer
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        for fasta_file_name in fasta_files:
            fasta_file_path = os.path.join(input_folder, fasta_file_name)
            records = list(SeqIO.parse(fasta_file_path, "fasta"))

            file_data = []

            for record in records:
                gene = str(record.seq)

                gene_length, total_motif_count, total_percentage = count_sequences(gene, sequences)
                cpg_count, cpg_percentage = count_cpg_motifs(gene)

                data = {
                    'File': fasta_file_name,
                    'Gene': record.id,
                    'Total Motif Count': total_motif_count,
                    'Total Length (nucleotides)': gene_length,
                    'Total Percentage': total_percentage,
                    'CpG Motif Count': cpg_count,
                    'CpG Percentage': cpg_percentage
                }
                file_data.append(data)
                all_data.append(data)

            # Create scatter plots for each file
            total_counts = [data['Total Motif Count'] for data in file_data]
            total_counts_chart_output_path = os.path.join(output_folder,
                                                          f"{os.path.splitext(fasta_file_name)[0]}_total_motif_counts_scatter.png")
            create_scatter_plot(total_counts, "Total Motif Counts vs Genes", "Total Motif Count",
                                total_counts_chart_output_path)

            total_percentages = [data['Total Percentage'] for data in file_data]
            total_percentages_chart_output_path = os.path.join(output_folder,
                                                               f"{os.path.splitext(fasta_file_name)[0]}_total_percentage_scatter.png")
            create_scatter_plot(total_percentages, "Total Percentage vs Genes", "Total Percentage",
                                total_percentages_chart_output_path)

            # Write each file's data to a separate sheet in the Excel file
            pd.DataFrame(file_data).to_excel(writer, sheet_name=f'{os.path.splitext(fasta_file_name)[0]} Analysis', index=False)

            # Add charts to Excel
            workbook = writer.book
            charts_sheet = workbook.create_sheet(title=f"{os.path.splitext(fasta_file_name)[0]} Charts")

            row = 1
            img = Image(total_counts_chart_output_path)
            charts_sheet.add_image(img, f"A{row}")
            row += 20

            img = Image(total_percentages_chart_output_path)
            charts_sheet.add_image(img, f"A{row}")

        # After processing all files, append the combined data into one sheet
        pd.DataFrame(all_data).to_excel(writer, sheet_name='Combined Analysis', index=False)

    print(f"Analysis results have been saved to {output_path}")
