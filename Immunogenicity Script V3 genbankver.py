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

sequences = ["TTGT", "TTTC", "TGTT", "CTGT", "TATT", "TTTA", "TTGC", "GTTT", "ATTT", "ATGT", "CTTT", "TCTT"]

input_folder = r"C:\Users\iulia\PycharmProjects\pythonProject"
output_folder = r"C:\Users\iulia\PycharmProjects\pythonProject"
os.makedirs(output_folder, exist_ok=True)

# Identify GenBank files (.gb or .gbk)
genbank_files = [f for f in os.listdir(input_folder) if f.endswith('.gb') or f.endswith('.gbk')]

if not genbank_files:
    print("No GenBank files found in the input folder.")
else:
    for gb_file_name in genbank_files:
        gb_file_path = os.path.join(input_folder, gb_file_name)
        records = list(SeqIO.parse(gb_file_path, "genbank"))

        all_initial_data = []

        summary_data = []

        for record in records:
            gene = str(record.seq)

            gene_length, total_motif_count, total_percentage = count_sequences(gene, sequences)
            cpg_count, cpg_percentage = count_cpg_motifs(gene)

            initial_data = {
                'Gene': record.id,
                'Total Motif Count': total_motif_count,
                'Total Length (nucleotides)': gene_length,
                'Total Percentage': total_percentage,
                'CpG Motif Count': cpg_count,
                'CpG Percentage': cpg_percentage
            }
            all_initial_data.append(initial_data)

        # Create scatter plots
        total_counts = [data['Total Motif Count'] for data in all_initial_data]
        total_counts_chart_output_path = os.path.join(output_folder,
                                                      f"{os.path.splitext(gb_file_name)[0]}_total_motif_counts_scatter.png")
        create_scatter_plot(total_counts, "Total Motif Counts vs Genes", "Total Motif Count",
                            total_counts_chart_output_path)

        total_percentages = [data['Total Percentage'] for data in all_initial_data]
        total_percentages_chart_output_path = os.path.join(output_folder,
                                                           f"{os.path.splitext(gb_file_name)[0]}_total_percentage_scatter.png")
        create_scatter_plot(total_percentages, "Total Percentage vs Genes", "Total Percentage",
                            total_percentages_chart_output_path)

        output_path = os.path.join(output_folder, "immunostimulatory_motif.xlsx")

        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            pd.DataFrame(all_initial_data).to_excel(writer, sheet_name='Initial Analysis', index=False)

            workbook = writer.book
            charts_sheet = workbook.create_sheet(title="Charts")

            row = 1
            img = Image(total_counts_chart_output_path)
            charts_sheet.add_image(img, f"A{row}")
            row += 20

            img = Image(total_percentages_chart_output_path)
            charts_sheet.add_image(img, f"A{row}")

    print(f"Analysis results have been saved to {output_folder}")

# Save the script
script_content = """
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

sequences = ["TTGT", "TTTC", "TGTT", "CTGT", "TATT", "TTTA", "TTGC", "GTTT", "ATTT", "ATGT", "CTTT", "TCTT"]

input_folder = "C:\\Users\\Jamie Mann\\Desktop\\Python submission folder"
output_folder = "C:\\Users\\Jamie Mann\\Desktop\\python excel folder"
os.makedirs(output_folder, exist_ok=True)

fasta_files = [f for f in os.listdir(input_folder) if f.endswith('.fasta')]

if not fasta_files:
    print("No FASTA files found in the input folder.")
else:
    for fasta_file_name in fasta_files:
        fasta_file_path = os.path.join(input_folder, fasta_file_name)
        records = list(SeqIO.parse(fasta_file_path, "fasta"))

        all_initial_data = []
        summary_data = []

        for record in records:
            gene = str(record.seq)

            gene_length, total_motif_count, total_percentage = count_sequences(gene, sequences)
            cpg_count, cpg_percentage = count_cpg_motifs(gene)

            initial_data = {
                'Gene': record.id,
                'Total Motif Count': total_motif_count,
                'Total Length (nucleotides)': gene_length,
                'Total Percentage': total_percentage,
                'CpG Motif Count': cpg_count,
                'CpG Percentage': cpg_percentage
            }
            all_initial_data.append(initial_data)

        # Create scatter plots
        total_counts = [data['Total Motif Count'] for data in all_initial_data]
        total_counts_chart_output_path = os.path.join(output_folder,
                                                      f"{os.path.splitext(fasta_file_name)[0]}_total_motif_counts_scatter.png")
        create_scatter_plot(total_counts, "Total Motif Counts vs Genes", "Total Motif Count",
                            total_counts_chart_output_path)

        total_percentages = [data['Total Percentage'] for data in all_initial_data]
        total_percentages_chart_output_path = os.path.join(output_folder,
                                                           f"{os.path.splitext(fasta_file_name)[0]}_total_percentage_scatter.png")
        create_scatter_plot(total_percentages, "Total Percentage vs Genes", "Total Percentage",
                            total_percentages_chart_output_path)

        output_path = os.path.join(output_folder, "immunostimulatory_motif.xlsx")

        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            pd.DataFrame(all_initial_data).to_excel(writer, sheet_name='Initial Analysis', index=False)

            workbook = writer.book
            charts_sheet = workbook.create_sheet(title="Charts")

            row = 1
            img = Image(total_counts_chart_output_path)
            charts_sheet.add_image(img, f"A{row}")
            row += 20

            img = Image(total_percentages_chart_output_path)
            charts_sheet.add_image(img, f"A{row}")

    print(f"Analysis results have been saved to {output_folder}")
"""

script_path = 'C:\\Users\\iulia\\PycharmProjects\\pythonProject\\output_script.py'
with open(script_path, 'w') as script_file:
    # Write content to the file
    script_file.write("Your script content here.")
    print(f"Script has been saved to {script_path}")  # Properly indented
