#Import the needed packages
import zipfile
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#Snakemake inputs
zip_species = snakemake.input.zipfile
output_fna = snakemake.output.fna   #download the fasta files
output_faa = snakemake.output.faa    #convert to protein files
species_name = snakemake.params.species

#Make sure the output directory exists
os.makedirs(os.path.dirname(output_fna), exist_ok=True)
os.makedirs(os.path.dirname(output_faa), exist_ok=True)

#Step 1: Extract .fna file from ZIP
try:
    with zipfile.ZipFile(zip_species, 'r') as zp:
        fasta_files = [
            f for f in zp.namelist()
            if f.startswith('ncbi_dataset/data/') and f.endswith('.fna')
        ]

        if not fasta_files:
            print(f'[{species_name}] No FASTA files found in {zip_species}')
            exit(1)

        fasta_path = fasta_files[0]
        with zp.open(fasta_path) as src, open(output_fna, 'wb') as dst:
            dst.write(src.read())
        print(f"[{species_name}] Extracted {os.path.basename(fasta_path)}")

except Exception as e:
    print(f'[{species_name}] Failed to extract {zip_species}:{e}')


#Step 2: Translate .fna to .faa
try:
    with open(output_faa, "w") as faa_handle:
        for record in SeqIO.parse(output_fna, "fasta"):
            try:
                protein_seq = record.seq.translate(to_stop=True)
                protein_record = SeqRecord(protein_seq, id=record.id, description="" )
                SeqIO.write(protein_record, faa_handle, "fasta")
            except Exception as e:
                print(f"[{species_name}] Error translating {record.id}: {e}")
    print(f"[{species_name}] Translation complete -> {output_faa}")
except Exception as e:
    print(f"[{species_name}] Failed to write protein FASTA: {e}")
    exit(1)
