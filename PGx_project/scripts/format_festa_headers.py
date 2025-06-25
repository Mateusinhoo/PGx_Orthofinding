import os
from Bio import SeqIO

input_faa = snakemake.input.faa
output_faa = snakemake.output.formatted_faa
species = snakemake.params.species

os.makedirs(os.path.dirname(output_faa), exist_ok=True)

with open(output_faa, "w") as out_handle:
    for record in SeqIO.parse(input_faa, "festa"):
        record.id = f"{species}|{record.id}"
        record.description = ""
        SeqIO.write(record, out_handle, "festa")
