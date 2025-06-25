configfile: "config.yaml"

SPECIES = config["species"]
ORTHOFINDER_OUTDIR = config["orthofinder_outdir"]
THREADS = config.get("threads", 4)

rule all:
    input:
        cds=expand("species_cds/{species}_cds.zip", species=SPECIES),
        fna=expand("fasta_output/{species}.fna", species=SPECIES),
        faa=expand("protein_output/{species}.faa", species=SPECIES),
        formatted_faa=expand("orthofinder_input/{species}.faa",species=SPECIES),
        orthogroups=f"{ORTHOFINDER_OUTDIR}/Results/Orthogroups/Orthogroups.tsv"

rule download_cds:
    output:
        cds="species_cds/{species}_cds.zip"
    conda:
        "envs/ncbi_env.yaml"
    shell:
        """
        datasets download genome taxon {wildcards.species} --reference --include cds --filename {output.cds}
        """

rule unzip_and_translate:
    input:
        zipfile="species_cds/{species}_cds.zip"
    output:
        fna= "fasta_output/{species}.fna",
        faa= "protein_output/{species}.faa"
    params:
        species="{species}"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/unzip_fasta.py"

rule format_fasta_headers:
    input:
        faa="protein_output/{species}.faa"
    output:
        formatted_faa="orthofinder_input/{species}.faa"
    params:
        species="{species}"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/unzip_fasta.py"

rule install_orthofinder:
    output:
        dir="tools/OrthoFinder/orthofinder"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        mkdir -p tools/OrthoFinder
        wget -q --show-progress https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz
        tar -xzf OrthoFinder.tar.gz -C tools/OrthoFinder --strip-components=1
        rm OrthoFinder.tar.gz
        """
rule run_orthofinder:
    input:
        dir="tools/OrthoFinder/orthofinder",
        faa=expand("orthofinder_input/{species}.faa", species=config["species"])
    output:
        orthogroups=f"{config['orthofinder_outdir']}/Results/Orthogroups/Orthogroups.tsv"
    threads: config.get("threads", 4)
    log:
        "logs/orthofinder.log"
    shell:
        """
        set -euo pipefail
        cd {input.dir}
        ./orthofinder -f {os.path.abspath('orthofinder_input')} -t {threads} -a {threads} -og -o {os.path.abspath(config['orthofinder_outdir'])} &> {log}
        """
