import os
from snakemake.utils import min_version

# Set minimum Snakemake version
min_version("6.0")

# Load configuration
configfile: "config.json"

# Define wildcards
SAMPLES = [os.path.basename(f).replace("_R1.fastq.gz", "") for f in config["reads"]]

# Define final output
rule all:
    input:
        "PROFILE.db",
        "splits_table.tsv"

# Index reference
rule index_reference:
    input:
        ref = config["reference"]
    output:
        bwt = config["reference"] + ".bwt",
        other = multiext(config["reference"], ".pac", ".ann", ".amb", ".sa")
    threads: config["threads"]["bwa_index"]
    shell:
        "mkdir -p logs && bwa index {input.ref} > logs/bwa_index.log 2> logs/bwa_index.err"

# Create BAM files
rule create_bam:
    input:
        ref = config["reference"],
        r1 = "reads/{sample}_R1.fastq.gz",
        r2 = "reads/{sample}_R2.fastq.gz",
        bwt = rules.index_reference.output.bwt
    output:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.csi"
    threads: config["threads"]["bwa_mem"]
    shell:
        "mkdir -p bam && "
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2}  2> logs/bwa_{wildcards.sample}.err | "
        "samtools view -bS -F4 - | "
        "samtools sort -o {output.bam} --write-index - > logs/sort_{wildcards.sample}.log 2> logs/sort_{wildcards.sample}.err"

# Create CONTIGS.db
rule create_contigs_db:
    input:
        ref = config["reference"]
    output:
        db = "CONTIGS.db"
    threads: config["threads"]["anvi_contigs_db"]
    shell:
        "mkdir -p logs && anvi-gen-contigs-database -f {input.ref} -o {output.db} 2> logs/anvi-gen.log > logs/anvi-gen.out  && "
        "anvi-run-hmms -c {output.db} -T {threads}  2> logs/anvi-hmms.log  > logs/anvi-hmms.out && "
        "anvi-run-ncbi-cogs -c {output.db} -T {threads}  2> logs/anvi-cogs.log > logs/anvi-cogs.out  && "
        "anvi-scan-trnas -c {output.db}  2> logs/anvi-trnas.log > logs/anvi-trnas.out && "
        "anvi-run-scg-taxonomy -c {output.db}  2> logs/anvi-tax.log > logs/anvi-tax.out "

# Create profile directories
rule create_profile:
    input:
        bam = "bam/{sample}.bam",
        db = rules.create_contigs_db.output.db
    output:
        profile = directory("bam/{sample}")
    threads: config["threads"]["anvi_profile"]
    shell:
        "anvi-profile -i {input.bam} -c {input.db} -T {threads} -o {output.profile}  > logs/anvi-profile.out 2> logs/anvi-profile.err"

# Export tables
rule export_tables:
    input:
        db = rules.create_contigs_db.output.db
    output:
        ctgtab = "contigs_table.tsv",
        spttab = "splits_table.tsv"
    shell:
        "anvi-export-table {input.db}  --table splits_basic_info -o splits_table.tsv >> logs/anvi-export.out 2>> logs/anvi-export.err && "
        "anvi-export-table {input.db}  --table genes_in_splits -o genes_splits_table.tsv >> logs/anvi-export.out 2>> logs/anvi-export.err && "
        "anvi-export-table {input.db}  --table genes_in_contigs -o genes_contigs.tsv  >> logs/anvi-export.out 2>> logs/anvi-export.err && "
        "anvi-export-table {input.db}  --table contig_sequences -o contigs_table.tsv  >> logs/anvi-export.out 2>> logs/anvi-export.err"
        



# Merge profiles
rule merge_profiles:
    input:
        profiles = expand("bam/{sample}", sample=SAMPLES),
        db = rules.create_contigs_db.output.db
    output:
        merged = "PROFILE.db"
    params:
        profile_paths = lambda wildcards, input: " ".join([d + "/PROFILE.db" for d in input.profiles])
    threads: config["threads"]["anvi_merge"]
    shell:
        "anvi-merge {params.profile_paths} -o ./MINITEST -c {input.db} > logs/anvi-merge.out 2> logs/anvi-merge.err && "
        "mv MINITEST/*.db . && "
        "rm -rf MINITEST"



print("When finished: ")
print("anvi-interactive -d extras/metadata.txt")