import re
import os

#Make sure you include the , after the variable or it'll be saved in an incorrect format
DBS, = glob_wildcards("databases/{database}.prepareref")
MLSTs, = glob_wildcards("databases/{MLST}.mlstdb")
sample_ids, = glob_wildcards("reads/{sample}_R1.fastq.gz")

#checks your files are being found
print(sample_ids,)
print(DBS,)
print(MLSTs,)

#rule all is needed to direct the outputs of other rules
rule all:
    input:
        expand(directory("ariba_out/{sample}.{MLST}.MLST_out"), sample=sample_ids, MLST=MLSTs),
        expand(directory("ariba_out/{sample}.{database}.out"), sample=sample_ids, database=DBS),
        expand("ariba_summaries/{MLST}_MLST.tsv", MLST=MLSTs)

#
rule ariba:    
    input:
        db = "databases/{database}.prepareref",
        r1 = "reads/{sample}_R1.fastq.gz",
        r2 = "reads/{sample}_R2.fastq.gz"
    output:
        directory("ariba_out/{sample}.{database}.out")
    threads: 80
    shell:
        "ariba run --threads 80 --noclean {input.db} {input.r1} {input.r2} {output}"
