import re
import os

configfile:
        "config/config.yaml"

sample_ids, = glob_wildcards(config['assembly_path']+"/{sample}.fasta")
dbs, = glob_wildcards(config['ref_dbs']+"/{db}.fa")

print(sample_ids)
print(dbs)

rule all:
        input:
                expand(config['ref_dbs']+"/{db}.fa.nhr", db=dbs),
                expand("blast/{db}/{sample}.csv", db=dbs, sample=sample_ids),
                expand("blast/{db}/{sample}.new.tsv", db=dbs, sample=sample_ids),

rule makeblastdb:
        input:
                config['ref_dbs']+"/{db}.fa"
        output:
                config['ref_dbs']+"/{db}.fa.nhr"
        log:
                "blast/logs/{db}.log"
        params:
                blast_bin = config['blast_bin']
        shell:
                """
                {params.blast_bin}/makeblastdb -dbtype nucl -in {input}
                """

rule blastn:
        input:
                db = config['ref_dbs']+"/{db}.fa",
                assemblies = config['assembly_path']+"/{sample}.fasta"
        output:
                "blast/{db}/{sample}.csv"
        log:
                "blast/{db}/{sample}.log"
        params:
                blast_bin = config['blast_bin']
        shell:
                """
                {params.blast_bin}/blastn -num_threads 2 -evalue 0.001 -db {input.db} -query {input.assemblies} -out {output} -outfmt "6 qseqid stitle sseqid pident length slen sstart send qstart qend qlen mismatch gapopen evalue bitscore"
                """

rule name_cat:
        input:
                "blast/{db}/{sample}.csv"
        output:
                "blast/{db}/{sample}.new.tsv"
        priority: 8
        shell:
                """
                for f in {input} ; do sed 's@$@\\t{wildcards.sample}@' ${{f}} > {output} ; done
                cat {output} >> {wildcards.db}.txt
                """
