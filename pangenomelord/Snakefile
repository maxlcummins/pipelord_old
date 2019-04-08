import re
import os

configfile:
	"config/config.yaml"

sample_ids, = glob_wildcards(config['assembly_path']+"/{sample}.fasta")

print(sample_ids)

rule all:
	input:
		expand("prokka/gffs/{sample}.gff", sample=sample_ids),
		"Roary.out"

rule prokka:
	input:
		config['assembly_path']+"/{sample}.fasta"
	output:
		"prokka/gffs/{sample}.gff"
	conda:
		"config/prokka.yaml"
	log:
		"prokka/logs/{sample}.log"
	threads:
		16
	shell:
		"""
		prokka --genus Escherichia --usegenus --cpus 0 --outdir prokka/{wildcards.sample}.out {input}
		cp prokka/{wildcards.sample}.out/*.gff {output}
		"""

rule roary:
	input:
		expand("prokka/gffs/{sample}.gff", sample=sample_ids)
	output:
		"Roary.out"
	conda:
		"config/roary.yaml"
	log:
		"roary/logs/roary.log"
	shell:
		"roary -p 16 -f Roary.out {input}"






