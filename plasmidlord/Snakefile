import re
import os

configfile:
	"config/config.yaml"

sample_ids, = glob_wildcards(config['read_path']+"/{sample}.R1.fastq.gz")
ref, = glob_wildcards(config['ref_path']+"/{reference}.fa")

print(sample_ids)
print(ref)

rule all:
	input:
		expand("ref/{reference}.fa.amb", reference = ref),
		expand("sorted_reads/{sample}.{reference}.bam", sample=sample_ids, reference=ref),
		expand("Coverage_files/{sample}.{reference}.bam_coverage.txt", sample=sample_ids, reference=ref),
		expand("heatmap_hclust.pdf", reference=ref)


rule index_ref:
	input:
		"ref/{reference}.fa"
	output:
		"ref/{reference}.fa.amb"
	priority: 10
	conda:
		"config/plasmidlord.yaml"
	shell:
		"bwa index {input}"
		
rule bwa_map:
   	input:
        	ref = "ref/{reference}.fa",
        	r1 = config['read_path']+"/{sample}.R1.fastq.gz",
		r2 = config['read_path']+"/{sample}.R2.fastq.gz"
	output:
		"sorted_reads/{sample}.{reference}.bam"
	priority:9
	log: 
		"logs/bwa_mem/{sample}.{reference}.log"
	conda: 
		"config/plasmidlord.yaml"
	shell:
        	"bwa mem {input.ref} {input.r1} {input.r2} | samtools view -ubS -F 0x904 - | samtools sort -o {output}"


#checkpoint samtools_sort:
#	input:
#		"mapped_reads/{sample}.{reference}.bam"
#	output:
#		"sorted_reads/{sample}.{reference}.bam"
#	priority:8
#	log: 
#		"logs/samtools/{sample}.{reference}.log"
#	conda: 
#		"config/plasmidlord.yaml"
#	shell:
#		"samtools sort -O bam {input} -o {output}"

rule samtools_depth:
	input:
		"sorted_reads/{sample}.{reference}.bam"
	output:
		"Coverage_files/{sample}.{reference}.bam_coverage.txt"
	priority:7
	conda: 
		"config/plasmidlord.yaml"
	shell:
		"samtools depth {input} > {output}"


rule heatmap:
	output:
		"heatmap_hclust.pdf"
	params:
		bins = config['bin_size'],
		tick = config['tick_size']
	conda:
		"config/plasmidlord.yaml"
	priority: 0
	shell:
		"python2 plasmidlord.py -b {params.bins} -t {params.tick} Coverage_files heatmap"





