# Pipelord
This repository contains a number of Snakemake pipelines that can be used for common bacterial genomic analyses workflows performed in the Djordjevic Genomics lab.

## abricatelord
[Abricate](https://github.com/tseemann/abricate) is a tool developed by Torsten Seeman for the genotyping of genomic datasets.

You can use abricate to screen for genes of interest from commonly used publically available databases (see the abricate readme for supported databases), or you can create custom databases. Note that it screens assemblies (i.e. fasta files) and not read sets (i.e. fastq files). If you want to screen read sets (such as for screening of genes that typically are less well detected using assembly based screening tools like IS elements) use ARIBA instead.

