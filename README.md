# Pipelord
This repository contains a number of Snakemake pipelines that can be used for common bacterial genomic analyses workflows performed in the Djordjevic Genomics lab.

If you need to learn about Snakemake installation and usage please check the documentation [here](https://snakemake.readthedocs.io/en/stable/).

## Troubleshooting and FAQ
If you have any difficulty with these workflows please click the issues tab on this repo and file an issue report with as much detail as you think might be necessary for us to help troubleshoot things. If you think you know a fix for the issue then please include this. Please make sure that nobody has filed an issue regarding the same issue that you have encountered prior to posting a new one.

## Will it -j?
If you know, you know. Currentlt, builds of Snakemake pipelines that are '-j'able include:
| Name  | -j status |
| ------------- | ------------- |
|  abricatelord  ||
|  aribalord  ||
|  BISGIlord  ||
|  BLASTlord  ||
|  pangenomelord  ||
|  phylosiftlord  ||
|  plasmidlord  ||
|  shuvlord  |&#9745;|
|  snplord  ||
|  SRAlord  ||




## abricate and abricatelord
[Abricate](https://github.com/tseemann/abricate) is a tool developed by Torsten Seeman for the genotyping of genomic datasets.

The tool is unpublished so there is no paper to cite, but if you use it for any manuscripts make sure you link the Github repo in your methodology. Also ensure you cite/provide access to any databases you use.

You can use abricate to screen for genes of interest from commonly used publically available databases (see the abricate readme for supported databases), or you can create custom databases. Note that it screens assemblies (i.e. fasta files) and not read sets (i.e. fastq files). If you want to screen read sets (such as for screening of genes that typically are less well detected using assembly based screening tools like IS elements) use ARIBA instead.

abricatelord can be used to streamline and parallelise the genotyping of collections of samples with multiple databases at the same time. It outputs individual abricate reports which currently require manual concatenation (using cat) or summary using (abricate --summary; see abricate wiki for details), though we will soon impliment a summary step in the Snakefile.

## aribalord
Antimicrobial Resistance Identification By Assembly-[ARIBA](https://github.com/sanger-pathogens/ariba) is a tool developed by Martin Hunt et al from the Sanger Institute which can be used for the genotyping of genomic datasets.

This tool is published - cite it appropriately. Also ensure you cite/provide access to any databases you use.

Similarly to abricate, you can use ARIBA to screen for the presence of genes of interest. ARIBA also integrates certain public databases and allows the use of custom databases. It is more computationally intensive than abricate is, so will take longer to run, but it can be performed on read sets rather than genomic assemblies and therefore does not necessitate genomic assembly prior to genotyping. It is also useful for identification of IS elements which don't tend to assemble very well (in our experience IS carriage is higher as determined by ARIBA in comparison with abricate and other BLAST based genotyping tools).

aribalord can be used to streamline and parallelise the genotyping of collections of samples with multiple databases at the same time. It outputs multiple ariba summary reports which can be combined using a python script inconveniently also . 

## BIGSI and BIGSIlord
BItsliced Genomic Signature Index [BIGSI](http://www.bigsi.io/) is a tool developed by Phelim Bradley et al for the screening of genomic datasets from the Sequence Read Archive [SRA](https://www.ncbi.nlm.nih.gov/sra) for query sequences of interest.

This tool is published - cite it appropriately. Also ensure you cite/provide access to any reference sequences you use.

BIGSIlord can be used to streamline and parallelise the searching of massive online databases for sequences of interest. Multiple sequence can be queried at once, and the output of BIGSIlord consists of a text file containing a list of accession numbers and their associated species identifications. This text file can then be fed into SRAlord to allow the downloading of the associated genomic datasets.

## BLASTn and BLASTlord
BLASTn is one of the most widely used sequence alignment tools available and has been around for decades. It has many uses but primarily we use it for genotyping of genomic datasets. There are command line versions of BLAST, as well as GUI based ones. BLASTn is for nucleotide sequence alignments, while BLASTp is used for amino acid sequence alignment. There are other types of BLAST also. You can find the online BLAST portal [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi), while the BLAST command line can be found (here)[https://www.ncbi.nlm.nih.gov/books/NBK279690/].

This tool is published - cite it appropriately. Some believe BLAST doesnt need to be cited because it is so well known - I feel it should be cited. Also ensure you cite/provide access to any reference sequences you use.

BLASTlord within this repo refers to a Snakefile that allows you to run BLAST on lots of samples with lots of different databases in a lazy fashion. It generates a text file for each of the reference databases used which are named according to the names of the databases.

**Note: Confusingly, maxlcummins has written two different tools with the same name. There is also an R package called [BLASTlord](https://github.com/maxlcummins/BLASTlord) which predates this repo by several years. Names will be changed to avoid confusion - I am thinking BLASTlordR for the R package. Apologies for this.**


## Authorship
[Max Cummins](https://github.com/maxlcummins/), [Cameron Reid](https://github.com/cjreid)

## Acknowledgements
[Dmitriy Li](https://github.com/Tu6ka)
