import sys, csv, pysam
#https://www.quora.com/How-does-one-convert-fasta-to-CSV-using-python
fasta_file = pysam.FastaFile(sys.argv[1])
with open(sys.argv[2], 'wb') as csvfile:
    fasta_writer = csv.writer(csvfile, delimiter=',')
    for ref in fasta_file.references:
        fasta_writer.writerow([ref,"http://api.bigsi.io/search?seq="+fasta_file.fetch(ref)+"&threshold=1"])
