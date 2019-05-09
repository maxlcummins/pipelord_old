import sys, csv, pysam

		fasta_file = pysam.FastaFile(sys.argv[1])
		with open(sys.argv[2], 'wb') as csvfile:
    		fasta_writer = csv.writer(csvfile, delimiter=',')
    		fasta_writer.writerow(['Ref','Length','Seq'])
    		for ref in fasta_file.references:
    			fasta_writer.writerow([ref, fasta_file.get_reference_length(ref), fasta_file.fetch(ref)])