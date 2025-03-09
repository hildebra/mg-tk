#!/usr/bin/env python
## program to reconstruct fasta sequence
import argparse
import sys
import re
import gzip
import numpy as np
from Bio import SeqIO
#from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC


## Returns nt and AA gene files from a contig file and gene annotation file (gff)
## Run like this:
## python contigs2fasta.py --gff [gff file] --contig [new SNP contig file] --outFNA [nt gene file] --outFAA [AA gene file]
## If input or output files ends with ".gz", the program assummes input or output should read/written as compressed

def main():

	if sys.version_info[0] < 3:
		sys.stderr.write("Require Python 3 or newer. You are currently using:\n" + str(sys.version_info) + "\n")
		sys.exit(1)

	parser = argparse.ArgumentParser(description='Returns nt and AA gene files from a contig file and gene annotation file')
	parser.add_argument("-g", "--gff", dest="gff_file", help="gff file",required=True)
	parser.add_argument("-c", "--contig", dest="contig_file", help="new SNP contig file", required=True)
	parser.add_argument("-n", "--outFNA", dest="fna_file", help="output nt gene file", required=True)
	parser.add_argument("-a", "--outFAA", dest="faa_file", help="output AA gene file", required=True)

	options = parser.parse_args()

	## Initialise
	gff_file = options.gff_file
	contig_file = options.contig_file
	fna_file = options.fna_file
	faa_file = options.faa_file

	print("gff file:", gff_file)
	print("contig file:", contig_file)
	print("output FNA:", fna_file)
	print("output FAA:", faa_file)

	transl_table = None
	reTransl_table = re.compile(r"^# Model Data:.+transl_table=(\d+);")
	reID = re.compile(r"^ID=\d+_(\d+);.+start_type=(.+);")
	posfreq = re.compile(r"POS=(.*) FR=(.*) FREQT")
	refreq = re.compile(r" COV=(\d+) .* FREQT=(.*) CONFL=")


	## Save contigs as dict
	if contig_file[-3:] == ".gz":
		handle = gzip.open(contig_file, 'rt')
	else:
		handle = open(contif_file, 'r')
	if handle.closed:
		print("Could not open handle:",handle)
		exit(5)

	contig_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	genekeys = contig_dict.keys()
	handle.close()


	## Open outfiles
	if fna_file[-3:] == ".gz":
		outFNA = gzip.open(fna_file, 'wt')
	else:
		outFNA = open(fna_file, 'w')
	if outFNA.closed:
		print("Could not open outFNA:",outFNA)
		exit(5)

	if faa_file[-3:] == ".gz":
		outFAA = gzip.open(faa_file, 'wt')
	else:
		outFAA = open(faa_file, 'w')
	if outFAA.closed:
		print("Could not open outFAA:",outFAA)
		exit(5)
	
	## Open gff file

	if gff_file[-3:] == ".gz":
		fp = gzip.open(gff_file, 'rt')
	else:
		fp = open(gff_file, 'r')
	if fp.closed:
		print("Could not open fp:",fp)
		exit(5)
		
	## calc preliminary chance of low freq second strain
	contig_conta = {}
	for rec in contig_dict:
		prob=0
		contig_conta[rec] = prob
		freReg = refreq.search(contig_dict[rec].description)
		freqs = np.array(freReg.group(2).split(","))
		freqs[freqs=='']='0'
		freqs=freqs.astype(float)
		sum = np.sum(freqs)
		if sum < 5:
			continue
		#Slen = len(contig_dict[rec].seq)
		Slen = int(freReg.group(1))
		midSum=np.sum(freqs[20:80])
		prob = midSum/Slen*100
		#print(sum,midSum,Slen,prob)
		contig_conta[rec] = prob


	## Pass through gff file
	for i, line in enumerate(fp):

		if line[0] == "#": ## Comment line
			transl_table=None
			searchTransl_table = reTransl_table.search(line)
			if searchTransl_table: ## Set IUPAC translation table for AA
				transl_table=int(searchTransl_table.group(1))

		else: ## Gene -> FNA and FAA

			if transl_table is None: ## Exit if didn't find translation table
				sys.stderr.write("Error! Did not find translation table for gene, line " + str(i) + "\n")
				sys.exit(1)

			line = line.split()

			if line[0] in genekeys:

				start = int(line[3])
				end = int(line[4])
				strand = line[6] == "-"
				geneid = reID.search(line[8])
				gene = line[0] + "_" + geneid.group(1)
				geneXtra = " P=[] F=[] oCSP=" + str(contig_conta[line[0]]) +" CSP=0"
				
				#print(contig_dict[line[0]].description )
				posReg = posfreq.search(contig_dict[line[0]].description )
				dna = contig_dict[line[0]].seq[start-1:end]
				DNAchars = dna.count("A") + dna.count("T") + dna.count("C") + dna.count("G")
				if (DNAchars == 0):
					#print(dna)
					continue
				if posReg.group(1) != "":
					posSNP = np.array(posReg.group(1).split(","),dtype=int)
					posF = np.array((posReg.group(2).split(",")),dtype=float)
					uses = (posSNP >= start) & (posSNP <= end)
					if any(uses):
						#print(posSNP, start, end)
						#print(uses)
						posF = posF[uses]
						posSNP = posSNP[uses]
						if strand: #turn pos around
							geneXtra = " P=" + np.array2string((end-start) - (posSNP-start),separator = ",") 
						else:
							geneXtra = " P=" + np.array2string(posSNP-start,separator = ",") 
						geneXtra +=  " F=" + np.array2string(posF,separator = ",",precision =3)
						#determine prob of conspecific strains
						CSP = (np.sum(posF<0.8)  * contig_conta[line[0]]) *100 /(DNAchars)
						geneXtra += " oCSP="+ str(contig_conta[line[0]])+" CSP=" + np.array2string(CSP,precision=4)
						#print(gene+geneXtra)
						
						#sys.exit()
				
				start_type = geneid.group(2)

				## Reverse complement for minus strand genes
				if strand: ## for the minus strand, get reverse complement
					dna = dna.reverse_complement()
				
				gene = gene+geneXtra
				## FNA
				outFNA.write(SeqRecord(dna, id=gene, description="", name="").format("fasta"))
				
				## FAA (Tranlate)
				if start_type == "Edge": # Do not use alternative codons
					stringFAA = SeqRecord(dna.translate(table=transl_table), id=gene, description="", name="").format("fasta")
				else:
					try: # Use alternative codons and add stop codon, which is removed
						stringFAA = SeqRecord(dna.translate(table=transl_table, cds=True) + Seq("*"), id=gene, description="", name="").format("fasta")
					except:
						try: # Test to see if translation failed do to 3' gene
							dna = dna + Seq("TAA", IUPAC.unambiguous_dna)
							stringFAA = SeqRecord(dna.translate(table=transl_table, cds=True), id=gene, description="", name="").format("fasta")
						except:
							stringFAA = SeqRecord(dna[:-3].translate(table=transl_table), id=gene, description="", name="").format("fasta")
			
				outFAA.write(stringFAA)


	## Close file handles
	fp.close()
	outFNA.close()
	outFAA.close()
    
	print ("Done rewriting\n")

if __name__ == "__main__":
	main()

