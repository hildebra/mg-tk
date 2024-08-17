This is a version that creates a non-redundant multi-kingdom gene catalog and runs on a GNU/Linux parallel processing system.

DISTRIBUTION:

This file				README.txt
Main shell script		genecatalog.sh
Helper script			bin/fastaNamesSizes.pl
						filtering_fromassemblies.pl
						filtering_fromformatching.pl
						finalgeneparsing.pl
						getAnnoFast.pl
						getgenes.pl
						match_nodes_bac.pl
						match_nodes_euk.pl
						overlap_nodes.pl
						process_nodes.pl


USAGE:

1. Download and install Bowtie2, Python2.7, khmer, R, seqtk, USEARCH7.0 or greater
2. Follow workflow provided in genecatalog.sh

Final files are a fasta file containing a non-redundant gene database and associated lengths file:
allgenecalled.centroids.nolines.fa
allgenecalled.lengths.txt

Or genes called only on contigs greater than 1000bp:
allgenecalled.1000.centroids.fa
allgenecalled.1000.lengths.txt


NOTES:
Gene prediction models differ not only from species to species, but kingdom to kingdom. Because of prior knowledge that Malassezia comprise the vast majority of fungi in the skin, we incorporated a phylogenetically near neighbor for our fungal gene calling model. This approach, therefore, will be sample-specific. We are also in the process of integrating a similar approach to viral gene calling, in which we target common members of the skin virome, such as phages endemic to common skin bacteria. 

The basic approach is for each contig, to call both bacterial and fungal genes. Where there are genes called from both kingdoms, gene assignment is made by BLAST against the nr database. Where there was no resolution in the kingdom, genes were marked as ambiguous or assigned to whichever caller generated a prediction. 

We then created a non-redundant catalog by clustering at a sequence identity cutoff of 0.95 and a minimum coverage cutoff of 0.9 for the shorter sequence. However, these are very stringent thresholds and it is possible that multiple fragments may in fact represent a single gene unit. 
