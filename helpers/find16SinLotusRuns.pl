#!/usr/bin/perl
#takes predicted 16S and scans real 16S experiments for given 16S subset

use warnings;
use strict;

my $mkBldbBin = "/g/bork5/hildebra/dev/lotus//bin//ncbi-blast-2.2.29+/bin/makeblastdb";
my $blastBin = "/g/bork5/hildebra/dev/lotus//bin//ncbi-blast-2.2.29+/bin/blastn";

my $lotusD = "/g/bork3/home/hildebra/results/lotus/HMP_35_sw2//";#HMP_35_SLV//";
my $otuTar = $lotusD."otus.fa";
my $tar16sDB = "/g/bork5/hildebra/results/TEC2/v5/Genomes/T2/6666666.214148.fna.16S";
my $taxblastf = "/g/bork5/hildebra/results/TEC2/v5/Genomes/T2/HMPhits.blast";
#my $TECdir = "/g/scb/bork/hildebra/SNP/GNMass2_singl/alien-11-374-0/Binning/";
#my $tar16sDB = $TECdir."TEC16.fa";
#my $taxblastf = $TECdir."HMPHits.blast";

my $cmd = "$mkBldbBin -in $tar16sDB -dbtype 'nucl'\n";
unless (-f $tar16sDB.".nhr"){	system($cmd);}
my $strand = "both";
#-perc_identity 75
$cmd = "$blastBin -query $otuTar -db $tar16sDB -out $taxblastf -perc_identity 95 -outfmt 6 -max_target_seqs 50 -evalue 0.001 -num_threads 60 -strand $strand \n"; #-strand plus both minus
system $cmd."\n";