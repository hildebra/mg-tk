#!/usr/bin/env perl
se warnings;
use strict;
use Mods::GenoMetaAss qw(readFasta );
use Mods::IO_Tamoc_progs qw( getProgPaths );

#input: contigs, tmp_path, gff for contigs, #threads
my ($ctgFile, $tmpPath, $gffF, $numCore) = @ARGV;

my $whok = getProgPaths("whokaryote");
my $cmd = "$whok --contigs $ctgFile --outdir $tmpD/whoK/ --gff $gffF --minsize 4000 --threads $cores\n";
systemW $cmd;


my $ratio = 0;
print "Bacs=$bacCnt, Euks=$EukCnt, Ratio:$ratio\n";


#Output
#1 file with bacterial contigs: $inputBac = "$GlbTmpPath/bact.kraken.fasta";
#1 file with euk contigs: "$GlbTmpPath/euk.kraken.fasta";

