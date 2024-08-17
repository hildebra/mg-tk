#!/usr/bin/env perl
#gets based on luis motu clustering script the genomes
#see also helpers/annotateMGwMotus.pl to annotate gene cat motus
#./motu_genome_integrate.pl
use warnings;
use strict;

use Mods::GenoMetaAss qw(readClstrRev );
use Mods::IO_Tamoc_progs qw(getProgPaths);

sub readMotuTax;
sub readGene2mlinkage;

my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";

my $GCd = $ARGV[0];
my $motuName = $ARGV[1];
my $motuClus = "X"; $motuClus = $ARGV[2] if (@ARGV > 2);


#get similar genes from motu catalog

#.. and compare to gene cat


#make an abundance matrix for the motu based clustered genes (Luis), to be later used in my R filtering scripts
if ($motuClus ne "X"){ 
	my $motuGenes = `cat $motuClus`; chomp $motuGenes; my @mGenes; 
	foreach( split /\n/,$motuGenes ){my @spl=split /\t/; next if ($spl[0] eq "Genes"); push @mGenes,$spl[0];}
	my $motuGF = $oDir."/motu.genes.tmp"; my $mOTUomat = "$oDir/$tarName.mOTU.mat";
	open O,">$motuGF";print O join("\n",@mGenes);close O;
	my $cmd = "$rarBin lineExtr -i $GCd/Matrix.mat -o $oDir/mOTU.mattmp -reference $motuGF\n";
	$cmd .= "head -n 1 $GCd/Matrix.mat > $mOTUomat;cat $oDir/mOTU.mattmp >> $mOTUomat; rm $oDir/mOTU.mattmp $motuGF";
	system $cmd;
	#die $cmd;
	print "Extracted ".@mGenes." genes (motu based) from $motuClus\n";
	#die;
}

