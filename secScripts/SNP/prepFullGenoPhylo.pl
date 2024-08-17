#!/usr/bin/env perl
#./prepFullGenoPhylo.pl /g/scb/bork/hildebra/SNP/GNMass3/GlbMap/Eubacterium_ramulus.genomes
#builds phylo from genome consensus sequence, created with 2nd map
use warnings;
use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(  systemW readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
my $bts = getProgPaths("buildTree_scr");

my $inD = $ARGV[0];

opendir my $dir, "$inD" or die "Cannot open directory: $!";
my @files = sort grep { /\.fna\.gz$/ && -f "$inD/$_" } readdir($dir);
rewinddir($dir);
my @ref = grep { /\.fa$/ && -f "$inD/$_" } readdir($dir);
closedir $dir;

#die "@files\n";
die "Too many refs: @ref\n" if (@ref > 1);
print "Using ". @files . " SNP files and ". @ref ." reference genomes, already aligned\n";

my $outD = "$inD/phylo/";
system "mkdir -p $outD" unless (-d $outD);

my $outFile = "$outD/MSA.fna";
open O, ">$outFile" or die "Can't open out file $outFile\n";
my $hr = readFasta("$inD/$ref[0]",1);my %FNA = %{$hr};
my $cnt=0; my $seq = "";
my @refCtgs = sort keys %FNA;
foreach my $k (@refCtgs){
	#print ".${k}. ".length($FNA{$k})."\n";
	$seq .= $FNA{$k};
}
print O ">ref\n$seq\n";
my $tLength = length($seq);
foreach my $sF (@files){
	$sF =~ m/_([^_]+)-0/;
	my $ID = $1;
	$seq = "";
	$hr = readFasta("$inD/$sF",1); %FNA=%{$hr};
	foreach my $k (@refCtgs){
#	foreach my $k (sort keys %FNA){
#		print ".${k}. ".length($FNA{$k})."\n";
		$seq .= $FNA{$k};
	}
	print O ">$ID\n$seq\n";
	die "newSeq length (".length($seq).") != expected length ($tLength) of fna: $ID $sF)" if ($tLength != length($seq));

}
close O;

my $numCores= 20;

my $cmd = "$bts -fna $outFile  -outD $outD -isAligned 1 -runIQtree 0 -runFastTree 1 -cores $numCores  -AAtree 0 -bootstrap 000 -NTfiltCount 300 -NTfilt 0.1 -MSAprogram 2 -AutoModel 1 -iqFast 0 \n";
die "$cmd\n";