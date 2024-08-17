#!/usr/bin/env perl
#deprecated
#gets all genes associated to an MGS and abundance matrix thereof
#can be used later for R postclustering of possible MGS
#./MGSonGeneFetch.pl /g/scb/bork/hildebra/SNP/GCs/T2_GNM3_ABR /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/TEC2Names.txt
#./MGSonGeneFetch.pl /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/specific_genes/Tec2.ABR.txt
use warnings;
use strict;

use Mods::GenoMetaAss qw(readClstrRev );
use Mods::IO_Tamoc_progs qw(getProgPaths);

sub readCanopyGenes;

my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";

my $GCd = $ARGV[0];
my $geneSet = $ARGV[1];
my $refSample = "X"; $refSample = $ARGV[2] if (@ARGV > 2);

#die "$motuClus\n";

my $tarName = "Unknown";
my $oDir = "";#$1$2"."MGS/";

if ($geneSet =~ m/(.*\/)([^\/]+)Names\.txt$/){
	$tarName = $2;
	$oDir = "$1$2"."MGS/";
}else{
	$geneSet =~ m/(.*\/)[^\/]+$/;
	$oDir = $1;
}
print "output dir= $oDir\n";
system "mkdir $oDir" unless (-d $oDir);
open LOG , ">$oDir/MGSonGeneFetch.log";
$GCd.="/" unless ($GCd =~ m/\/$/);
#dir with gene catalog
#my $GCd = "$bigF/GeneCatalog/";



my ($hr1,$hr2) = readCanopyGenes($GCd."Canopy2/clusters.txt");

my %Cans = %{$hr1}; my %revCans =  %{$hr2};

my @genes;
open I,"<$geneSet" or die "can't open gene set $geneSet\n"; while (<I>){chomp; my @spl=split /\s/;push(@genes, $spl[0]);} close I;
if (@genes ==0 ){die "No genes in input geneset!\n";}

my %asCans; my $gcnt = 0; my $gmiss=0; my @glist;
foreach my $ge (@genes){
	unless (exists ( $revCans{$ge} )){
		#print "Can't find gene $ge\n" ; 
		$gmiss++;
		next;
	}
	push(@glist,$ge);
	$gcnt++;
	if (exists($asCans{ $revCans{$ge} } )){  $asCans{ $revCans{$ge} } ++;
	} else {$asCans{ $revCans{$ge} } = 1;}
}
my @kks = keys %asCans;
print "Found " . @kks . " Canopies associated using $gcnt/".($gcnt+$gmiss)." genes\n:@kks\nGenes:@glist\n";

if (@kks ==0 ){die "NO Canopy detected\n";}

#get genes 
for (my $k=0;$k<@kks;$k++){
	my $tarCan = $kks[$k];
	#die $Cans{$tarCan}."\n";
	open O ,">$oDir/$tarCan.can.genes";
	print O join("\n",split(/,/,$Cans{$tarCan}) );
	close O;
	my $cmd = "$rarBin lineExtr -i $GCd/Matrix.mat -o $oDir/$tarCan.can.mattmp -reference $oDir/$tarCan.can.genes\n";
	$cmd .= "head -n 1 $GCd/Matrix.mat > $oDir/$tarCan.can.mat;cat $oDir/$tarCan.can.mattmp >> $oDir/$tarCan.can.mat; rm $oDir/$tarCan.can.mattmp";
	system $cmd;
#get the assemblies of these genes
	if ($refSample ne "X"){
		my $cmd = "./geneListSameAssembly.pl print $oDir/$tarCan.can.genes $GCd $oDir";
		system "$cmd 2>> $oDir/MGSonGeneFetch.log";
	}
}
#lineExtr /g/scb/bork/hildebra/SNP/GCs/T2_GNM3_ABR/Matrix.mat /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2//TEC2.MaxBin.genes.matX /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2//TEC2.MaxBin.gene.list
#head -n 1 /g/scb/bork/hildebra/SNP/GCs/T2_GNM3_ABR/Matrix.mat >  /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2//TEC2.MaxBin.genes.mat;cat  /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2//TEC2.MaxBin.genes.matX >>  /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2//TEC2.MaxBin.genes.mat; rm  /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2//TEC2.MaxBin.genes.matX;

close LOG;



print "Finished\n"; exit(0);

sub readCanopyGenes{
	my ($canF) = @_;
	open I,"<",$canF or die "can't open canopy clusteres file $canF\n";
	my %Cans; my %revCans; print "Reading Canopys..";
	while (<I>){
		chomp; my @spl = split /\t/;
		if (exists ($Cans{$spl[0]})){
			$Cans{$spl[0]} .= ",$spl[1]";
		} else {$Cans{$spl[0]} = $spl[1];}
		$revCans{$spl[1]} = $spl[0];
	}
	close I;
	print "Done\n";
	return (\%Cans,\%revCans);
}

