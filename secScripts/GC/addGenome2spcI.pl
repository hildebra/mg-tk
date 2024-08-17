#!/usr/bin/env perl
#takes as input a genome (with / without genes, with / without FMGs) and adds these to an existing specI database 

use warnings;
use strict;

use Mods::GenoMetaAss qw( readFasta writeFasta);
use Mods::IO_Tamoc_progs qw(getProgPaths);

use Mods::phyloTools qw( getGenoGenes getFMG);


my $tarDir = "/g/bork5/hildebra/results/TEC2/v5/Genomes/T2/";
my $genoN = "6666666.214148.fna";
my $name = "TEC2";
my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI_f11/"; my $freeze11=0;
#my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI_2017";my $freeze11=0;

my %FMGcutoffs = (COG0012=>94.8,COG0016=>95.8,COG0018=>94.2,COG0172=>94.4,COG0215=>95.4,COG0495=>96.4,COG0525=>95.3,COG0533=>93.1,COG0541=>96.1,
COG0552=>94.5,COG0048=>98.4,COG0049=>98.7,COG0052=>97.2,COG0080=>98.6,COG0081=>98,COG0085=>97,COG0087=>99,COG0088=>99,COG0090=>98.8,COG0091=>99,
COG0092=>99,COG0093=>99,COG0094=>99,COG0096=>98.6,COG0097=>98.4,COG0098=>98.7,COG0099=>98.9,COG0100=>99,COG0102=>99,COG0103=>98.4,
COG0124=>94.5,COG0184=>98.2,COG0185=>99,COG0186=>99,COG0197=>99,COG0200=>98.4,COG0201=>97.2,COG0202=>98.4,COG0256=>99,COG0522=>98.6);

my $tarG = "$tarDir/$genoN";
my $ncore = 10;
my ($genes,$prots) = getGenoGenes($tarG);
#die "$prots\n";
my $FMGdir = getFMG("",$prots,$genes,$ncore);
#die $FMGdir."\n";
my @mvtars;
my $xtrLab= "";$xtrLab= ".rep" if ($freeze11);
my $geneNameConsistent = "";
foreach my $COG (keys %FMGcutoffs){
	my $DBnt = "$SpecID/$COG${xtrLab}.fna";
	my $DBaa = "$SpecID/$COG${xtrLab}.faa";
	die "SpecI DB nt or aa missing: $DBnt\n" if (!-e $DBnt || !-e $DBaa);
	my $hrFNA = readFasta($DBnt); my $hrFAA = readFasta($DBaa); 
	my $thn = readFasta("$FMGdir/$COG.fna");my $tha = readFasta("$FMGdir/$COG.faa");
	my $idx=0; 
	my @thK = keys %{$thn};
	
	if (@thK > 1){my $maxL=0; 
		my $cnt=0;foreach (@thK){if (length(${$thn}{$_}) > $maxL){$maxL = length(${$thn}{$_});$idx=$cnt;}$cnt++;}
	}
	my $tarKey = $thK[$idx];
	#keep specI format: 6666666.214148{_  .}COG0090 / 1410658.PRJNA223501
	my $tarKeyS = $tarKey; $tarKeyS =~ s/_COG/\.COG/;
	#die "$tarKeyS\n";
	if ($tarKeyS !~  m/^(.*\..*)\./){die "$tarKey does not follow format required by specIs\n";}
	if ($geneNameConsistent eq ""){$geneNameConsistent = $1 ;
	}elsif ($geneNameConsistent ne $1){
		die "Conflicting consisten gene names:\n$1  ..  $geneNameConsistent\n";
	}
	if (exists(${$hrFNA}{$tarKey}) || exists(${$hrFNA}{$tarKey}) ){
		die "$tarKey is already in $DBnt / $DBaa\n";
	}
	#add DNA seq to COG database
	${$hrFNA}{$tarKeyS} = ${$thn}{$tarKey};
	${$hrFAA}{$tarKeyS} = ${$tha}{$tarKey};
	#write COG database
	writeFasta($hrFNA,$DBnt.".bk");writeFasta($hrFAA,$DBaa.".bk");
	push(@mvtars,$DBnt.".bk",$DBaa.".bk");
}

#open OC,">>$SpecID/progenomes.specIv2_2";
open OC,">>$SpecID/prok-refdb-v11.0.0_specI-v2_clusters-v1.map";
print OC "$name\t$geneNameConsistent\trepresentative\n";
close OC;
open OC,">>$SpecID/specI.tax2";
print OC "$name\tBacteria\tFirmicutes\tClostridia\tClostridiales\tEriseae\tBorkfalki\tBorkfalki ceftriaxensis\n";
close OC;

foreach my $mvv (@mvtars){
	my $mvt = $mvv;$mvt =~ s/\.bk$//;
	system "rm $mvt;mv $mvv $mvt\n";
}
print "added $name to specI database $SpecID\n";
exit (0);