#!/usr/bin/perl
#script to get TEC2, TEC6 gene names for sample A374, A377 to plot on graph for paper

use strict; use warnings;
use Mods::GenoMetaAss qw(readClstrRev);


my $oname = "T6";
my $odir = "/g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/$oname/R_filt/";
my $refF = "$odir/${oname}_filt.txt";
my $GCdir = "/g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/";
my %genes ;
open I,"<$refF" or die "can;t open $refF\n"; while(<I>){chomp; next if ($_ eq "");$genes{$_}=1;} close I;

my @ks = sort(keys %genes);
#die "$ks[0]X\n";

#my ($geneIdxH,$numGenes) = readGeneIdx($GCdir."Matrix.genes2rows.txt");
#my %geneIdx = %{$geneIdxH};

#my ($hr1,$hr2) = readClstrRev("$GCdir/compl.incompl.95.fna.clstr.idx");
#my %gene2cl = %{$hr1}; my %cl2gene = %{$hr2};

my $cnt =0;
my %rep;
open I,"<$GCdir/compl.incompl.95.fna.clstr.idx" or die "can t open GC idx\n";
while (<I>){
	$cnt++;
	chomp; my @spl = split /\t/; next if (@spl < 2);
	#print $spl[0]."\n";
	
	if (exists($genes{$spl[0]})){
		my @spl2 = split /,>/,$spl[1];
		$spl2[0] =~ s/^>//;
		#die "@spl2\n";
		foreach my $g (@spl2){
			$g =~ m/^(.*)__/;
			$rep{$1}{$g} = 1;
		}
	}
}
close I;

my @smpls = sort keys %rep;
open O,">$odir/$oname.genes.per.smpl.txt" or die "can't open $odir/$oname.genes.per.smpl.txt\n";
foreach my $s (@smpls){
	my @ke = keys %{$rep{$s}};
	print O "$s\t".join(',',@ke)."\n";
}
close O;
print "$odir/$oname.genes.per.smpl.txt\n";

