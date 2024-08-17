#!/usr/bin/perl
#uses MGS to cluster genes and kraken tax assignments
#./taxPerMGS.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//MB2.clusters.ext.can.Rhcl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/
use warnings; use strict;

use Mods::Binning qw(readMGSrev );


my $refMGf = $ARGV[0];
my $GCd = $ARGV[1];
my $outF = $ARGV[2];



my $LCAout="$outF.LCA";
my $taxout = "$outF.tax";

if (-e $LCAout && -e $taxout){exit(0);}

my $hr = readMGSrev($refMGf);
my %MGs = %{$hr};
my %tCnt;
my $krakF = "$GCd/Anno/Tax/krak2.txt"; #krak_0.01.txt
#my $krakF = "$GCd/Anno/Tax/krak_0.01.txt"; #
open I,"<$krakF" or die "Can't open kraken input $krakF\n";
while (<I>){
	chomp;my @spl=split /\t/;
	next unless (exists($MGs{$spl[0]}));
	#my @s2 = split /;/,$spl[1];
	unless (defined $spl[1]){
		next;
	}
	my $i=0;
	foreach my $t (split /;/,$spl[1]){
		$tCnt{$MGs{$spl[0]}}{$i}{$t}++;
		$i++;
	}
}
close I;
my %tStat;
open OL,">$LCAout";
open OC,">$taxout";
foreach my $mgs (sort keys %tCnt){
	my $tax=""; my $maxD=0;
	print OC "$mgs\t";
	print OL "$mgs\t";
	for (my $i=0;$i<8;$i++){
		next unless (exists( $tCnt{$mgs}{$i} ));
		my %curT = %{$tCnt{$mgs}{$i}}; my $tSum=0;
		my $max=0;my$maxT="";
		print OC "\t$i";
		foreach my $t (keys %curT){
			$tSum+=$curT{$t};
			print OC "$t:$curT{$t};";
			if ($curT{$t} > $max){
				$max = $curT{$t};
				$maxT = $t;
			}
		}
		next if ($tSum < 100);
		if (($max/$tSum) > 0.8){
			$tax .= "${maxT};";
			$maxD=$i;
		} else {
			$tax .= "?;";
		}
	}
	$tStat{$maxD}++;
	print OC "\n";
	print OL "$tax\n";
}
#print stats how many LCA levels were hit by MGS:
foreach my $d (sort(keys%tStat)){
	print "$d:$tStat{$d}\t";
}
print "\n";


