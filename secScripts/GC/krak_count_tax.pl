#!/usr/bin/env perl
#./krak_count_tax.pl /scratch/bork/hildebra/SimuB/simulated_metaG3SA_4/tmp/krak/krak0.out /scratch/bork/hildebra/SimuB/simulated_metaG3SA_4/tmp/krak/krak0.tax
use warnings;
use strict;

my $inF = $ARGV[0];
my $outF = $ARGV[1];

my %cnts;
open I,"<$inF" or die "can't open ifile $inF\n"; 
while (my $lin = <I>){
	chomp $lin;
	my @spl = split /\t/,$lin;
	if ($spl[1] !~ /^d_/){
		#die "$lin\n";
		next;
	}
	my @spl2 = split /\|/,$spl[1];
	#die @spl2. "gg\n";
	my $ctax = "$spl2[0]";
	for (my $i=1;$i<7;$i++){
		if ($i >= @spl2){
			$ctax .= ";?";
		} else {
			$spl2[$i] =~ s/^\S__//;
			$ctax .= ";".$spl2[$i]
		}
	}
		#die "$ctax\n$lin\n";
	if (!exists $cnts{$ctax}){
		$cnts{$ctax} = 1;
	} else {
		$cnts{$ctax} ++;
	}

}
close I;

open O,">$outF";
foreach my $k (keys %cnts){
	print O $k."\t".$cnts{$k}."\n";
}
close O;