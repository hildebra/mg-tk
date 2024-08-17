#!/usr/bin/perl
#./unsplit_krak.pl /scratch/bork/hildebra/GNMass3/peacemaker-29-918-0/tmp/rawRds/x1.fq  /scratch/bork/hildebra/GNMass3/peacemaker-29-918-0/tmp/rawRds/y1.fq /scratch/bork/hildebra/GNMass3/peacemaker-29-918-0/tmp/rawRds/y2.fq
use warnings;
use strict;

my $inF = $ARGV[0];
my $outF1 = $ARGV[1];
my $outF2 = $ARGV[2];

open O2,">$outF2" or die "Can;t open $outF2\b"; open O1,">$outF1" or die "Can;t open $outF1\b";
open I,"<$inF" or die "Can;t open $inF\b";

my $NU = "";my $seq =""; my $qual = "";
while (<I>){
	if (m/(^@.*);(\d+);/){
		my $len = $2 ;
		print O2 $1."/2\n"; print O1 $1."/1\n";
		$seq = <I>;$NU = <I>;$qual = <I>;
		chomp($seq);  chomp($qual);
		print O1 substr($seq,0,$len)."\n+\n".substr($qual,0,$len)."\n";
		print O2 substr($seq,$len+1,$len)."\n+\n".substr($qual,$len+1,$len)."\n";
	} else {
		die "out of fq cycle: $_\n";
	}
}



close I; close O2; close O1;
