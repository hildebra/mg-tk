#!/usr/bin/perl
use strict;
use warnings;

my $inMMS = $ARGV[0];
my $outCLSTR = $ARGV[1];

die "not enough input args\n" if (@ARGV < 2);

open I,"<$inMMS" or die "Can't open input $inMMS\n";
open O,">$outCLSTR" or die "Can't open output $outCLSTR\n";

my $cnt =1 ; my $wiCnt=0; my $lcore = "";
my $totMem=0;
while (my $line = <I>){
	chomp $line;
	my @spl = split /\s+/,$line;
	my $ccore=  shift @spl;
	if ($lcore ne $ccore){
		print O ">Cluster $cnt\n"; $cnt ++;
		$lcore = $ccore;
		$totMem += $wiCnt;
		$wiCnt = 0;
	}
	my $hit = " at +/99.1%";
	$hit = " *" if ($ccore eq $spl[0]);
	print O "$wiCnt\t200nt, >$spl[0]...$hit\n";
	$wiCnt++;
}

print "\n\nFound $cnt clusters, $totMem members, ". $totMem/$cnt ." average size\n\n";

close I; close O;