#!/usr/bin/env perl
use warnings;
use strict;
#./maxBin_rewrk.pl /g/bork1/hildebra/SNP/GNMassSimu/simulated_metaG2SI_3/Binning/MaxBin MB

my $inD = $ARGV[0];
my $tag = $ARGV[1];
my $num = 0;

my $inS = "$inD/$tag.summary"; 
my @bins = ("$tag.tooshort","$tag.noclass");
open I,"<$inS" or die "Can't open $inS\n"; my $cnt =0;
while (my $line = <I>){
	$cnt ++;
	next if ($cnt == 1);
	my @ss = split (/\t/,$line);
	push(@bins,$ss[0]);
}
close I;


my $binTxt = ""; my $tooshort="";  my $noclass="";
foreach my $inFas (@bins){
	print $inFas."\n";
	my $inF = "$inD/$inFas";
	open I,"<$inF" or die "Can't open $inF\n";
	while (my $line = <I>){
		if ($line =~ m/>(.*)$/){
			if ($num > 1 ){
				$binTxt .= $1."\t".$num."\n";
			} elsif($num==1) {
				$noclass .= $1."\n";
			} else {
				$tooshort .= $1."\n";
			}
		}
	}
	close I;
	$num++;
}

open O,">$inD/$tag.ctg.bin.txt"or die "Can't open output $inD/$tag.ctg.bin.txt\n";
print O $binTxt;
close O;

open O,">$inD/$tag.ctg.noclass"or die "Can't open output $inD/$tag.ctg.noclass\n";
print O $noclass;
close O;

open O,">$inD/$tag.ctg.tooshort"or die "Can't open output $inD/$tag.ctg.tooshort\n";
print O $tooshort;
close O;

foreach my $inFas (@bins){
	system "rm -f $inD/$inFas";
}

