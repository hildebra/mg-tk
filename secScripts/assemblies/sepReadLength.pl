#!/usr/bin/perl
#takes a list of FNA's and prints every fna longer than X bp to seperate file
# ./sepReadLength.pl 7000 /tmp/hildebra/GC/35Pcompl.fna

use strict; use warnings;

if (@ARGV != 2){die "requires 2 args!!!\n";}
my ($len,$inF) = @ARGV;

open I,"<$inF";
open OL,">$inF.long";
open OS,">$inF.short";
my $fna =""; my $bps =0;
while (<I>){
	if ( $_ !~ m/^>/){
		$bps += length($_)
	} else {
		#decide where to write
		if ($bps > $len){
			print OL $fna;
		} else {
			print OS $fna;
		}
		$fna = ""; $bps =0;
	}
	$fna .= $_;
}
if ($bps > $len){
	print OL $fna;
} else {
	print OS $fna;
}

close OS; close OL; close I;

my $iniFNA = `wc -l $inF | cut -f1 -d ' ' `;
my $lFNA = `wc -l $inF.long | cut -f1 -d ' ' `;
my $sFNA = `wc -l $inF.short | cut -f1 -d ' '`;
print "$iniFNA $lFNA $sFNA\n";

#system etc
system "rm $inF;mv $inF.short $inF";