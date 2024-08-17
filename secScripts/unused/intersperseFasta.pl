#!/usr/bin/env perl
#parseBlastFunct.pl /g/bork5/hildebra/data/Ach_proteins/res/DiaAssignment.txt
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipopen);


my $fasF1 = $ARGV[0]; my $fasF2 = $ARGV[1];

my $I1 = gzipopen $fasF1,"fasta"; my $I2 = gzipopen $fasF2,"fasta";


my $curI = $I2; my $storeLine=<$I1>; my $mode=2; 
while (my $line = <$curI>){
	#chomp $l2;
	if ($line =~ m/^>/){
		print $storeLine;
		$storeLine = $line;
		if ($mode == 1){
			$curI = $I2;$mode=2;
		} else {
			$curI = $I1;$mode=1;
		}
	} else { #not  a head, just print
		print $line;
	}
}


close I1; close I2;

