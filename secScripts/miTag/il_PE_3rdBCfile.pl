#!/usr/bin/perl

use warnings;
use strict;

my $ifile = "/g/bork5/hildebra/transfer/A88BF.1.fq";
my $ofile = "/g/bork5/hildebra/transfer/A88BF.mid.fq";
open I,"<",$ifile;
open O,">",$ofile;
my $cnt =-1;
while (my $line = <I>){
	$cnt++;
	if ($cnt%4 !=0){
		next;
	}
	chomp $line;
	#print $line."\n";
	die ("Wrong fastq header $cnt:\n$line\n") unless ($line=~m/^@/);
	$line =~ m/.*:([ACGTN]+)/;
	my $BC = $1;
	#die $BC."\n";
	print O $line."\n".$BC."\n"."+\n".( "A" x length($BC) )."\n";
	
	#die ($ofile) if ($cnt > 100);
}


close I; close O;