#!/usr/bin/perl
#takes a list of FNA's and prints every fna longer than X bp to seperate file
# ./sepReadLength.pl 7000 /tmp/hildebra/GC/35Pcompl.fna

use strict; use warnings;

sub splitFNA($ $ $ $){
	my ($HD,$fna,$len,$OL) = @_;
	my $pos=0; my $tlen = length($fna);
	while($pos < $tlen){
		my $subFNA= substr($fna,$pos,$len);
		$pos += $len;
		print $OL $HD."\n".$subFNA."\n";
	}
}



if (@ARGV != 2){die "requires 2 args!!!\n";}
my ($inF,$len) = @ARGV;

open I,"<$inF";
open my $OL,">$inF.tmp";
my $fna ="";my $HD; my $bps =0; my $cnt =0; my $isFQ=0;
my $splCnt=0;
while (<I>){
	if ($cnt==0){#detect fastq
		if ($_ =~ m/^@/){$isFQ=1; print "Fastq file, cannot process\n";last;}
	}
	chomp;
	if ( $_ !~ m/^>/){
		$fna .= $_;
	} else {
		$bps = length($fna);
		#decide where to write
		if ($bps> 0){
			if ($bps < $len){
				$splCnt++;
				print $OL $HD."\n".$fna."\n";
			} else {
				splitFNA($HD,$fna,$len,$OL);;
			}
		}
		$fna = ""; 
		$HD = $_;
	}
	$cnt++;
}
if ($isFQ){
	close $OL; system "rm $inF.tmp"; exit(0);
}


$bps = length($fna);
if ($bps < $len){
	print $OL $HD."\n".$fna."\n";
} else {
	splitFNA($HD,$fna,$len,$OL);;
}

close $OL; close I;
print "Found $splCnt sequences longer than $len\n"; 
system "rm $inF;mv $inF.tmp $inF";