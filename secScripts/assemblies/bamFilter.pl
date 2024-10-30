#!/usr/bin/perl
use warnings;
use strict;


my $id = 0.05; my $lengthC = 0.8; my $MAQ=20;
$id = $ARGV[0] if (@ARGV >= 1);
$lengthC = $ARGV[1] if (@ARGV >= 2);
$MAQ = $ARGV[2] if (@ARGV >= 3);
my $idFail=0; my $mapscoreFail=0; my $coverFail=0; my $totalFail=0;
my $finalCnt=0; #really the number of entries written..
my $failMask = 1 << 2;
my $bwt2sam=1; #switch to mini2, if not matching found..
#die "$failMask\n";
my $totalCnt=0; my $totalNotMapped=0;
my $prevLine="";my $line ="";
while ($line = <STDIN>){
	if ($line =~ m/^@/){print $line; next;}
	chomp $line;
	if ($totalCnt==0 && $line =~ m/^(ERR):/){
		die "Error in bamfile:\n\"$line\"\n";
	}
	my @sam = split/\t/,$line;
	if (@sam < 11){
		my $line2 = <STDIN>; print STDERR "Error:: malformed SAM (L${totalCnt}):\n$prevLine\n ** $line\n$line2\nAborting..\n";
		last;
	}
	my $refL = length($sam[9]);
	my $qualL = length($sam[10]);
	if ($refL != $qualL){
		print STDERR "DNA length != QUAL length!! \n(L" . $totalCnt-1 . ")  $prevLine\n(L${totalCnt}) ** $line\nAborting\n";
		last;
	}
	$totalCnt++;
	if ($sam[1] & 0x4 ){print "$line\n"; $totalNotMapped++; next;} #fail map
	my $fail=0; #assumme innocence
	my $xtrField = join("\t",@sam[11..$#sam]);
	my $mismatches = -1; my $gapopen=0;my$gapext=0; my $gapL=0;my $pid = -1;
	if ($bwt2sam){
		if ($xtrField =~ m/XM:i:(\d+)\s.*XO:i:(\d+)\s.*XG:i:(\d+)/){
			$mismatches = $1; $gapopen=$2; $gapext=$3;
			$gapL=$gapext+$gapopen;
			$pid = $mismatches/($refL-$gapL);
		} else {$bwt2sam=0;}
	} 
	if ($bwt2sam==0) {#mini2
		if ($xtrField =~ m/NM:i:(\d+).*de:f:([0-9\.]+)/){
		$pid = $2; $mismatches = $1;
		} else {$bwt2sam=1;}
		my $cig = $sam[5];chomp $cig;
		my @len = split (/\D+/,$cig); # storing the length per operation
		my @ops = split (/\d+/,$cig); # storing the operation
		shift @ops;
		foreach my $index(0..$#len){ #count insertions and delections from cigar string
			if($ops[$index] eq 'D' || $ops[$index] eq 'I'){ # deletions do affect the end position
				$gapL += $len[$index];
				#print STDERR "$ops[$index] $len[$index]\n";
			}
		}
		#die "$gapL\n\n$cig\n\n".scalar(@len)."  ".scalar(@ops)." $ops[0] $len[0]\n"; #ops emtpy
	}
	#print "$1 $2 $3 $refL ".$1/$refL." ".($2+$3)/$refL."\n";
	#95% id || 90% seq length
	
	
	
	if ( ($gapL)/$refL > $lengthC){$fail=1;$coverFail++;}
	if ($pid > $id){$fail=1;$idFail++;}
	if ($sam[4] < $MAQ){$mapscoreFail++;$fail=1;}
	#print $sam[4]." $pid\n";
	if ($fail){
		#set field to unmapped (4)
		#my $tmp =$sam[1];
		$sam[1] |= $failMask;
		$totalFail++;
		#die "$tmp   $sam[1]";
		#if ($sam[1] & 0x4 ){die "yes";} else {die "no\n";}
	}
	#join("\t",@sam)
	print "$line\n";
	$finalCnt++;
	$prevLine=$line;
	
}


#print STDERR "\n$prevLine\n$line\n";

print STDERR "BamFilter\nInentries: $totalCnt\nTotalRetained: $finalCnt\nTotalRm: $totalFail\nDue to\n  Mapping Qual (<$MAQ): $mapscoreFail\n  Coverage (<$lengthC): $coverFail\n  ID (>$id): $idFail\n  notMapped: $totalNotMapped\n";
#print STDERR "BamFilter: $totalCnt\nTotalFilter: $totalFail\nMapping Qual (<$MAQ): $mapscoreFail\nCoverage (<$lengthC): $coverFail\nID (>$id): $idFail\n";
sleep(1);
print STDERR "Done\n";

exit(0);