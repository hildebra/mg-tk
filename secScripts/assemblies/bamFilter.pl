#!/usr/bin/perl
use warnings;
use strict;


my $id = 0.05; my $queryCovMin = 0.8; my $MAQ=20; my $rmClpdRds=0;
$id = $ARGV[0] if (@ARGV >= 1); #max % dif nt id of alignment
$queryCovMin = $ARGV[1] if (@ARGV >= 2); #min length of alignment as % of query
$MAQ = $ARGV[2] if (@ARGV >= 3); #discard reads below MAQ threshhold (mapping qual)
$rmClpdRds = $ARGV[3] if (@ARGV >= 4); #1=remove reads clipped at both start and end; 0= don't;
my $idFail=0; my $mapscoreFail=0; my $coverFail=0; my $totalFail=0; my $clipFail=0;
my $finalCnt=0; #really the number of entries written..
my $failMask = 1 << 2;
my $bwt2sam=0; #update: don't use this any longer.. #switch to mini2, if not matching found..
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
	my $querL = length($sam[9]);
	my $qualL = length($sam[10]);
	if ($querL != $qualL){
		print STDERR "DNA length != QUAL length!! \n(L" . $totalCnt-1 . ")  $prevLine\n(L${totalCnt}) ** $line\nAborting\n";
		last;
	}
	$totalCnt++;
	if ($sam[1] & 0x4 ){print "$line\n"; $totalNotMapped++; next;} #fail map
	my $fail=0; #assumme innocence
	my $xtrField = join("\t",@sam[11..$#sam]);
	my $mismatches = -1;  my $gapDL=0;my $gapIL=0; my $matchL=0; my $clipL=0; #mismatches, gaps, matches and clips (nt) in alignments
	my $pid = -1; #inverse % nt id
	#my $queryCov = 0;
	if ($bwt2sam){
		my $gapopen=0;my$gapext=0;
		if ($xtrField =~ m/XM:i:(\d+)\s.*XO:i:(\d+)\s.*XG:i:(\d+)/){
			#$mismatches = $1; $gapopen=$2; $gapext=$3;$gapL=$gapext+$gapopen;$pid = $mismatches/($querL-$gapL);
		} else {$bwt2sam=0;}
	} 
	if ($bwt2sam==0) {#mini2
		#if ($xtrField =~ m/NM:i:(\d+).*de:f:([0-9\.]+)/){$pid = $2; $mismatches = $1;}# else {$bwt2sam=1;}
		my $cig = $sam[5];chomp $cig;
		my @len = split (/\D+/,$cig); # storing the length per operation
		my @ops = split (/\d+/,$cig); # storing the operation
		shift @ops;
		foreach my $index(0..$#len){ #count insertions and delections from cigar string
			if($ops[$index] eq 'M' ){
				$matchL += $len[$index];
			} elsif($ops[$index] eq 'D' ){ # Deletion (gap in the target sequence)
				$gapDL += $len[$index];
				#print STDERR "$ops[$index] $len[$index]\n";
			} elsif ($ops[$index] eq 'I'){ #	Insertion (gap in the query sequence)
				$gapIL += $len[$index];
			} elsif($ops[$index] eq 'H' || $ops[$index] eq 'S'){ #hard or soft clip
				$clipL += $len[$index];
			}
		}
			
		if ($pid < 0){
			
			if ($xtrField =~ m/NM:i:(\d+)/){$mismatches=$1;} #mismatches + gaps
			#if ($xtrField =~ m/de:f:([0-9\.]+)/){$pid = $1; }
			
			
			#$pid = $mismatches/($querL-$gapL);
			$pid = $mismatches/($querL-$clipL); #-$gapIL - $gapDL ?? -> gap counts as mismatch?
		} 
		#$queryCov = $matchL/$querL;
		#die "$gapL\n\n$cig\n\n".scalar(@len)."  ".scalar(@ops)." $ops[0] $len[0]\n"; #ops emtpy
	}
	#print "$1 $2 $3 $querL ".$1/$querL." ".($2+$3)/$querL."\n";
	#95% id || 90% seq length
	
	if ($rmClpdRds && $sam[5] =~ m/^(\d+)[HS].*\D(\d+)[HS]$/){#clipped alignment..
		#print STDERR " $1 $2 $sam[5] XX";
		if ($1 > $rmClpdRds && $2 > $rmClpdRds){
			$fail =1; $clipFail++;
			#print STDERR "$sam[5] ";
		}
	
	}
	
	
	#reasons for failing mapping:
	#if ( ($gapL)/$querL > $queryCov){$fail=1;$coverFail++;}
	if ( (1-($gapIL+$clipL)/$querL) < $queryCovMin){#new algo that looks only at clipped + gaps (I) regions
		$fail=1; $coverFail++;
		#print STDERR "(1-($gapIL+$clipL)/$querL) < $queryCovMin          ";
	}
	if ($pid > $id){
		$fail=1; $idFail++;
		#print STDERR "$sam[5] $pid $mismatches/($querL-$clipL)     ";
	}
	if ($sam[4] < $MAQ){
		$fail=1; $mapscoreFail++;
	}
	
	
	#print $sam[4]." $pid\n";
	if ($fail){
		#set field to unmapped (4)
		#my $tmp =$sam[1];
		$sam[1] |= $failMask;
		$totalFail++;
		#die "$tmp   $sam[1]";
		#if ($sam[1] & 0x4 ){die "yes";} else {die "no\n";}
	} else {$finalCnt++;}
	#join("\t",@sam)
	print "$line\n";
	
	$prevLine=$line;
}


#print STDERR "\n$prevLine\n$line\n";

print STDERR "BamFilter\nInentries: $totalCnt\nTotalRetained: $finalCnt\nTotalRm: $totalFail\nDue to\n  Mapping Qual (<$MAQ): $mapscoreFail\n  Coverage (<$queryCovMin): $coverFail\n  ID (>$id): $idFail\n  clipped reads rm ($rmClpdRds): $clipFail\n  notMapped: $totalNotMapped\n";
#print STDERR "BamFilter: $totalCnt\nTotalFilter: $totalFail\nMapping Qual (<$MAQ): $mapscoreFail\nCoverage (<$queryCovMin): $coverFail\nID (>$id): $idFail\n";
sleep(1);
print STDERR "Done\n";

exit(0);