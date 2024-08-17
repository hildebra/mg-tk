#!/usr/bin/perl
#./calcGC.pl /g/bork1/hildebra/SNP/GNMassSimu2/AssmblGrp_1/metag/scaffolds.fasta.filt  /g/bork1/hildebra/SNP/GNMassSimu2/AssmblGrp_1/metag/ContigStats/sgc.2
use warnings;
use strict;

use Mods::GenoMetaAss qw(readGFF);

sub evalGC;

my $inF = $ARGV[0];
my $outF = $ARGV[1];
my $isGenes = 0;
$isGenes = 1 if (@ARGV > 2);
#readGFF($ARGV[2]);


open I,"<$inF" or die "Can't open $inF";
open O,">$outF" or die "Can't open $outF";
print O "contig\tGC\n";

if ($isGenes){
	open O3,">${outF}3" or die "Can't open ${outF}3";
	print O3 "contig\tGC\n";
	my $ctgF = "${outF}3";
	$ctgF =~ s/\.pergene//;
	open OC3,">${ctgF}" or die "Can't open ${ctgF}";
	print OC3 "contig\tGC\n";
}

my $curTag = ""; my $GC=0; my $AT=0;
my $GC3=0; my $AT3=0; my $line = "";
my $GC3c=0; my $AT3c=0;my $curCtg="";
my $cnt=0;
while (my $lin = <I>){
	if ($lin =~ m/^>(.*)/){
		my $ct1= $1;
		if ($cnt > 0){
			my $isSameCtg =0;
			$isSameCtg =1 if ($ct1 !~ m/$curCtg/);
			#die "$isSameCtg  $curCtg $ct1\n" if ($isSameCtg);
			evalGC($isSameCtg);
			if ($isSameCtg){
				$curCtg = $ct1;
				$curCtg =~ s/_\d+$//;
			}
		} else {#ini contig
			$curCtg = $ct1;
			$curCtg =~ s/_\d+$//;
		}
		$curTag = $ct1; $GC=0;$AT=0; $line = "";
		
		
		
		$cnt++; 
		next;
	}
	chomp $lin; $line .= $lin;
}
evalGC(1);
close O; close I;
if ($isGenes){
close O3 ;close OC3 ;
}



exit(0);




sub evalGC($){
	my ($difCtg) = @_;
	$GC += ($line =~ tr/G//);	$GC += ($line =~ tr/C//);
	$AT += ($line =~ tr/A//);	$AT += ($line =~ tr/T//);
	if ($GC + $AT == 0){
		print O "$curTag\t-1";
	} else {
		print O "$curTag\t". sprintf('%.3f', ($GC/($GC+$AT)*100)) ."\n";
	}
	#get every third char
	if ($isGenes){
		my $line2="";
		for (my $c=2;$c<length($line);$c+=3){$line2 .= substr $line, $c,1;}
		if ((length($line)%3) != 0){die "$curTag $line $line2\n";}
#		die "$line\n\n$line2\n" ."\n";
		$GC3 += ($line2 =~ tr/G//);	$GC3 += ($line2 =~ tr/C//);
		$AT3 += ($line2 =~ tr/A//);	$AT3 += ($line2 =~ tr/T//);		
		if ($GC3 + $AT3 == 0){
			print O3 "$curTag\t-1";
		} else {
			print O3 "$curTag\t". sprintf('%.3f', ($GC3/($GC3+$AT3)*100)) ."\n";
		}
		$GC3c+=$GC3;$AT3c+=$AT3;
		$GC3=0;$AT3=0;
	
		if ($difCtg){#different contig
			if ($GC3c + $AT3c == 0){
				print OC3 "$curTag\t-1";
			} else {
				print OC3 "$curCtg\t". sprintf('%.3f', ($GC3c/($GC3c+$AT3c)*100)) ."\n";
			}
			$GC3c=0;$AT3c=0;
		}
	}

}


















