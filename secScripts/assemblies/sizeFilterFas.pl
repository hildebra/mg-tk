#!/usr/bin/perl

use warnings;
use strict;

my $sizRestr = 500;my $sizRestr2 = 200;

if (@ARGV < 1){die"not enough args!\n";}

my $outFile1 = $ARGV[0].".filt";
my $outFile2 = $ARGV[0].".filt2";

if (@ARGV >1){
	$sizRestr = $ARGV[1];
}
if (@ARGV >2){
	$sizRestr2 = $ARGV[2];
}
if (@ARGV >3){
	$outFile1 = $ARGV[3]; $outFile2 = $outFile1."2";
}
my @singleFiles = split /,/,$ARGV[0];
open O,">",$outFile1;
open O2,">",$outFile2;
for (my $j=0;$j<@singleFiles;$j++){
	open I,"<",$singleFiles[$j];
	my $cnt =-1;
	my $head = ""; my $fas; my $fqMode = 0; my $curSize = 0; my $wrMode=0;#my $wrMode2=0;
	while (my $line = <I>){
		$cnt ++;
		if ($cnt == 0 && $line =~ m/^@/){$fqMode = 1;}
		
		if ($fqMode){#fastq reader
			if ($cnt % 4 == 0){#header
				$head = ">".substr($line,1);
			} elsif ($cnt % 4 == 1){$fas = $line;
				#check here, if size requirements are met
				if (length($fas) >= $sizRestr){
					print O $head.$fas; $wrMode=1; 
				} 
				#die $head.$fas;
				#and clean up
				$head=""; $fas="";
			}
			next;
		}
		
		
		if ($line =~ m/^>/){
			if ($sizRestr2 > 0 && !$wrMode && $curSize >= $sizRestr2){#check if falls into filter2
				print O2 $head.$fas; 
			}
			$head = $line; $curSize = 0; $wrMode=0;$fas="";
			next;
		}
		if ($wrMode){
			print O $line;
		} else {
			$fas .= $line;
			$curSize += length($line)-1;
			if ($curSize >= $sizRestr){
				print O $head.$fas; $wrMode=1; $head=""; $fas="";
			} 
		}
	}
	close I;
}

close O;close O2;
