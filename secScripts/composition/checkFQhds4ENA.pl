#!/usr/bin/perl
#USAGE: checkFQhds4ENA.pl [fq[.gz]] [1,2,3]
#checks if the headers in 1st arg fastq files are >256 chars (not allowed for ENA); 2nd arg is whether this is read 1/2 or singleton (3)

use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipopen gzipwrite);

my $inF = $ARGV[0];
my $readP = $ARGV[1];

my $LTHR = 255;
my $LCUT = $LTHR -5;
my $doRewrite = 0;


my ($IN,$OK) = gzipopen($inF,"fastq file 4 upload");
if (!$OK){die "something wrong when opening $inF\n";}
my $lcnt  = 0;

#first round: check if any header needs to be changed?? -> faster than rewriting..
while (<$IN>){
	if ($lcnt % 4 == 0){
		if (substr($_,0,1) ne "@"){die "Something wrong with fastq format: Expected \"@\" at line $lcnt, found $_\n";} 
		if (length($_) >= $LTHR){
			print "$inF needs rewriting, found fastq header with length ". length($_) . "\n";
			$doRewrite = 1;
			last;;
		}
	}
	$lcnt++;
}
close $IN;

#nothing needs to be done..
if (!$doRewrite){
	print "file $inF passed header check\n";
	exit (0) ;
}

#something needs to be done, rewriting file..
my ($IN2,$OK2) = gzipopen($inF,"fastq file 4 upload");
if (!$OK2){die "something wrong when opening $inF\n";}
my $tmpF = $inF.".tmp.gz";
my $OUT = gzipwrite($tmpF,"temp fastq file 4 upload");


$lcnt  = 0;
while (<$IN2>){
	if ($lcnt % 4 == 0){
		if (substr($_,0,1)  ne "@"){die "Something wrong with fastq format: Expected \"@\" at line $lcnt, found $_\n";} 
		if (length($_) >= $LTHR){
			my $newHD = substr($_,0,$LCUT);
			if ($readP == 1){$newHD.= "/1";
			} elsif ($readP == 2){$newHD.= "/2";
			}
			print $OUT $newHD."\n";
		} else {
			print $OUT $_;
		}
	} else {
		print $OUT $_;
	}
	$lcnt++;
}


close $IN2;close $OUT;

system "rm -f $inF; mv $tmpF $inF;";

print "Done rewriting\n";
