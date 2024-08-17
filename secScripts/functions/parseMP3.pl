#!/usr/bin/perl
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipwrite gzipopen);

my $inF = $ARGV[0];
my $outF = $ARGV[1];
my $mode = $ARGV[2];


my ($I,$OK) = gzipopen($inF,"mp3 output file [in]",1); 
my($O,$OK2) = gzipwrite($outF,"mp3 [out]",1);
my $cnt=0;
while (my $line = <$I>){
	next unless ($line =~ m/\d+\./);
	chomp $line;
	#if (length($line>1000)){{#newline not correctly set.. fucking mp3}
#	my @spl = split /\d+\./,$line2;
#	my $cnt=0;
#	foreach my $line (@spl){
		#print $line."\n";
	my @spl = split /\t/,$line;
	next unless (@spl > 1); 
	my $rep=0;
	my $len = scalar(@spl);
	while ($rep+6 < $len){
		#print "/";
		my $gID = $spl[$rep+1];
		chomp $gID; $gID =~ s/\s//g;
		#my $classi = "";
		#if ($spl[3] eq "Pathogenic" || $spl
		print $O "$gID\t$spl[$rep+6]\n";
		$rep+=8;
		$cnt++;
		#die if ($cnt==20);
	}
}
print "totalEntries=$cnt\n";
close $I; close $O;
