#!/usr/bin/perl
#script to filter out all MGS that have > XX genes attached to them

use strict;
use warnings;


sub readCluster{
	my ($cF) = @_;
	my %canos;
	open I,"<$cF" or die "can;t open $cF\n";
	while (my $line = <I>){
		chomp $line;
		my @spl = split /\t/,$line;
		push(@{$canos{$spl[0]}}, $spl[1]);
	}
	close I;
	return \%canos;
}

sub readProfile{
	my ($pF) = @_;
	my %simplePro;
	open I,"<$pF" or die $!;
	while (my $line = <I>){
		chomp $line;
		$line =~ m/^(\S+)\s/;
		$simplePro{$1} = $line;
		#print "$1\n"; 
	}
	close I;
	return \%simplePro;
}

my ($clusterF, $profileF, $newProfF, $XX)  = @ARGV;

my $hr = readCluster($clusterF);
my %canos = %{$hr};
$clusterF  =~ m/(^.*\/)[^\/]+$/;
my $idir = $1;


my $numSize=0;
my @selPros;
open O,">$idir/MGS_size.txt" or die $!;
foreach my $can (sort keys %canos){
	if (scalar(@{$canos{$can}}) >= $XX){
		$numSize++ ;
		push(@selPros, $can);
	}
	print O "$can\t".scalar(@{$canos{$can}})."\n";
}
close O;

#since profiles are ordered by size, I only need to get the first X profiles
#system "head -n$numSize $profileF > $newProfF";
$hr = readProfile($profileF);
my %Pro = %{$hr};

open O,">$newProfF" or die $!;
foreach my $pr (@selPros){
	die "can't find profile $pr" unless (exists($Pro{$pr}));
	print O "$Pro{$pr}\n";
}
close O;

print "Extracted $numSize profiles as guides into $newProfF\n";

