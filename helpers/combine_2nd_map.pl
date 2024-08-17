#!/usr/bin/env perl
#script that takes gene abundance from second mapping and creates abudnance matrix from this
#./combine_2nd_map.pl /ei/workarea/users/hildebra/projects/SaltTranscript/GlbMap/molly/ /ei/workarea/users/hildebra/projects/SaltTranscript/GlbMap/molly.mat

use warnings; use strict;
my $inD = $ARGV[0];
my $outM = $ARGV[1];
my %ABmat;my %smpls_all;
opendir (DIR, "$inD") or die $!;
while (my $f = readdir(DIR)){
	next unless (-e "$inD/$f" && $f =~ m/gz\.median\.pergene/);
	$f =~ m/(^.*)-0-smd\.bam\.coverage\.gz\.median\.pergene/;
	my $smpl = $1;
	$smpls_all{$smpl} = 1;
	open I,"<$inD/$f" or die $!;
	while (<I>){
		chomp;my @spl = split/\t/;
		$ABmat{$spl[0]}{$smpl} = $spl[1];
	}
	close I;
}
close DIR;

my @smpls = sort keys %smpls_all;
open O,">$outM";
print O "MapCmb";
foreach my $smpl (@smpls){
	print O "\t".$smpl;
}
print O "\n";
foreach my $gen (sort keys %ABmat){
	print O $gen;
	foreach my $smpl (@smpls){
		if (exists($ABmat{$gen}{$smpl})){
			print O "\t".$ABmat{$gen}{$smpl};
		} else {
			print O "\t0";
		}
	}
	print O "\n";
}
close O;