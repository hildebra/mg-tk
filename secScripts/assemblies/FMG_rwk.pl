#!/usr/bin/env perl
#creates FMGids.txt file with all fna names from the marker genes
use warnings;
use strict;
#./FMG_rwk.pl /g/scb/bork/hildebra/SNP/GNMass2_singl/alien-11-374-0/assemblies/metag/ContigStats/FMG/

my $inD1 = $ARGV[0];
my $inD = $inD1."/temp/";

opendir(DIR, $inD) || die "can't opendir $inD: $!";
my @FMGs = grep { /\.IDs\.txt/ && -f "$inD/$_" } readdir(DIR);
closedir DIR;
open O,">$inD1/FMGids.txt";
foreach my $FM ( @FMGs){
	#print $FM."\n";
	my $ID = $FM; $ID =~ s/\.IDs\.txt//;
	open I,"<$inD$FM" or die "Cant open $inD$FM";
	while (my $l=<I>){chomp($l);print O "$l $ID\n";}
	close I;
}
close O;

