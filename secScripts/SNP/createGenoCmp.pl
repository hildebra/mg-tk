#!/usr/bin/env perl
#compares two (long) sequences / genomes and finds insertions / deletions
#further also links to 
#usage: ./eval_assembly.pl [assembly.fa] [ref.fa] [num cores] [id_of_cmp]
use warnings;
use strict;
use Mods::GenoMetaAss qw(readFasta);
use Mods::IO_Tamoc_progs qw(getProgPaths);

#my $blatBin = getProgPaths("blat");
my $pigzBin = getProgPaths("pigz");
my $nucmP = "/g/bork3/home/hildebra/bin/MUMmer3.23/";# getProgPaths("nucmer");
my $inD1 = "/g/bork3/home/hildebra/Collab/AlexSNP/RefHMPmock";
my $inD2 = "/g/bork3/home/hildebra/Collab/AlexSNP/RefHMPwrongStrain";
my $outD = "/g/bork3/home/hildebra/Collab/AlexSNP/RefHMPmock/refGdiff/";
system "mkdir -p $outD" unless (-d $outD);

my @inRefG = ("Msmithii35061/640427121.fna","ABaumannii17978/640069301.fna","Aodontolyticus17982/640963058.fna",
			"Bcereus10987/637000016.fna","Bvulgatus8492/640753008.fna","Efaecalis47077/2511231133.fna","Sepidermis12228/637000281.fna");
my @inQueG = ("Msmithii2375/643886144.fna","ABaumaniiiMDRTJ/2513237250.fna","AodontolyticusF0309/647000206.fna",
			"BcereusNC7401/2511231081.fna","BvulgatusPC510/647000214.fna","EfaecalisERV68/2519103141.fna","SepidermidisRP62A/637000282.fna");
my @outNames = ("Msmithii35061","ABaumannii17978","Aodontolyticus17982","Bcereus10987","Bvulgatus8492","Efaecalis47077","Sepidermis12228");
			

for (my $i=0;$i<@outNames;$i++){
	my $ref = "$inD1/$inRefG[$i]";
	my $query = "$inD2/$inQueG[$i]";
	my $outF = $outNames[$i];
	my $tmpDelta = "$outD/$outF.delta";
	my $cmd = "$nucmP/nucmer -maxmatch -c 60 -p $outD/$outF $ref $query \n";

	$cmd .= "$nucmP/show-coords -T -r -c -l $tmpDelta > $outD/$outF.coords\n";
	$cmd .= "$nucmP/show-snps -T -C $tmpDelta > $outD/$outF.snps\n";
	$cmd .= "$nucmP/show-tiling $tmpDelta > $outD/$outF.tiling\n";

	system $cmd;
}