#!/usr/bin/env perl
#creates GTDB tax per single gene
#perl reformatProGenomes.pl /hpc-home/hildebra/DB/MarkerG/proGenomes3/markerGenes/ /hpc-home/hildebra/DB/MarkerG/proGenomes3/proGenomes3_specI_clustering.tab /hpc-home/hildebra/DB/MarkerG/proGenomes3/proGenomes3_specI_lineageGTDB.tab /hpc-home/hildebra/DB/MarkerG/proGenomes3/specI.pergene.tax
#$SpecID $inSImap $GTDBspecI $taxPerGene

use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( readClstrRev systemW median readFastHD);

my $SpecID = $ARGV[0];
my $inSImap = $ARGV[1];
my $GTDBspecI = $ARGV[2];
my $taxPerGene = $ARGV[3];

my %FMGcutoffs = (COG0012=>94.8,COG0016=>95.8,COG0018=>94.2,COG0172=>94.4,COG0215=>95.4,COG0495=>96.4,COG0525=>95.3,COG0533=>93.1,COG0541=>96.1,
COG0552=>94.5,COG0048=>98.4,COG0049=>98.7,COG0052=>97.2,COG0080=>98.6,COG0081=>98,COG0085=>97,COG0087=>99,COG0088=>99,COG0090=>98.8,COG0091=>99,
COG0092=>99,COG0093=>99,COG0094=>99,COG0096=>98.6,COG0097=>98.4,COG0098=>98.7,COG0099=>98.9,COG0100=>99,COG0102=>99,COG0103=>98.4,
COG0124=>94.5,COG0184=>98.2,COG0185=>99,COG0186=>99,COG0197=>99,COG0200=>98.4,COG0201=>97.2,COG0202=>98.4,COG0256=>99,COG0522=>98.6);

die "no dir given" unless (-d $SpecID);

my @allgenes;
foreach my $cog (keys %FMGcutoffs){
	my $ar = readFastHD("$SpecID/$cog.fna");
	push(@allgenes, @{$ar});
}

if (1){
	print "Reading GTDB taxids\n";
	open I,"<$GTDBspecI" or die $!;
	if (0){ #; separated GTDB tax
		while (<I>){
			chomp; my @spl = split /\t/;
			my $gt1 = $spl[0]; $gt1 =~ s/[dpcofgs]__//g;
			my @gt = split /;/,$gt1;
			#die "@gt\n";
			$GTDB{$spl[1]}=\@gt;
		}
	} else {
		while (<I>){
			chomp; my @spl = split /\t/;
			my $gt1 = shift @spl; shift @spl; #remove second entry
			#$gt1 =~ s/[dpcofgs]__//g;
			#my @gt = split /;/,$gt1;
			#die "@gt\n";
			$GTDB{$gt1}=\@spl;
		}
	}
	close I;
}

my %gene2SI;
open I,"<$inSImap" or die $!;
while (<I>){
	chomp; my @spl = split /\t/;
	#$spl[1] =~ m/^([^;]+)/;
	my @s2 = split /;/,$spl[1];
	foreach (@s2){
		$gene2SI{$_} = $spl[0];
	}
}
close I;

my @taxAbbr = ("k__","p__","c__","o__","f__","g__","s__");

my $noHit=0;my$hasHit=0;

open O,">$taxPerGene" or die $!;
foreach my $gene (@allgenes){
	my $gen1 = $gene;
	$gen1 =~ s/\.[^\.]+$//;
	unless (exists($gene2SI{$gen1})){
		print "can't find specI: $gene $gen1\n" ;
		$noHit++;
		next;
	}
	die "can't find Genome: ".$gene2SI{$gen1} ." $gene $gen1\n" unless (exists($GTDB{$gene2SI{$gen1}} ));
	my $taxStr = ""; my @locTax = @{$GTDB{$gene2SI{$gen1}}};
	for (my $i=0; $i<7;$i++){
		$taxStr .= $taxAbbr[$i].$locTax[$i].";";
	}
	#die "$gene\t$taxStr\n";
	print O "$gene\t$taxStr\n";
	$hasHit++;
}
close O;


print "Finished: $hasHit with hits, $noHit no hits\n";


