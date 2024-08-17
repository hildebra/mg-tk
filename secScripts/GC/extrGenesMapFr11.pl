#!/usr/bin/env perl
#takes mapping of reads to a gene catalog and identifies genes that are mapped to, i.e. candidates to include in gene catalog

use warnings;
use strict;
#use Mods::GenoMetaAss qw(readMap qsubSystem emptyQsubOpt );
#use Mods::IO_Tamoc_progs qw(jgi_depth_cmd inputFmtSpades createGapFillopt getProgPaths);
my %collctGns;

my $doStep1=0; #extraction of high coverage genes from bam mapping files
my $doStep2=1; #sdm based extraction of genes, that had high coverage
my $at = 0;
my @nms = ("frz11","Tara","IGC");
my @fnaRefs = ("/g/bork3/home/hildebra/DB/freeze11/freeze11.genes.representatives.fa","/g/bork3/home/hildebra/DB/GeneCats/Tara/Tara.fna","/g/bork3/home/hildebra/DB/GeneCats/IGC/1000RefGeneCat.fna");
my @faaRefs = ("/g/bork3/home/hildebra/DB/freeze11/freeze11.proteins.representatives.fa","/g/bork3/home/hildebra/DB/GeneCats/Tara/Tara.faa.gz","/g/bork3/home/hildebra/DB/GeneCats/IGC/1000RefGeneCat.faa");

my $inD = "/g/scb/bork/hildebra/Tamoc/FinSoil/GlbMap/$nms[$at]/";#Tara/";#frz11/ IGC/
my $fnaRef = $fnaRefs[$at];
my $faaRef = $faaRefs[$at];
my (@FList) = glob($inD."/*.jgi.depth.txt.gz");
if ($doStep1){
	open L,">$inD/geneLogs.txt";
	opendir D,$inD or die "Can't open dir $inD\n";
	my $cnt=0;
	while (my $dd = readdir D){
		if (-e $inD.$dd  && $dd =~ m/.jgi.depth.txt.gz$/){
			my $gnHit = 0; my $lncnt=0; my $bpCovTot=0;
			#print $dd."   ";
			open I,"gunzip -c $inD/$dd | " or die "can't open $inD/$dd\n";
			while (my $line = <I>){
				$lncnt++;
				next if ($lncnt == 1);
				my @spl = split /\t/,$line;
				if (0.8*250/$spl[1] < $spl[2]){
					$collctGns{$spl[0]} ++;
					$gnHit++;
					my $bpcov = $spl[2]*$spl[1];
					$bpCovTot += $bpcov;
				}
			}
			
			close I;
			print "$dd\t$gnHit\t$bpCovTot\n";
			print L "$dd\t$gnHit\t$bpCovTot\n";
			
		}
		$cnt++;
		if ($cnt > 10000){die "Abort while1\n";}
	}
	closedir D;
	close L;

	open O,">$inD/CoveredGenes.txt";
	foreach my $gn (keys %collctGns){
		print O $gn."\n";
	}
	close O;

}

if ($doStep2){
	my $cmd = "";
	#$cmd .= "/g/bork3/home/hildebra/dev/C++/sdm/./sdm -i_fna $fnaRef -o_fna $inD/extr.rds -paired 1 -specificReads $inD/CoveredGenes.txt -log nolog\n";
	$cmd .= "/g/bork3/home/hildebra/dev/C++/sdm/./sdm -i_fna $faaRef -o_fna $inD/extr.AA -paired 1 -specificReads $inD/CoveredGenes.txt -log nolog\n";
	print $cmd."\n";
	system $cmd."\n";
}
