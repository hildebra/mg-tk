#!/usr/bin/env perl
#takes a list of genes and finds contigs in each sample, on which these genes occur
#USAGE: ./geneListSameAssembly.pl [mode] [geneListFile.txt] [path to gene catalog] [output path]
#./geneListSameAssembly.pl  stat /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/TEC2MGS/TC0044.can.genes /g/scb/bork/hildebra/SNP/GCs/T2_GNM3_ABR MM3
#./geneListSameAssembly.pl  stat /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/R_filt/filter.genes.txt /g/scb/bork/hildebra/SNP/GCs/T2_GNM3_ABR /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/R_filt/contigs MM7
use warnings;
use strict;

use Mods::GenoMetaAss qw(readClstrRev readMapS);
my $samBin = "/g/bork5/hildebra/bin/samtools-1.2/samtools";
use Mods::GenoMetaAss qw(readTabByKey systemW);





die "Not enough input args given!\n" if (@ARGV < 4);
my $mode = $ARGV[0]; my $inF = $ARGV[1]; 
if (!-f $inF){die "Can;t find input gene file $inF\nAborting\n";}
my $GCd = $ARGV[2]; 
my $oDir = $ARGV[3];
my $refSmpl ="";$refSmpl = $ARGV[4] if (@ARGV>4); 
my $refSmplDir = "";

my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;
my ($hr,$hr2) = readMapS($mapF);
my %map = %{$hr};
my $tamocDir = $map{outDir}; $tamocDir =~ s/,.*//;

if ($refSmpl ne ""){
	if (!exists $map{$refSmpl}{wrdir}){die "Can't find refsample $refSmpl in map!\n";}
	$refSmplDir = $map{$refSmpl}{wrdir};
	print "Processing contigs of ".$refSmplDir."\n";
}
$oDir .= "/" unless ($oDir =~ /\/$/);
system "mkdir -p $oDir";

my %Conts;
my @tgenes;#genes in GC that are somehow binned
open I,"<$inF";
while (my $line = <I>){
	$line =~ s/\r\n?/\n/; chomp $line; 
	push (@tgenes,$line);
}
close I;
if (@tgenes == 1){#only one entry, check for , separation
	@tgenes = split /[,\s+]/,$tgenes[0];
}
#read gene list


#die "@tgenes\n";
my ($hr1s,$hr2s) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx");
my %gene2cl = %{$hr1s}; my %cl2gene = %{$hr2s};

my @tmp = keys %cl2gene; #print "$tmp[1] $tmp[1024]\n";
#print $cl2gene{1058}."\n".$cl2gene{"1058"}."\n";

print "Processing genes.. \n"; my $geneCnt=0;my $geneCnt1=0;
foreach my $genNum (@tgenes){
	#MM3__C10_L=109019;_51
	#$genNum = int $genNum;
	die "can't find gene \"$genNum\"\n" unless (exists $cl2gene{$genNum});
	my $genLst = $cl2gene{$genNum};
	foreach my $gen (split (/,/,$genLst)){
		$gen =~ s/^>//; $geneCnt++;
		#print $gen."\t$genNum\n";
		my @splX = split(/__/,$gen); #this is the sample
		my @spl = split(/;_/,$splX[1]); #this is the contig
		#print $splX[0]."\n" unless (exists ($Conts{$splX[0]}));
		if (exists ($Conts{$splX[0]}{$spl[0]})){	push (@{$Conts{$splX[0]}{$spl[0]}},$spl[1]);
		} else { $Conts{$splX[0]}{$spl[0]} = [$spl[1]]; }
	}
	$geneCnt1++;
}
print "Done with $geneCnt / $geneCnt1 genes\n";
die "can;t find refsmpl $refSmpl\n" if ($refSmpl ne "" && !exists $Conts{$refSmpl});

#grand summary file
open S,">$oDir/SummaryStats.txt";
#"$refSamp\t".@allCtgs."\t$ctgCngCnt\t$ctgCngCnt2\t$missedGenes\t$missedGenes2\n";
print S "Smpl\t#Ctgs\t#genes\t#genesFilt\t#missed\t#missedFilt\n";
#extract contig fastas
if ($refSmpl ne ""){
	my $scaffs = $refSmplDir . "assemblies/metag/scaffolds.fasta";
	my @allCtgs = keys %{$Conts{$refSmpl}};
	for (my $i=0;$i<@allCtgs;$i++){$allCtgs[$i] = "'".$refSmpl."__".$allCtgs[$i].";'";}
	my $cmd = "$samBin faidx $scaffs ". join (" ", @allCtgs) . " > $oDir/$refSmpl.ctgs.fna";
	#print $cmd."\n";
	system $cmd;
} else {
	foreach my $refSamp (sort keys %Conts){
		my @allCtgs = keys %{$Conts{$refSamp}};
		my $refSamp2 = $refSamp;
		if (exists(  $map{altNms}{$refSamp}  )){$refSamp2 = $map{altNms}{$refSamp}; print "Nx "}
		if (!exists($map{$refSamp2})){die "can't find sample $refSamp\n";}
		my $SmplDir = $map{$refSamp2}{wrdir};
		my $metaGD = `cat $SmplDir/assemblies/metag/assembly.txt`; chomp $metaGD;

		my $scaffs = $metaGD . "/scaffolds.fasta";
		
		
		#stats and determine if this contig is even taken..
		#count number of genes in Sample associated to MGS
		my %geneCnts = readTabByKey("$metaGD/genePred/genes.per.ctg");
		open O,">$oDir/ContigList_$refSamp.txt";
		my $ctgCngCnt=0; my $missedGenes = 0;my $ctgCngCnt2=0; my $missedGenes2 = 0; my $frac=0;
		#foreach my $spCtg (@allCtgs){
		my $idx =0;
		while ($idx <= $#allCtgs){
			my $spCtg = $allCtgs[$idx]; 
		#die @tmp . " $geneCnts{$tmp[0]} $tmp[0] $spCtg\n";  MM350
			my $ctgKey = $refSamp."__".$spCtg.";";
			@{$Conts{$refSamp}{$spCtg}} = sort {$a <=> $b}(@{$Conts{$refSamp}{$spCtg}});
			my @tmp = @{$Conts{$refSamp}{$spCtg}};
			die "Can't find contig $ctgKey in length ref\n" if (!-exists $geneCnts{$ctgKey});
			my $maxGenes = $geneCnts{$ctgKey};
			$frac = @tmp / ($maxGenes+1);
			#if (@tmp != ($maxGenes) ){print STDERR @tmp ." != $maxGenes\n";}
			$ctgCngCnt += @tmp;
			$missedGenes += ($maxGenes+1) - @tmp ;
			
			#$allCtgs[$idx] = "'".$ctgKey.";'"; 
			my $rmTag = "";
			if ($frac < 0.5){splice @allCtgs,$idx,1; $rmTag="X";
			} else {$idx ++; $ctgCngCnt2 += @tmp; $missedGenes2 += ($maxGenes+1) - @tmp;}
			print O "$rmTag$spCtg\t".@tmp."\t" . ($maxGenes+1) ."\t".$frac."\t". join(",",@tmp) . "\n";
		}
		print  "".@allCtgs." Contigs, $ctgCngCnt($ctgCngCnt2) genes, $missedGenes($missedGenes2) missed genes in $refSamp\n";
		print S "$refSamp\t".@allCtgs."\t$ctgCngCnt\t$ctgCngCnt2\t$missedGenes\t$missedGenes2\n";
		close O;
		
		#actual contig extraction from fasta
		for (my $i=0;$i<@allCtgs;$i++){$allCtgs[$i] = "'".$refSamp."__".$allCtgs[$i].";'";}
		my $cmd = "$samBin faidx $scaffs ". join (" ", @allCtgs) . " > $oDir/$refSamp.ctgs.fna";
		#print $cmd."\n";
		system $cmd;

		
		
	}
}
close S;

foreach my $refSamp (sort keys %Conts){
last; #already done, moved to upper loop
	my @allCtgs = keys %{$Conts{$refSamp}};
	open O,">$oDir/ContigList_$refSamp.txt";
	#count number of genes in Sample associated to MGS
	my $ctgCngCnt=0; my $missedGenes = 0;
	#guilty by assembly
	my $refSamp2 = $refSamp;
	if (exists(  $map{altNms}{$refSamp}  )){$refSamp2 = $map{altNms}{$refSamp}; print "Nx "}
	if (!exists($map{$refSamp2})){die "can't find sample $refSamp\n";}
	my $SmplDir = $map{$refSamp2}{wrdir};
	my $metaGD = `cat $SmplDir/assemblies/metag/assembly.txt`; chomp $metaGD;

	my %geneCnts = readTabByKey("$metaGD/genePred/genes.per.ctg");
	#my @tmp = keys %geneCnts;
	#my $geneCnts = `cat $SmplDir/assemblies/metag/genePre/genes.per.ctg`;
	foreach my $spCtg (@allCtgs){
	#die @tmp . " $geneCnts{$tmp[0]} $tmp[0] $spCtg\n";  MM350
		my $ctgKey = $refSamp."__".$spCtg.";";
		@{$Conts{$refSamp}{$spCtg}} = sort {$a <=> $b}(@{$Conts{$refSamp}{$spCtg}});
		my @tmp = @{$Conts{$refSamp}{$spCtg}};
		die "Can't find contig $ctgKey in length ref\n" if (!-exists $geneCnts{$ctgKey});
		my $maxGenes = $geneCnts{$ctgKey};
		#if (@tmp != ($maxGenes) ){print STDERR @tmp ." != $maxGenes\n";}
		print O "$spCtg\t".@tmp."\t" . ($maxGenes+1) ."\t". join(",",@tmp) . "\t". join(",",@allCtgs). "\n";
		$ctgCngCnt += @tmp;
		$missedGenes += ($maxGenes+1) - @tmp ;
	}
	print  "".@allCtgs." Contigs, $ctgCngCnt genes, $missedGenes missed genes in $refSamp\n";
	
	close O;
	#die();
}















