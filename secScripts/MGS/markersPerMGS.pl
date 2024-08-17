#!/usr/bin/env perl
#script to extract all marker genes associated to single MGS. Used to later make sure each marker gene is only used once in abundance calculations
use warnings;
use strict;
use Data::Dumper; 
use Getopt::Long qw( GetOptions );
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(systemW gzipopen);
use Mods::Binning qw (readMGSrevRed);


my $MAGrep = ""; my $GCdir = "";
my $useGTDBmg = "GTDB"; my $Binner = "SB";
my $numCor = 4; my $MGStaxF = "";
my $GC_bin_dir = "";

GetOptions(
	"GCd=s"      => \$GCdir,
	"binD=s"     => \$GC_bin_dir,
	"MAGlogFile=s" => \$MAGrep, 
	"MGset=s" => \$useGTDBmg,#GTDB or FMG
	"Binner=s" => \$Binner,
	"cores=i" => \$numCor,
	"MGStaxFile=s" => \$MGStaxF,
	
);

die "-GCd needs to be specified!\n" if ($GCdir eq "");


#dir in GCdir where marker genes are stored..
my $COGdir = "FMG";
if ($useGTDBmg eq "GTDB"){ 	$COGdir = "GTDBmg";}
$GC_bin_dir = "$GCdir/Bin_SB" if ($GC_bin_dir eq ""); 
my $markerGdir  = "$GCdir/$COGdir/";
my $outDir = "$GC_bin_dir/Annotation/";
my $outFile = "$outDir/marker2MGS.txt";
my $outFile2 = "$outDir/marker2MGS.LCA.txt";
my $primaryClusF = "$GC_bin_dir/SB.clusters.core"; 
my $matrOutDir = "$outDir/Abundance/";
$MGStaxF = "$outDir/GTDBTK.tax" if ($MGStaxF eq "");
system "mkdir -p $matrOutDir" unless (-d $matrOutDir);


if ($MAGrep eq ""){$MAGrep = "$GC_bin_dir/LOGandSUB/MAGvsGC.txt.gz"; print "Inferring MAGMGS file at $MAGrep\n";}
#my $outFile = "$markerGdir/marker2MGS.txt";


my $clusHR = readMGSrevRed($primaryClusF); my %oriG2MGS  = %{$clusHR};

#first read in MGs / MGS
my $lcnt=0;my $mgsCnt=0;my @MGcats; my %MGinMGS;
my ($I,$OK) = gzipopen($MAGrep,"MAG bin report",1);

while (my $line = <$I>){
	chomp $line; 
	$lcnt++;
	my @spl = split /\t/,$line;
	if ($lcnt==1){
		die "Not correct header in $MAGrep:\n$line\n" unless ($spl[0] eq "MAG" && $spl[7] eq "Contamination");
		#remove: MAG     MGS     Representative4MGS      Match2MGS       Uniqueness      AssociatedMGS   Completeness    Contamination   LCAcompleteness N50     N_Genes  GC      CodingDensity   CentreScore     CompoundScore   Domain  Phylum  Class   Order   Family  Genus   Species   & other_genes
		for (my $i=0;$i<22;$i++){shift @spl;} pop @spl; 
		@MGcats = @spl;
		print "@MGcats\n\n";
		next;
	}
	die "header not found! \n" if (@MGcats == 0);
	my $MGS = $spl[1];
	next if ($MGS eq "?");
	$mgsCnt++;
	#die "$spl[0]  $spl[5]\n";
	my $maxCat = scalar(@MGcats)+22;
	for (my $i=22; $i < $maxCat; $i++){
		next unless defined $spl[$i];
		my @spl2 = split/,/,$spl[$i];
		foreach my $gene (@spl2){
			$MGinMGS{$MGS}{$i}{$gene} ++;
		}
		#print "$spl[$i]  ";
	}
	#die;
}
close $I;

my @MGSids = sort(keys(%MGinMGS));

print "Found ". ($lcnt-1) ." MAGs, $mgsCnt associated to ". @MGSids . " MGS.\n";
#die "$outFile\n";

#second: collate into single file
my $totgenes = 0; my $totMGS = 0; my $catPres=0; my $cntCatPres=0;
my $geneAss2MGS =0; my $geneNotAss2MGS=0;
my %gene2MGS;  #marks genes that were clearly in the unique core to an MGS
my %geneDirtyMGS;#marks genes that are in principle in MGS, but were unclear/not unique -> should be removed!
my @cntPerCat;
for (my $i=0; $i < (@MGcats); $i++){ $cntPerCat[$i]=0;}
open O,">$outFile" or die "Can't open outfile $outFile\n";
print O "MGS\t" . join("\t",@MGcats)."\n";
foreach my $MGS (@MGSids){
	print O "$MGS";
	for (my $i=0; $i < (@MGcats); $i++){ 
		my @genes = keys %{$MGinMGS{$MGS}{$i}};
		next unless(scalar(@genes));
		$totgenes += scalar(@genes);
		$catPres++;$cntCatPres+= scalar(@genes);
		$cntPerCat[$i] += scalar(@genes);
		foreach (@genes){ 
			#needs to control for what is found in reference MGS clustering to attain clean profile..
			if (exists($oriG2MGS{$_}{$MGS})){
				push(@{$gene2MGS{$_}}, $MGS);  $geneAss2MGS++;
			} else {
				$geneDirtyMGS{$_} = 1;
				$geneNotAss2MGS++;
				#die "$_  $MGS not found!\n";
			}
		}
		#die "@genes  $totgenes\n";
		print O "\t".join(",",@genes); 
	}
	print O "\n";
	$totMGS++;
}
close O;
for (my $i=0; $i < (@MGcats); $i++){ 
	$cntPerCat[$i] /=scalar(@MGSids);
}

open Ox,">$outDir/MarkerGeneStats.txt" or die $!;
print Ox "GeneCat	avgOccPerMGS\n";
for (my $i=0; $i < (@MGcats); $i++){ 
	print Ox "$MGcats[$i]	$cntPerCat[$i]\n";
}
close Ox;

@cntPerCat = sort(@cntPerCat);
for (my $i=0; $i < (@MGcats); $i++){ 
	print int(10*$cntPerCat[$i])/10 . " ";
}

print "Done, on average ". int($totgenes/$totMGS+0.5) ." genes in $totMGS MGS, ". int(100*$totgenes/$totMGS/scalar(@MGcats)+0.5)/100 ." genes per category, " . int(100*$cntCatPres/$catPres+0.5)/100 . " if gene present. Wrote MGs in MGS to $outFile\n";
print "Found $geneAss2MGS genes assigned to MGS clustering, $geneNotAss2MGS genes ambigous and removed\n";


#combine with LCA.. from per MG .LCA files..
#readMG_LCA($markerGdir,\@MGcats);
my $LCAfiles=0; my $LCAentries=0; my %LCA; my $LCAhead="";
foreach my $GTcat (@MGcats){
	my $GTLCA = "$markerGdir/$GTcat.LCA";
	next unless (-e $GTLCA);
	$LCAfiles++; my $lcnt=0;
	open I ,"<$GTLCA" or die $!;
	while (<I>){
		$lcnt++;
		chomp; 
		my @spl = split /\t/;
		my $id = shift @spl;
		if ($lcnt==1){
			$LCAhead = join(";",@spl);next;
		}
		#m/(^\S+)\s/;
		while (@spl < 7){push(@spl,"?");}
		$LCA{$id} = join(";",@spl);
		$LCAentries++;
	}
	close I;
}

#read MGS tax
my %MGStax;
open I,"<$MGStaxF" or die $!;
while (<I>){
	chomp; my @spl = split /\t/;
	#$MGStax{$spl[0]} = $spl[1];
	my $id = shift @spl;
	$spl[0] = "?" if ($spl[0] =~ m/Unclassified Bacteria/);
	@spl = split /;/,$spl[0];
	while (@spl < 7){push(@spl,"?");}
	$MGStax{$id} = join(";",@spl);
}
close I;
print "Read " . scalar(keys %MGStax) . " MGS tax annoations\n"; 


#write a file with the genes associated to MGS and their taxonomy. Also add other marker genes for completeness
my $wrLCA=0; my $link2MGS=0; my $ambGenes=0;
my @mgkeys = sort(keys %LCA);


#now compare how many counts of an LCA based tax were found..
my %newMGStax;
foreach my $mark (@mgkeys){
	if (exists($gene2MGS{$mark})){
		my @MGSs = @{$gene2MGS{$mark}};
		foreach my $MGSl (@MGSs){
			$newMGStax{$MGSl}{$LCA{$mark}}++;
		}
	}
}


my $replTax =0;
#check if MGStax needs update..
foreach my $MGSl (keys %newMGStax){
	#$MGSl = "MGS.1";
	#print "$MGSl  : $MGStax{$MGSl}  \n";
	if ($MGStax{$MGSl} =~ m/\?$/){
		#print "Yep\n";
		my @cnts = values %{$newMGStax{$MGSl}};
		my @keys = keys %{$newMGStax{$MGSl}};
		my $sum =0; foreach (@cnts){$sum += $_;}
		my $sumS =0; for (my $i=0;$i<@cnts; $i++){if ($keys[$i] !~ m/\?$/){$sumS += $cnts[$i] ;}}
		#print "$sumS/$sum\n";
		next if (($sumS/$sum) < 0.7); #at least 60% should be defined at species level..
		my $hiHit =0; my $hiVal="";
		for (my $i=0;$i<@cnts; $i++){
			#print "$keys[$i]\n";
			if ($cnts[$i] > $hiHit){$hiHit = $cnts[$i]; $hiVal = $keys[$i];}
		}
		#print "$hiHit/$sumS $hiVal \n";
		my $replRate = ($hiHit/$sumS);
		if ($replRate > 0.75 && $hiVal !~ m/\?$/){
			print int($replRate*100)."%: Replacing $MGStax{$MGSl} with $hiVal\n";
			$MGStax{$MGSl} = $hiVal;
			$replTax ++;
		}
	}
}
print "\n\nReplaced $replTax taxa to species level based on LCA assignments\n";

my %Tax2MGS;  my $MGSatSpecLvl=0;
foreach my $MGS (keys %MGStax){
	$Tax2MGS{$MGStax{$MGS}} = $MGS;
	$MGSatSpecLvl++ if ($MGStax{$MGS} !~ m/\?$/);
}

my $dirtGene = 0; my $TaxTakenGene=0;my $LCAused=0;
#my %newGTtax;
open O,">$outFile2" or die $!; #this file will be used for guiding abundance estimates..
open O2,">$outFile2.nonGTDBtk" or die $!;
print O "geneID\t$LCAhead;MGS\n";
print O2 "geneID\t$LCAhead;MGS\n";
foreach my $mark (@mgkeys){
	$wrLCA++;
	if (exists($gene2MGS{$mark})){
		my @MGSs = @{$gene2MGS{$mark}};
		$ambGenes ++ if (@MGSs > 1);
		foreach my $MGSl (@MGSs){
			die "Couldn't find MGS tax for $MGSl!\n" unless (exists($MGStax{$MGSl}));
			print O $mark."\t".$MGStax{$MGSl}.";$MGSl\n";
			print O2 $mark."\t".$LCA{$mark}.";$MGSl\n";
		}
		$link2MGS++;
		next;
	} elsif (exists($geneDirtyMGS{$mark})){
		$dirtGene++;
		next;
	} elsif (exists($Tax2MGS{$LCA{$mark}})){
		$TaxTakenGene++;
		next;
	}
	print O2 $mark."\t".$LCA{$mark}.";?\n";
	print O $mark."\t".$LCA{$mark}.";?\n";
	$LCAused++;
}
close O;
close O2;


#write out new GTDB tax file with corrected, tab delim tax
open O,">$outDir/MGS.GTDB.LCA.tax" or die $!; 
foreach my $MGS (sort(keys %MGStax)){
	my @spl = split /;/,$MGStax{$MGS};
	print O "$MGS\t". join("\t",@spl) . "\n";
}
close O;



print "\nWrote $link2MGS / $wrLCA (" . int($link2MGS / $wrLCA * 1000+0.5)/10 . "%) LCA keys with MGS link, $ambGenes ambigous, $dirtGene fuzzy.\n$MGSatSpecLvl/". scalar(keys%MGStax) . " MGS at species level. \nFound $LCAentries genes total in $LCAfiles .LCA files, $TaxTakenGene skipped due to conflicts with MGS tax, used $LCAused LCA assigned genes.\n$outFile2\n";


print "Starting rtk tax sumup..\n";

my $rarBin = getProgPaths("rare");#

my $cmd = "";
#highFromLow: take higher lvl abundance by summing lowest level abundance (that are based on mean)
$cmd .= "$rarBin sumMat -i $GCdir/Matrix.mat.gz -o $matrOutDir/MGS.mat -t $numCor -refD $outFile2  -filterOccPerSmpl 8 -filterMinOcc 15 -extHiera -funcHieraSep \";\" -mean -highFromLow \n";# -mean -filterOccPerSmpl 5 -filterMinOcc 5 -extHiera
print $cmd;
systemW $cmd;
#die "$cmd\n";

