#!/usr/bin/perl
#script that takes a selection of MGS genes (canopy format) and sorts them based on a) marker genes b) copy unmber c) overall occurence
#./resortMGSgenes4importance.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/ /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can.Rhcl.mgs
use warnings;
use strict;
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(gzipopen readFasta writeFasta systemW);
use Mods::TamocFunc qw(readTabbed);
use Mods::math qw(meanArray medianArray quantileArray);

sub evalCurMGS;


#v0.1: adopt .core MGS files to get additional info for sorting genes by importance
#v0.11: 9.2.24: retain more genes/MGS
#v0.12: 11.2.24: adopted to weighted multiBin scores; more subs to make script more modifiable
my $version = 0.12;


# set up some base variables
my $rareBin = getProgPaths("rare");
my $GCd = $ARGV[0];
my $MGSfile = $ARGV[1];
my $useGTDBmg = $ARGV[2];
my $obsFile = $MGSfile; #$ARGV[3];# if (@ARGV > 3);
$obsFile =~ s/\.core$//; $obsFile.=".obs";
die "ARG 2 option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");

#main output file
my $finout = "$MGSfile.srt";
if (-e $finout && -s $finout){
	#print $finout."\n"; exit(0);
	print "Overwriting $finout\n";
}


print "\n--------------------------------------------------\nResorting MGS genes for importance in strain phylo ver $version\n--------------------------------------------------\n";
#my @FMG40 = ("COG0012","COG0016","COG0018","COG0048","COG0049","COG0052","COG0080","COG0081","COG0085","COG0087","COG0088","COG0090","COG0091","COG0092","COG0093","COG0094","COG0096","COG0097","COG0098","COG0099","COG0100","COG0102","COG0103","COG0124","COG0172","COG0184","COG0185","COG0186","COG0197","COG0200","COG0201","COG0202","COG0215","COG0256","COG0495","COG0522","COG0525","COG0533","COG0541","COG0552");
#my %FMG40 = map { $_ => 1 } @FMG40;


die "Can't find main infile $MGSfile\n" unless (-s $MGSfile);

#read MGS occurrence to understand distribution
my $hr = readTabbed($obsFile);
my %MGSobs = %{$hr};

#load GTDB/FMG genes directly..
my $inMGFile="$GCd/FMG.subset.cats";
if ($useGTDBmg eq "GTDB"){
	$inMGFile="$GCd/GTDBmg.subset.cats";
}

#load list of reference marker genes (to mark these later as important genes)
my %MGset=();
open I,"<$inMGFile" or die "resortMGSgenes4importance.pl: Couldn't open $inMGFile\n"; 
my $totMGSgenes=0;
while (<I>){
	my @spl1 = split /\t/;
	my @spl2 = split /,/,$spl1[2];
	foreach my $gene (@spl2){$MGset{$gene} = 1;$totMGSgenes++;}
}
close I;
print STDERR "Loaded $totMGSgenes $useGTDBmg marker genes from \n";



#alt: go with compl.incompl.95.fna.clstr.idx to calc gene occurrences..
my %geneOcc;
my ($I,$ST) = gzipopen("$GCd/compl.incompl.95.fna.clstr.idx","gene cat index file");
#my $maxOcc = 0;
while (<$I>){
	chomp; my @spl= split /\t/;
	if (@spl<2){next; }#die $1." no tab char\n";}
	my $count = $spl[1] =~ tr/,//;
	$geneOcc{$spl[0]} = $count;
	#$maxOcc = $spl[1] if ($maxOcc < $spl[1]);
}
close $I;
print STDERR "Counted gene occurrence for " . scalar(keys(%geneOcc)) . " genes\n";

#my $gene2taxF = createGene2MGS($MGSfile,$GCd);
#my %gen2Bin ; my %gene2COG;
#open I,"<$gene2taxF" or die "Cant' open $gene2taxF\n";
#while (<I>){
#	chomp; my @spl = split /\t/;
#	push (@{$gen2Bin{$spl[1]}}, $spl[0]);
#	$gene2COG{$spl[0]} = $spl[2] if (defined ($spl[2]));
#}
#close I;


#my $sortedOutFile = "$MGSfile.srt";
open O,">$finout" or die "can't open outfile $finout\n";
open I,"<$MGSfile" or die "cant open infil $MGSfile\n";
my $cn=0; my $MGScnt = 0; my $geneCnt=0;
my $curMGS=""; 
my %occ; my %multiCp; my %markers; my %multiBin; 
#foreach my $mg (keys %gen2Bin){
while (my $line = <I>){
	chomp $line;
	next if ($line =~ m/^#/); #commented line
	my @spl = split /\t/,$line;
	my $MGS = $spl[0];  
	$curMGS = $MGS if ($curMGS eq ""); 
	if ($MGS eq ""){die "Undefine MGS on line $line\n";next;}

	if ($MGS ne $curMGS){
		my $retS = evalCurMGS($MGS);
		print O $retS;
	}
	my $gene = $spl[1];
	#push @genes,$gene; 
	next if (exists($markers{$gene}));
	if ($spl[5] == 1){
		$markers{$gene}=$spl[2];
	} else {
		die "$gene in both markers and not!\n" if (exists($markers{$gene}));
		$occ{$gene} = $spl[2]; 
	}		
	$multiCp{$gene} = $spl[3]; $multiBin{$gene} = $spl[4];
	$cn ++;
	

}
close O;
close I;
#report that all went fine
print  "Finished \nProcessed $cn genes, used $geneCnt genes in $MGScnt MGS\nSaved in $finout\n";


exit(0);



sub evalCurMGS{
	#this routine decides which MGS genes (already pre-filtered for core genes) will be handed on to strain phylo construction.. should be "certain" cutoffs for removing genes (intra phylo will do another round of filtering)
	my ($MGS) = @_;
	die "Can't find observed val for MGS $MGSobs{$curMGS}\n" unless (exists($MGSobs{$curMGS}));
	my %finalList; my $mrkCnt=0; my $avgOcc=0;
	my $maxMocc=0; my $minMocc=10000000;
	my $MGSob = $MGSobs{$curMGS}; #my $MGSob001 = ((0.01*$MGSob)+1 );
	my @mrks = keys %markers;
	my %mBinMrks = %multiBin{@mrks};
	#my $medMCp = medianArray(values %multiCp);	my $avgMCp = meanArray([values %multiCp]);
	my $medMBi = medianArray(values %multiBin);	my $avgMBi = meanArray([values %multiBin]);
	my $q75MBi =quantileArray(0.75,values %multiBin); 
	my $q75MBiM =quantileArray(0.75,values %mBinMrks);
	my $q75MBiMF = $q75MBiM;
	if ($q75MBiMF > 2.2){$q75MBiMF=2.2;print "Warning: very high markerG multibin: $q75MBiM\n";}
	foreach my $gn (@mrks){
		if ($multiBin{$gn} > $q75MBiMF #|| $multiCp{$gn} > $MGSob001 
		){
			#print "$gn :: $multiBin{$gn} > 1 || $multiCp{$gn} > 0\n";
			next;
		} 
		$finalList{$gn} =  scalar(keys %finalList);
		$avgOcc += $geneOcc{$gn};
		$mrkCnt++;
		$maxMocc = $markers{$gn} if ($markers{$gn} > $maxMocc);
		$minMocc = $markers{$gn} if ($markers{$gn} < $minMocc);
	}
	$avgOcc /= $mrkCnt; 
	my $avgOcc2 = $avgOcc; $avgOcc2 = 1 if ($avgOcc < 1);
	my @srtedGenes = sort {$occ{$b} <=> $occ{$a}} keys %occ;
	#make sure lower and upper bound is not unreasonable..
	my $avgOccH = $avgOcc2 *1.5; if ($avgOccH < 4){$avgOccH = 4;} 
	my $avgOccL = int($avgOcc2 *0.5); #if ($avgOccL > 4){$avgOccL = 4;} 
	foreach my $gn (@srtedGenes){
		if ($multiBin{$gn} > ($q75MBi*1.1) #|| $multiCp{$gn} > $MGSob001
				|| $geneOcc{$gn} > ($avgOccH) || $geneOcc{$gn} < $avgOccL
				){
			next;
		}
		$finalList{$gn} = scalar(keys %finalList);
	}
	
	
	
	print "${curMGS} (".scalar(keys %multiBin)."):: " ;
	print scalar(keys %finalList) ." genes, $mrkCnt markerGs used, avgOcc: " . int($avgOcc*100)/100 . ", MAG occ: $MGSob, median/avg/q75/q75Mark multiBin $medMBi/" . int($avgMBi*100)/100  . "/". int($q75MBi*100)/100  ."/". int($q75MBiM*100)/100  ."\n";
	
	$geneCnt+= scalar(keys %finalList);
	
	my @finalGs = sort { $finalList{$a} <=> $finalList{$b} } keys(%finalList);
	#my @finalGs = @finalList{@keys};
	#die "\n@finalGs\n";

	my $retStr= $curMGS."\t".join(",", @finalGs)."\n";
	
	
	$curMGS = $MGS; 
	%occ = (); %multiCp= (); %markers= (); %multiBin= ();
	$MGScnt++;
	return $retStr;
}











