#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long qw( GetOptions );
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(readMap readClstrRev readClstrRevContigSubset readClstrRevSmplCtgGenSubset getDirsPerAssmblGrp  getAssemblPath systemW gzipopen parse_duration);
use Mods::Binning qw (runMetaBat runCheckM runCheckM2 createBinFAA readMGS filterMGS_CM MB2assigns minQualFilter calcLCAcompl readCMquals);
use Mods::geneCat qw(readMG_LCA);


sub countUpBin;
sub countStats;
sub gene2GCg;
sub MGuniqStats;

#.12: LCA integration
#.13: .cm2 write (rep MAG)
#.14: integrate canopies as well..
#.15: tightened substantially cluster addition, reduced num cycles
#.16: extended decluter step, using binocc as primary guide
#.17: complete rewrite, MAG clustering now prioritizes low contamination and is much faster
#.18: switched to reading sample_contig_gene from idx
#.20: switched to C++ implementation of clusterMAGs.. this is only a hull script now, mostly disused.

my $version = 0.20;

my $startTime = time ;

my $inD="";my $binSpeciesMG = 0; 
my $outD ="";
my $useGTDBmg = "FMG";
my $useCheckM2 = 1;
my $useCheckM1 = 0;
my $tmpD = ""; my $numCore = 1;
my $ignoIncomplMAGs = 1; my $logDir = "";
my $redo = 0;
my $BinTerm = "MGS.";
my $legacyV=0;
my $camoIn = "";

my $ph1flag = 1; #sets up for using binnings..
my %gen2Bin;#structure: {gene}{Bin}=cnt
my %genesInMAG; #struc: {gene} = Bin
my %reprMAGpMGS; #1 MAG that best represents MGS..
my %binName; #name BinNum->MGS given to Bin
my %binName2; #name MAGname->MGS


#MAG specific qual, genes etc, using uniqBinID;
my %MAGc; my %MAGq; #overall storage of qual and conta
my %MAGlcaq; #overall LCA completeness
my %MAgene; #store genes per MAG.. large hash!



#-GCd $inD -binSpeciesMG $binSpeciesMG -logDir $logDir -MGset $useGTDBmg -cores $numCore -useCheckM1 $useCheckM1 -useCheckM2 $useCheckM2 
#options to pipeline..
GetOptions(
	"GCd=s"      => \$inD,
	"BinDir=s"   => \$outD,
	"tmp=s" => \$tmpD,
	#"submit=i" => \$doSubmit,
	#"canopies=s" => \$canopyF,
	"cores=i" => \$numCore,
	"MGset=s" => \$useGTDBmg,#GTDB or FMG
	#"mem=i" => \$memG,
	#"strains=i" => \$doStrains,
	"logDir=s"  => \$logDir,
	"canopies=s" => \$camoIn,
	"useCheckM2=i" => \$useCheckM2,
	"useCheckM1=i" => \$useCheckM1,
	"binSpeciesMG=i" => \$binSpeciesMG, #0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder
	"ignoreIncompleteMAGs=i" => \$ignoIncomplMAGs,
	"redo=i" => \$redo,
	"legacy=i" => \$legacyV,
);

print "#######################################################\nclusterMAGs algorithm, v$version\n#######################################################\n\n";
print "Running in legacy mode\n" if ($legacyV);

my $cmSuffix = ".cm"; $cmSuffix = ".cm2" if ($useCheckM2); 

my $BinnerShrt = "MB2";
if ($binSpeciesMG == 2){$BinnerShrt = "SB";}#SemiBin
if ($binSpeciesMG == 3){$BinnerShrt = "MD";}
if ($outD eq ""){$outD = $inD."/Bin_$BinnerShrt/";}
$outD .= "/" unless ($outD =~ m/\/$/);
my $ctg2gen = {};
my $ctg2gen2 = {};

#dir in GCdir where marker genes are stored..
my $COGdir = "FMG";
if ($useGTDBmg eq "GTDB"){ 	$COGdir = "GTDBmg";}

system "rm $outD/$BinnerShrt.clusters"  if ($redo);


#read map to get assembly groups..
my $mapF="";my $GCd = "";
#die Dumper($hrm);	
if (-e "$inD/LOGandSUB/GCmaps.inf"){
	my $tmp = `cat $inD/LOGandSUB/GCmaps.inf`; chomp $tmp;
	$mapF = $tmp;
	$GCd = $inD;
} else{
	die "can't find indir $inD\n";
	$mapF = $inD."LOGandSUB/inmap.txt";
	#($hrm,$asGrpObj) = readMap($inD."LOGandSUB/inmap.txt");
}
die "Couldn't find map file in clusterMAGs.pl:: $mapF\n" unless (-e $mapF|| $mapF =~ m/,/);

my $clMAGsBin = getProgPaths("clusterMAGs");

my $cmd = "";
my $canoFlag = ""; $canoFlag = "-canopyDir $camoIn " if ($camoIn ne "");
#-FILEtag SBx -MGtag MM2 -geneCatIdx C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/compl.incompl.95.fna.clstr.idx -MGdir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/MGs/ -outDir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/out/ -map C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/map.0.txt,C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/map.1.txt -canopyDir C:\Users\hildebra\OneDrive\science\data\test\clusterMAGsMock/Cano/
$cmd .= "$clMAGsBin  -CMsuffix .cm2 -FILEtag $BinnerShrt -MGStag MGS. -geneCatIdx $GCd/compl.incompl.95.fna.clstr.idx -LCAdir $GCd/${COGdir} ";
$cmd .= "-outDir $outD -map $mapF $canoFlag -MGfile $GCd/GTDBmg.subset.cats ;\n";
$cmd .= "gzip -c $outD/MAGvsGC.txt > $logDir/MAGvsGC.txt.gz;\nrm $outD/MAGvsGC.txt;\n";

if (1){ #C++ path.. better unless for testing something
	print "$cmd\n";
	system $cmd;
	print "\n\nDone with binary-based MAG clustering.\n\n";
	exit;
}



my ($hrD,$hrM) = getDirsPerAssmblGrp($mapF);
my %map = %{$hrM};
my %DOs = %{$hrD};
my @DoosD = sort keys %DOs; #dirs of assembly groups
print "Found ".scalar(@DoosD) ." assembly groups\n";
system "mkdir -p $logDir\n";



my %gene2COG;  my %uniCOGs;
#my %FMGlist;
print "Reading COG to MGs\n";
open IC,"<$GCd/${COGdir}.subset.cats" or die "Can't open $GCd/${COGdir}.subset.cats\n";
while (<IC>){
	chomp;	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my @genes = split(/,/,$spl[2]);
	my $curCOG = $spl[0];
	#$gen2COG{$curCOG} = \@genes;
	$uniCOGs{$curCOG} = 1;
	foreach (@genes){$gene2COG{$_} = $curCOG;}
}
close IC;

my @uniqCOGs = keys %uniCOGs;
my $LCAhref = readMG_LCA("$GCd/${COGdir}",\@uniqCOGs,6);
#my %LCA = %{$LCAhref};

my %Canos ; my %CanQual;
if ($camoIn ne ""){
	#my ($canosHR,$canQualHR) = MB2assigns($camoIn,$camoIn.$cmSuffix); 
	my $canosHR = readMGS($camoIn);
	my $canQualHR = readCMquals($camoIn.$cmSuffix);
	%Canos = %{$canosHR}; %CanQual = %{$canQualHR};
	print "Read ". scalar(keys(%Canos)) ." canopy MGS\n";
}

clusterMB2();

print "Primary MAG clustering done \n";

print "Elapsed time: ". parse_duration((time - $startTime)) . "\n";

exit(0);










sub decluter($){
	#my ($g2bhr,$bghr,$BOhr) = @_;
	my ($BOhr) = @_;
	print "Further declutering gene->MGS assignments to increase uniquely assigned genes\n";
	#my %gen2Bin = %{$g2bhr};
	#my %BinGrps = %{$bghr};
	#my %BinObs = %{$BOhr};
	my @genes = keys %gen2Bin;
	
	my $report = "";
	#assign genes to single bins based on where they occur most often..
	my %GeneMultiBin;	 my %bin2gene;  #my %cntsPerBin;
	my $clearGene=0; my $amb1Gene =0 ; my $multAss=0; my $unAss=0;
	my $hiFracUniqAssi = 0;my $hiFracMultiAssi = 0;my $hiFracFailAssi = 0;
	foreach my $gen (@genes){
		my @bins = keys %{$gen2Bin{$gen}};
		if (@bins == 1 ){ #simple, unambigous case...
			$bin2gene{$bins[0]}{$gen} = $gen2Bin{$gen}{$bins[0]};
			$clearGene++;
			#$GeneMultiBin{$gen} = 1;
			next;
		}
		#stats on gene..
		my $sum=0; my $hiCnt=0; my $hiB;
		my $hiFrac = 0; my @hiFB;
		foreach my $b (@bins){
			#bin where this gene has highest count in.. will draw most genes to high prevalence bins
			my $ccnt = $gen2Bin{$gen}{$b};
			$sum += $ccnt;
			if ($ccnt > $hiCnt){$hiCnt = $ccnt;	$hiB = $b;}
			#instead calc expected fraction for each gene: which Bin does it fit best to, based on prevalence?
			my $bFrac = ($ccnt / $BOhr->{$b});	if ($bFrac > 1){$bFrac -= (1-$bFrac);}
			if ($bFrac>0.7){
				$hiFrac = $bFrac; push(@hiFB , $b);
			}
		}
		#stats on effectiveness
		my $numFB = @hiFB;
		if ($numFB == 0){ $hiFracFailAssi ++;$GeneMultiBin{$gen} = @bins;
		} elsif ($numFB == 1) {$hiFracUniqAssi ++;$GeneMultiBin{$gen} = $numFB;
		} else {$hiFracMultiAssi++;$GeneMultiBin{$gen} = $numFB;}
		
		
		
		
		foreach my $fBin (@hiFB){
			$bin2gene{$fBin}{$gen} = $gen2Bin{$gen}{$fBin};
		}
		
		#outdated..
		next;
		my $hiRat = 0; $hiRat = ($hiCnt/$sum) if ($sum>0);
		$GeneMultiBin{$gen}=$hiRat;
		if ( ($hiRat) > 0.9){
			$clearGene++;
			#$cntsPerBin{$hiB} += $hiCnt;
			#$bin2gene{$hiB}{$gen} = $hiCnt;
			#$gen2Bin{$gen}={};	$gen2Bin{$gen}{$hiB} = $hiCnt;
			$bin2gene{$hiB}{$gen} = $hiCnt;
			
		} elsif ( ($hiRat) > 0.6 ){
			$amb1Gene++;
			#$cntsPerBin{$hiB} += $hiCnt if ( $hiB > 30);
			#$gen2Bin{$gen}={};			$gen2Bin{$gen}{$hiB} = $hiCnt;
			$bin2gene{$hiB}{$gen} = $hiCnt;
		} elsif ( ($hiRat) > 0.3) { #just assign everything above 0.3
			#my %newHash;
			foreach my $b (@bins){
				if ( ($gen2Bin{$gen}{$b} / $sum ) > 0.1 ){
					#$cntsPerBin{$b} += $gen2Bin{$gen}{$b};
					$bin2gene{$gen}{$b} = $gen2Bin{$gen}{$b};
				} #else {delete $gen2Bin{$gen}{$b};}
			} 
			#$gen2Bin{$gen}={};
			$multAss++;
		} else { #seems to be all over hte place gene.. don't use
			$unAss++;
		}
		
	}
	$report .=  "Genes assigned to bins:\n";
	#$report .=  "Unique assigned(>90%): $clearGene, ambivalent(>60%): $amb1Gene multiBin(>30%): $multAss unclear(<30%): $unAss \n" ;
	$report .=  "Unique assigned: $clearGene, ambivalent (>70% uniq): $hiFracUniqAssi multiBin(>70% multi): $hiFracMultiAssi unclear(<70%): $hiFracFailAssi \n" ;

	

	
	#$report .= MGuniqStats();
	return (\%GeneMultiBin, \%bin2gene, $report);
}



sub writeGene2MGS{
	my ($ofile,$ofile2,$mainStatStr,$BinObsHR,$BinGrpsHR,
		$bin2geneHR,$GeneMultCopyHR,$GeneMultiBinHR,$binQualTierHR) = @_;
	print "Writing Gene->MGS links ($ofile)\n";
	#write out genes
	#now create the names for each Bin (MGS.xx)
	my @BinOrd = sort {$BinObsHR->{$b} <=> $BinObsHR->{$a}|| $b cmp $a	} keys %{$BinObsHR};
	my $bcnt=1; 
	foreach my $bi (@BinOrd){
		my $bN = "$BinTerm$bcnt";
		$binName{$bi} = $bN;
		$bcnt++;
		foreach my $MAGn (@{$BinGrpsHR->{$bi}}){
			$binName2{$MAGn} = $bN;
		}
	}
	
	my %binMem;
	open O,">$ofile" or die "Can't open $ofile\n";
	print O "Bin\tGene\tOcc\tMultiCopy\tMultiBin\tisMarkerGene\n";
	#sort by both abundance and bin name
	foreach my $bi (@BinOrd){
		my $bN = $binName{$bi} ;
		my $laterStr="";
		$binMem{$bi} = join ",",@{$BinGrpsHR->{$bi}};
		my @srtdGenes = sort {$bin2geneHR->{$bi}{$b} <=> $bin2geneHR->{$bi}{$a}} keys %{$bin2geneHR->{$bi}};
		foreach my $g (@srtdGenes){
			my $isMG=0;		if (exists($gene2COG{$g})){$isMG=1};
			my $str= "$bN\t$g\t$bin2geneHR->{$bi}{$g}\t";
			if (exists($GeneMultCopyHR->{$g})){$str .= "$GeneMultCopyHR->{$g}\t";} else { $str .= "0\t"; }
			if (exists($GeneMultiBinHR->{$g})){$str .= "$GeneMultiBinHR->{$g}";} else {$str .= "1"; }
			$str .= "\t$isMG\n";
			if ($isMG){
				print O $str;
			} else {
				$laterStr .= $str;
			}
		}
		print O $laterStr;
	}
	close O;


	#
	print "Writing MGS occurrences ($ofile2)\n";
	if (@BinOrd > 1){
		open O,">$ofile2" or die "Can't open out $ofile2\n";;
		print O "Bin\tObservations\tQualTier\tMembers\n";
		foreach my $bi (@BinOrd){
			print O "$binName{$bi}\t$BinObsHR->{$bi}\t$binQualTierHR->{$bi}\t$binMem{$bi}\n";
		}
		close O;
		open O,">$logDir/mergMAG.log" or die "mergeMag can't open $logDir/mergMAG.log";
		print O $mainStatStr;
		close O;
	}
}


sub clusterMB2{ 
	if (-e "$outD/$BinnerShrt.clusters"){print "$outD/$BinnerShrt.clusters already exists.. nothing to do here\n"; return; }
	#converges MB2 based on quality of bins and gene catalog..
	print "Reading gene index reverse\n";
	my $hr2;
	#($ctg2gen,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",2); $hr2 = {};
	#%ctg2gen = %{$hr1};
	#die;

	
	
	#just use one log file..
	#open LOG,">$outD/MB2.log";
	#die $ctg2gen{">MM1__C12020_L=3295=_1"}."\n";
	print "Cluster Bins\n\n";
	#die "DONE\n";

	my $BinNumDeNovo = 0; my $scnt=0;
	my $binQualHR = {};
	my %ambBin ; 
	my $novelBin=0; my $matchBin=0;my $unknwnBin=0;
	my $CanoAdd=0; my $CanoNovel=0;
	#too lenient
	#	my @ComplTiers = (95,90,80,90,60,50); my @ContaTiers = (5,5,5,10,5,10); my @AllowNovelMGS = (1,1,1,1,1,0);
	#stricter so I have less to do later
	my @ComplTiers = (98,95,80,60); my @ContaTiers = (1,3,5,10); 
	 my @LCAcompl = (0.97, 0.95, 0.92, 0.9);my @AllowNovelMGS = (1,1,1,0); 
	 my $maxRounds = 3; my @OverlapTiers = (0.8,0.8,0.9,0.9); my @unknwnTiers = (0.75,0.75,0.8,1.1);
	if ($legacyV){
		@ComplTiers = (95,90,80);  @ContaTiers = (5,5,5);  @AllowNovelMGS = (1,1,1); @LCAcompl = (0,0,0);
		$maxRounds=3;
	}
	my @createPerTier = (0,0,0,0,0);
	
	my %regBins;
	my %BinObs; #saves how often bin is observed..
	my %GeneMultCopy; #saves how often a gene occurrs in bins
	my %binQualTier;
	my %BinGrps;
	my %MAG2Bin; my %Bin2MAG;
	
	my @DirSet = sort keys %DOs;
	push(@DirSet,"__CanopyMGSset__") if (scalar(keys(%Canos)));
	
	
	#DEBUG
	#$DirSet[1] = "__CanopyMGSset__";
	#1: go over all dirs, process MAGs, store in RAM
	my %BinsGlobal; my %contigsGlobal;#to be deleted later..
	foreach my $Doo (@DirSet){
		my $hr1;my @smplIDs;my $uniqMBid ="XX";my $isCanopy=0;
		if ($Doo eq "__CanopyMGSset__"){ #get canopies.. only activate at the end of a run
			$hr1 = \%Canos; $hr2 = \%CanQual;
			@smplIDs = ("Cano"); $isCanopy = 1;
			print "\nAdding MGS from Canopies (N=". scalar(keys(%Canos)) . ").";
			#next;
		} else {
			@smplIDs = @{$DOs{$Doo}{SmplID}}; my @paths = @{$DOs{$Doo}{wrdir}};
			my $metaGD = getAssemblPath($paths[-1]);
			my $MBf = $metaGD."/Binning/$BinnerShrt/$smplIDs[-1]";	my $MBfQual = $MBf.$cmSuffix;
			next unless (-e $MBfQual);
			print "$MBf\n";
			($hr1,$hr2) = MB2assigns ($MBf,$MBfQual);
		}
		#some stats on bin quality
		foreach my $bin (keys %{$hr1}){
			$binQualHR = countStats($hr2,$bin,$binQualHR) ;
		}
		#only keep high qual bins..
		($hr1,$hr2) = minQualFilter($hr1,$hr2,$ComplTiers[-1], $ContaTiers[-1], $LCAcompl[-1]);
		#%MB = %{$hr1}; %MBQ = %{$hr2};
		#log overall stats for bins..
		foreach my $bin (keys %{$hr1}){
			$uniqMBid = "$smplIDs[-1].$bin";
			$MAGq{$uniqMBid} = ${$hr2}{$bin}{compl};$MAGc{$uniqMBid} = ${$hr2}{$bin}{conta};
			$BinsGlobal{$uniqMBid} = ${$hr1}{$bin};
			foreach my $kk (@{$BinsGlobal{$uniqMBid}}){
				$contigsGlobal{$kk} = 1;
				#print "$kk\n";
			}
			#$MAgene{$uniqMBid} = gene2GCg(${$hr1}{$bin});
			#$MAGlcaq{$uniqMBid} = calcLCAcompl($MAgene{$uniqMBid},$LCAhref);
			#print "LCA: $MAGlcaq{$uniqMBid} ";
		}
	}
#	$ctg2gen = readClstrRevContigSubset("$GCd/compl.incompl.95.fna.clstr.idx",\%contigsGlobal); $hr2 = {};
	print parse_duration((time - $startTime)). " - ";
	$ctg2gen2 = readClstrRevSmplCtgGenSubset("$GCd/compl.incompl.95.fna.clstr.idx",\%contigsGlobal); $hr2 = {};

	print parse_duration((time - $startTime)) . " - Converting Bins to gene cat MAGs..\n";
	my $missed =0; my $incompleteCtgMatch=0; 
	my %Tiers; my @TierCnt; 
	for (my $T=0;$T<$maxRounds;$T++){$Tiers{$T} = [];$TierCnt[$T] = 0;}
	foreach my $uBin (keys %BinsGlobal){
		($MAgene{$uBin},$missed) = gene2GCg2($BinsGlobal{$uBin});
		$MAGlcaq{$uBin} = calcLCAcompl($MAgene{$uBin},$LCAhref);
		if ($missed){ print "$uBin ";$incompleteCtgMatch ++;}
		#now I have all info.. add MAGs to tiers..
		for (my $T=0;$T<$maxRounds;$T++){ my $conti=1;
			if ( $MAGq{$uBin} >= $ComplTiers[$T] && $MAGc{$uBin} <= $ContaTiers[$T] 
				&& scalar(@{$MAgene{$uBin}}) != 0 && $MAGlcaq{$uBin} >= $LCAcompl[$T]){
					push(@{$Tiers{$T}}, $uBin);
					$TierCnt[$T] ++;
					$conti=0;
				}
				last if ($conti == 0);
		}
	}
	print parse_duration((time - $startTime)) . " - Finished converting Bins to Gene Cat genes\n";
	print "Bins in Qual Tiers: @TierCnt \n";
	print "$incompleteCtgMatch/" . scalar(keys %BinsGlobal) . " MAGs with incomplete contig assignments\n" if ($incompleteCtgMatch);
	#delete vars not required any longer..
	%BinsGlobal = (); %contigsGlobal = ();
	$ctg2gen2 = {};
			
	
	#primary data collection (MAGs + quals + gene cat translation) done..
	#next: sort by qual
	print "Loaded ". scalar(keys %MAgene) . "/" . scalar(keys %MAGc) . " unique high-qual MAGs. Sorting.. \n";
	
	my @order ; my %loggedBins; my $ocnt = 0; my @tierSw;
	#my @orderPre = sort { $MAGc{$a} <=> $MAGc{$b} } keys(%MAGc);
	for (my $T=0;$T<$maxRounds;$T++){
		#foreach my $bin (@orderPre ){
		#	if (exists($loggedBins{$bin}) || $MAGq{$bin}< $ComplTiers[$T] || $MAGc{$bin} > $ContaTiers[$T] 
		#		|| scalar(@{$MAgene{$bin}}) == 0 || $MAGlcaq{$bin} < $LCAcompl[$T]){next;}
		#	$loggedBins{$bin} = 1;push (@order,$bin); $ocnt++;
		#}
		next unless (@{$Tiers{$T}} > 0);
		my %subH = %MAGc{@{$Tiers{$T}}};
		#print " subhash size: " . scalar(keys %subH); #DEBUG
		push(@order, sort { $subH{$a} <=> $subH{$b} } keys(%subH));
		push (@tierSw,$order[-1]);
		print parse_duration((time - $startTime)) .  " - Sorted T$T " . scalar(keys %subH) . " MAGs\n";
	}
	%loggedBins = (); %Tiers = ();#@orderPre = ();
	
	print parse_duration((time - $startTime)) . " - Sorted $ocnt MAGs.. starting clustering\n";
	
	#sorted MAGs by qual. Next: cluster MAGs based on marker gene overlaps
	
	my $locT = 0; #needs to increase with bins..
	foreach my $bin (@order){
		my $isCanopy = 0;
		my ($BinHR,$BinMGHR,$ar,$arMG, $totGenes,$assignedGenes,$nonMGgene) = countUpBin($MAgene{$bin},1);
		my %BinCnt = %{$BinHR};my %BinMGcnt = %{$BinMGHR};
		#list of marker genes in Bin
		my @MGsG = @{$arMG};
		my $MGgenes =  @MGsG ;
		my $BinNumL = $BinNumDeNovo; my $HiCnt=0; my $HiBin= "-1";
		
		foreach my $bin (keys %BinCnt){ 
			next if ($bin eq "-1");
			if ($BinCnt{$bin} > $HiCnt){$HiCnt = $BinCnt{$bin}; $HiBin = $bin;}
		}
		my $unknwnFrac = 0;	if (exists $BinCnt{"-1"}){ $unknwnFrac = ($BinCnt{"-1"}/$nonMGgene) ;}
		my $HiFrac = $HiCnt / $nonMGgene;#best match to other MG
		#same numbers, but focus on marker genes (more important measure)
		my $HiMGCnt=0;my $HiBinMG= "-1";
		foreach my $bin (keys %BinMGcnt){ 
			next if ($bin eq "-1");
			if ($BinMGcnt{$bin} > $HiMGCnt){$HiMGCnt = $BinMGcnt{$bin}; $HiBinMG = $bin;}
		}
		my $MGunknwnFrac = 0;	if (exists $BinMGcnt{"-1"}){ $MGunknwnFrac = ($BinMGcnt{"-1"}/$MGgenes) ;}
		my $MGHiFrac = $HiMGCnt / $MGgenes;#best match to other MG
		my $MGHiFrac2= $MGHiFrac; my $MGHiFrac3 = $MGHiFrac;
		
		if ($MGHiFrac > 0.6){#test if these are "core core" genes being hit
			my $sumHit=0;my $sumHitOcc=0; my $Nhit=0; my $BinObs = $BinObs{$HiBinMG};
			foreach my $gen (@MGsG){
				my $locSc=0;
				if (exists($gen2Bin{$gen})){#just catalogue how often these core MGs are occurring..
					foreach my $sc (values%{$gen2Bin{$gen}}){ $locSc +=$sc;}
					if (exists($gen2Bin{$gen}{$HiBinMG})){
						my $HiBinCnt = $gen2Bin{$gen}{$HiBinMG};
						$sumHit += ($HiBinCnt/$locSc); 
						$sumHitOcc += $HiBinCnt;
						$Nhit++; 
					} 
				}
			}
			$MGHiFrac2 = ($sumHitOcc/($BinObs*$MGgenes)); #occurrence in core
			$MGHiFrac3 = ($sumHit/$MGgenes); #not used by other MGS
			#print "" . ($sumHit/$Nhit) . " $MGHiFrac $BinObs{$HiBinMG} 2:$MGHiFrac2 3:$MGHiFrac3\n$sumHitOcc/($BinObs*$MGgenes)\n$sumHit/$MGgenes\n";
		}
		
		#print "$uniqMBid $HiMGCnt $MGHiFrac $MGHiFrac2 H: $HiFrac unkn: $MGunknwnFrac $unknwnFrac\n";
		if ( $MGunknwnFrac > $unknwnTiers[$locT]){##too many unknown, make a new bin
			if ($AllowNovelMGS[$locT]==0){$unknwnBin++; print " U ";next;}
			$createPerTier[$locT] ++;$binQualTier{$BinNumL} = $locT;
			$novelBin++; $CanoNovel++ if ($isCanopy);
			$regBins{$bin}=1;$BinNumDeNovo ++ ; 
			print "\nCreating from $bin new bin $BinNumL (". int(100*$MGHiFrac) . "%, ". int(100*$MGunknwnFrac) . "%) (Q:$MAGq{$bin}, C:$MAGc{$bin})";
		} elsif ( $MGHiFrac2 > $OverlapTiers[$locT] && $MGHiFrac3 > $OverlapTiers[$locT] ) { #matching to known bin..
			$BinNumL = $HiBinMG; 
			$matchBin++; $CanoAdd++ if ($isCanopy);
			$regBins{$bin}=1;
			print "\nMatch:  bin $BinNumL (". int(100*$MGHiFrac) . "%, ".  int(100*$MGHiFrac2) . "%, " . int(100*$MGHiFrac3) . "%)";
		} else {
			$ambBin{$bin}=1;print " A ";next;
		}
		push(@{$BinGrps{$BinNumL}}, $bin);$MAG2Bin{$bin} = $BinNumL;
		push(@{$Bin2MAG{$BinNumL}}, $bin);$BinObs{$BinNumL} ++;
		#my $genesInseted=0;my $newGene=0;
		my %genesInBin;
		foreach my $curGene (@{$MAgene{$bin}}){
			$gen2Bin{$curGene}{$BinNumL}++;$genesInBin{$curGene} ++;
		}
		#print "\n$genesInseted : $HiCnt : $newGene : $unknwnFrac\n" if ($unknwnFrac < 0.7);
		foreach my $k (keys %genesInBin){ #save which genes were present in multi copy
			if ($genesInBin{$k} > 1){
				$GeneMultCopy{$k} ++;
			}
		}
		#time to increase tier??
		if ($bin eq $tierSw[$locT]){
			 $locT ++;
			 print "\n". parse_duration((time - $startTime)) . "  ------------------------ At Tier $locT ------------------------\n";
		}
	}
	$scnt++;		
	

	
	
	#die;
	

	
	print "\n\n";
	
	#make sure that ambigous bins that were added after all, are removed from the list now..
	my %ambBin2;
	foreach my $k (keys %ambBin){		next if (exists($regBins{$k}));		$ambBin2{$k} = 1;	}
	
	#----- statistics on the MGS grouping --------
	my $mainStatStr="";
	my @genes = keys %gen2Bin;	
	my %binQual = %{$binQualHR};
	$mainStatStr .=  "\n\n-------- de novo MGS creation --------\n";
	$mainStatStr .=  "Completeness cutoffs:".join(",",@ComplTiers)."\n";
	$mainStatStr .=  "Contamination cutoffs:".join(",",@ContaTiers)."\n";
	$mainStatStr .=  "NovelBins	matchingBins	ambigous(rm)	novelAtLowQual(rm)\n";
	$mainStatStr .=  "$novelBin	$matchBin	" . scalar(keys %ambBin2) ."	$unknwnBin\n";
	$mainStatStr .=  "Bins at Qual Tiers: ". join(" ",@createPerTier) . "\n";
	$mainStatStr .=  "Bins above quality:\n";
	$mainStatStr .=  "Completeness (99,95,90,80,60,0): $binQual{g99} $binQual{g95} $binQual{g90} $binQual{g80} $binQual{g60} $binQual{g}\n";
	$mainStatStr .=  "Contamination (1,2,5,10,20,X): $binQual{con1} $binQual{con2} $binQual{con5} $binQual{con10} $binQual{con20} $binQual{con}\n";
	$mainStatStr .=  "Total dataset Bins / Genes: $BinNumDeNovo / ".scalar @genes . "\n";
	$mainStatStr .=  "Canopies de novo MGS: $CanoNovel; added to MGS: $CanoAdd  (out of " . scalar(keys %Canos) . ")\n" if (scalar(keys(%Canos)));

	
	#print "$mainStatStr\n";
	
	#select representative MAG per MGS
	print "Selecting representative MAG per MGS..\n";
	foreach my $MGS (keys %Bin2MAG){
		my @listMAGs = @{$Bin2MAG{$MGS}}; my $hiScore=0; my $hiMAG="";
		foreach my $uniqMBid (@listMAGs){
			my @CCs = @{$MAgene{$uniqMBid}};
			my $N=0; my $sum=0;
			foreach my $gen (@CCs){
				my $fac=1;
				$fac=2 if (exists($gene2COG{$gen}));#bit higher weight on marker genes..
				$sum+=$gen2Bin{$gen}{$MGS}*$fac;$N++;
			}
			my $sco = $sum/$N / $BinObs{$MGS};# if (exists($BinObs{$MGS}));
			$sco += (.01 * ($MAGq{$uniqMBid} - 2*$MAGc{$uniqMBid}));
			if ($sco > $hiScore){
				$hiScore = $sco; $hiMAG = $uniqMBid;
			}
			#recent bins to MAG with most often used genes..
			#$MAG2Bin{$uniqMBid} = $BinNumL;
			#
			#
		}
		$reprMAGpMGS{$MGS} = $hiMAG;
	}
	
	#MGuniqStats();
	
	$mainStatStr .= MGuniqStats();
	
	my ($GMBhr, $b2ghr, $rep) = decluter(\%BinObs);#\%gen2Bin,\%BinGrps,
	$mainStatStr .= $rep;
	#my %GeneMultiBin = %{$GMBhr};
	#my %bin2gene = %{$b2ghr};
	#my @BinOrd = @{$BOar};
	
	#writes multiple reports..
	writeGene2MGS("$outD/$BinnerShrt.clusters","$outD/$BinnerShrt.clusters.obs",$mainStatStr,\%BinObs,
			\%BinGrps,$b2ghr,\%GeneMultCopy,$GMBhr,\%binQualTier);


	
	#some log stats about the MAGs & Bins..
	my $qualStr = summarizeMAGcontent(\%MAG2Bin);
	open OX,">$outD/$BinnerShrt.clusters$cmSuffix" or die $!;
	print OX $qualStr;
	close OX;

	print $mainStatStr."\n";
	print "-------- MAG merging end --------\n\n";
	
	
	
	#return(\%BinGrps);

}










sub MGuniqStats{
	#determine how redundant MG genes are 
	print "Calculating average uniqueness of genes assigned to MGS (pre- postfilter)\n";
	my %BinStatsMGN; my %BinStatsMGbins; my %BinStatsMGOcc; my %BinStatsMGOccOth;
	foreach my $MG (keys %gene2COG){
		if (exists($gen2Bin{$MG})){
			my $lsc = 0; foreach my $sc (values %{$gen2Bin{$MG}}){$lsc += $sc;}
			my @kGen = keys %{$gen2Bin{$MG}};
			foreach my $BN (@kGen){
				$BinStatsMGOcc{$BN} += $gen2Bin{$MG}{$BN};
				$BinStatsMGOccOth{$BN} += $lsc - $gen2Bin{$MG}{$BN};
				$BinStatsMGN{$BN} ++;
				$BinStatsMGbins{$BN} += scalar(@kGen);
			}
		}
	}
	my @uniqA; my $uniqAvg=0;
	foreach my $BN (keys %BinStatsMGN){
		#my $uniqueness = $BinStatsMGOcc{$BN}/($BinStatsMGOcc{$BN}+$BinStatsMGOccOth{$BN});
		my $uniqueness = $BinStatsMGN{$BN}/($BinStatsMGbins{$BN});
		push(@uniqA,$uniqueness);
		$uniqAvg += $uniqueness;
		#print int($uniqueness * 100) . "% ";
	}
	my $mainRep = "Num Bins: " . scalar(keys %BinStatsMGN) . ". Uniqueness MG assiments Avg: " . int($uniqAvg/scalar(keys %BinStatsMGN)*10000)/100 . "%\n";
	print $mainRep;
	foreach my $va (sort @uniqA){
		print int($va*100)."% ";
	}
	print"\n";
	
	return $mainRep;
}



sub countUpBin{
	my ($ar,$isGCgene) = @_;
	my @CCs = @{$ar};
	my %BinCnt; #genes shared between Bins
	my %BinMGcnt; my @genes; my @MGgenes;
	my $totGenes=0; my $nonMGgene=0; my $MGgenesC=0;
	my $seenGeneCnt=0; my $seenMGcount=0;
	if (@CCs == 0){
		return (\%BinCnt,\%BinMGcnt,\@genes,\@MGgenes,$totGenes,$seenGeneCnt,$nonMGgene,$MGgenesC) ;
	}
	foreach my $cc (@CCs){ #$cc are contigs, therefore needs to test for all genes potentially on contig
		#die "$cc\n";
		my $cnt=0; my $curGene = "";my $miscnt=0;my $testKey="";
		my $oneGene = 0;
		if (!$isGCgene){ #case for normal binning..
			$testKey = ">${cc}_$cnt";
			$curGene = $$ctg2gen{$testKey} if exists($$ctg2gen{$testKey});
		} else {#case for MGS
			$curGene = $cc; $oneGene=1;
		}
		#count up how often each gene is found in already existing Bins (from other samples)..
		while ( $miscnt < 4 || $curGene ne ""){ ##$cnt < 4 ||
			if ($curGene ne ""){
				push(@genes, $curGene);
				$totGenes++;$miscnt=0; my $isMG=0;
				if (exists($gene2COG{$curGene})){
					$MGgenesC++;$isMG=1;
					push(@MGgenes,$curGene);
				} else {
					$nonMGgene++;
				}

				if (exists($gen2Bin{$curGene})){
					$seenGeneCnt++;
					my $HiCnt=0; my $HiBin= "-1";
					foreach my $k (keys %{$gen2Bin{$curGene}}){
						if ($isMG){
							$BinMGcnt{$k}++;
							$seenMGcount++;
						}else {
							$BinCnt{$k} ++;# if ($HiBin ne "-1");
						}
					}
				} else {
					if ($isMG){
						$BinMGcnt{"-1"}++;
					} else {
						$BinCnt{"-1"} ++;
					}
				}
			} else {
				$miscnt++;
				#print "$testKey ";
			}
			$cnt++;
			last if ($oneGene);
			#look for next gene..
			$testKey = ">${cc}_$cnt";
			if (exists($$ctg2gen{$testKey})){$curGene = $$ctg2gen{$testKey} ;} else {$curGene = "";}
		}
	}
	#print "$totGenes $MGgenesC $nonMGgene $seenGeneCnt $seenMGcount\n";
	#print "$seenGeneCnt ";
	#die;
	if ($totGenes ==0 ){die "0 count\n@CCs\n@genes; \n @MGgenes\n";}#";}
	return (\%BinCnt,\%BinMGcnt,\@genes,\@MGgenes,$totGenes,$seenGeneCnt,$nonMGgene);
}




sub gene2GCg2{#
	my ($ar) = @_; #,$MAGn
	my @CCs = @{$ar};
	my @ret;
	my $missedCtg = 0;
	foreach my $cc (@CCs){ #goes over every contig, checks which genes might be present..
		#die "$cc\n";
		if ($cc =~ m/^\d+$/){;
			push (@ret,$cc);
			next;
		} 
		my @spl2 = split(/__/,$cc);
		
		unless (exists($$ctg2gen2{$spl2[0]}{$spl2[1]})){
			next;
		}
		#next unless (exists(
		my @val =  values( %{$$ctg2gen2{$spl2[0]}{$spl2[1]}} );
		push (@ret, @val);
		#foreach my $kk (keys %{$$ctg2gen2{$spl2[0]}{$spl2[1]}}){
			#push (@ret, $$ctg2gen2{$cc}{$kk});
			#print " $$ctg2gen2{$cc}{$kk}";
		#}
		push(@ret, "");
	}
	print "Missed $missedCtg/". scalar (@CCs) . " MAG contigs in gene cat\n" if ($missedCtg);
	return (\@ret,$missedCtg);
}

sub gene2GCg{#
	my ($ar) = @_; #,$MAGn
	my @CCs = @{$ar};
	my @ret;
	foreach my $cc (@CCs){ #goes over every contig, checks which genes might be present..
		#die "$cc\n";
		my $cnt=0; my $curGene = "";my $miscnt=0;
		my $testKey = ">${cc}_$cnt";
		$curGene = $$ctg2gen{$testKey} if exists($$ctg2gen{$testKey});
		while ($cnt < 4 || $miscnt < 10 || $curGene ne ""){
			if ($curGene ne ""){
				push (@ret, $curGene) ;
				#push(@{$gene2MAG{$curGene}},$MAGn);
			} else {
				$miscnt++ ;
			}
			$cnt++;
			$testKey = ">${cc}_$cnt";
			if (exists($$ctg2gen{$testKey})){$curGene = $$ctg2gen{$testKey} ;} else {$curGene = "";}
		}
		#separates contigs by empty entry..
		push(@ret, "");
	}
	#print "$seenGeneCnt ";
	#die;
	return (\@ret);
}


sub countStats{
	my ($hr,$bin,$BQr) = @_;
	my %MBQ = %{$hr};
	my %binQual = %{$BQr};
	if ($MBQ{$bin}{compl}>99){$binQual{g99}++;
	}elsif ($MBQ{$bin}{compl}>95){$binQual{g95}++;
	}elsif ($MBQ{$bin}{compl}>90){$binQual{g90}++;
	}elsif ($MBQ{$bin}{compl}>80){$binQual{g80}++;
	}elsif ($MBQ{$bin}{compl}>60){$binQual{g60}++;
	} else {$binQual{g}++;}
	
	if ($MBQ{$bin}{conta}<1){$binQual{con1}++;
	}elsif ($MBQ{$bin}{conta}<2){$binQual{con2}++;
	}elsif ($MBQ{$bin}{conta}<5){$binQual{con5}++;
	}elsif ($MBQ{$bin}{conta}<10){$binQual{con10}++;
	}elsif ($MBQ{$bin}{conta}<20){$binQual{con20}++;
	}else {$binQual{con}++;}
	return \%binQual;
}



sub summarizeMAGcontent{
	my ($hr) =  @_;

	print "Writing detailed MAG->MGS and MAG->genecat report ($logDir/MAGvsGC.txt.gz)..\n";

	my %MAG2Bin = %{$hr};
	my $MAGreprep = "";
	$MAGreprep .= "Name\tCompleteness\tContamination\tCheckMmodel\tOrigin\n";
	my $MAG2MGScnt=0;
	#create file that reports bin stats and genes (almost all bins)
	open OX,">$logDir/MAGvsGC.txt";
	#HEADER
	print OX "MAG\tMGS\tRepresentative4MGS\tCompleteness\tContamination\tLCAcompleteness\t";
	foreach my $COG (keys %uniCOGs){ print OX $COG."\t";}
	print OX "other_genes\n";
	#fill file..
	foreach my $Doo (sort keys %DOs){
		my @smplIDs = @{$DOs{$Doo}{SmplID}}; my @paths = @{$DOs{$Doo}{wrdir}};
		my $metaGD = getAssemblPath($paths[-1]);
		my $MBf = $metaGD."/Binning/$BinnerShrt/$smplIDs[-1]";	my $MBfQual = $MBf.$cmSuffix;
		next unless (-e $MBfQual);
		my ($hr1,$hr2) = MB2assigns ($MBf,$MBfQual);
		my %MB = %{$hr1};my %MBQ = %{$hr2};
		foreach my $bin (keys %MB){
			my $uniqMBid = "$smplIDs[-1].$bin";
			next unless (exists($MAGlcaq{$uniqMBid}));
			my $MGSid = "?";
			if (exists($binName2{$uniqMBid})){ #$MAG2Bin{$uniqMBid})){
				$MGSid = $binName2{$uniqMBid};#$BinTerm.$MAG2Bin{$uniqMBid} ;
				$MAG2MGScnt++;
			}
			my $star = ""; $star = "*" if (exists($MAG2Bin{$uniqMBid}) && $reprMAGpMGS{$MAG2Bin{$uniqMBid}} eq $uniqMBid);
			print OX "$uniqMBid\t$MGSid\t$star\t$MBQ{$bin}{compl}\t$MBQ{$bin}{conta}\t$MAGlcaq{$uniqMBid}\t";
			#my $ar = gene2GCg($MB{$bin});#,$uniqMBid);
			print "Warning: unidentified MAG: $uniqMBid\n" unless (exists($MAgene{$uniqMBid}));
			my @genes = @{$MAgene{$uniqMBid}}; my @remGene;
			my %G2C; 
			foreach my $gene (@genes){
				if (exists($gene2COG{$gene})){
					push(@{$G2C{$gene2COG{$gene}}}, $gene);
				} 
				#still write to complete gene list,to keep contigs consistent
				push(@remGene,$gene);
				
			}
			
			foreach my $COG (keys %uniCOGs){
				if (exists($G2C{$COG})){
					print OX join(",",@{$G2C{$COG}} ) . "\t";
				} else {
					print OX "\t";
				}
			}
			print OX  join(",",@remGene)."\n";
			
			#report on rep MAG
			if ($star eq "*"){
				$MAGreprep .= "$MGSid\t$MBQ{$bin}{compl}\t$MBQ{$bin}{conta}\t$cmSuffix\t$uniqMBid\n";
			}
		}
	}
	close OX;
	system "rm -f $logDir/MAGvsGC.txt.gz; gzip $logDir/MAGvsGC.txt";
	print "Found $MAG2MGScnt MAGs with MGS assignment\n";
	return $MAGreprep;
}



#original code, kept for backtracking ..
	# for (my $T=0;$T<$maxRounds;$T++){
		# ### deactivated now!
		# last; 
		# ### deactivated now!
		# last unless ($ph1flag);
		# print "\n------------------------ At Tier $T ------------------------\n";
		# foreach my $Doo (@DirSet){
			# #key vars used throughout the loop
			# my @smplIDs;my $uniqMBid ="XX";my %MB ;my %MBQ; my $isCanopy=0;
			
			# #decide on data source (canopies or regular binning)
			# if ($Doo eq "__CanopyMGSset__"){ #get canopies.. only activate at the end of a run
				# #next unless ($T == (@ComplTiers-1));
				# %MB = %Canos; %MBQ = %CanQual;
				# @smplIDs = ("Cano");
				# print "\nAdding MGS from Canopies (N=". scalar(keys(%MB)) . ").";
				# $isCanopy = 1;
			# } else { #get Bin from SemiBin, MetaBat2 etc
				# @smplIDs = @{$DOs{$Doo}{SmplID}}; my @paths = @{$DOs{$Doo}{wrdir}};
				# my $metaGD = getAssemblPath($paths[-1]);
				# my $MBf = $metaGD."/Binning/$BinnerShrt/$smplIDs[-1]";	my $MBfQual = $MBf.$cmSuffix;
				# next unless (-e $MBfQual);
				# my $hr1;
				# print "$MBf\n";
				# ($hr1,$hr2) = MB2assigns ($MBf,$MBfQual);
				# #some stats on bin quality
				# if ($T==0){
					# foreach my $bin (keys %{$hr1}){
						# $binQualHR = countStats($hr2,$bin,$binQualHR) ;
						# $uniqMBid = "$smplIDs[-1].$bin";
						
					# }
				# }
				# #only keep high qual bins..
				# ($hr1,$hr2) = minQualFilter($hr1,$hr2,$ComplTiers[-1], $ContaTiers[-1], $LCAcompl[-1]);
				# %MB = %{$hr1}; %MBQ = %{$hr2};
				# if ($T==0){ #log overall stats for bins..
					# foreach my $bin (keys %{$hr1}){
						# $MAGq{$uniqMBid} = $MBQ{$bin}{compl};$MAGc{$uniqMBid} = $MBQ{$bin}{conta};
					# }
				# }
			# }
			
			# #go through each bin, and identify genes associated to it..
			# #first round: try to find best match
			# foreach my $bin (keys %MB){
				# my %genesInBin;
				# $uniqMBid = "$smplIDs[-1].$bin";
				# #exclusion criteria
				# next if ( exists($regBins{$uniqMBid} ) || exists($ambBin{$uniqMBid}) );
				# #this line decides if the bin qual goes into the next lower tier..
				# if ($MBQ{$bin}{compl}< $ComplTiers[$T] || $MBQ{$bin}{conta}> $ContaTiers[$T] || scalar(@{$MB{$bin}}) == 0){next;}
				# #my @CCs = @{$MB{$bin}};
				# #this routine calculates overlapping genes between already existing bins
				# #	return (\%BinCnt,\%BinMGcnt,\@genes,\@MGgenes,$totGenes,$seenGeneCnt,$nonMGgene,$MGgenesC);
				# my ($BinHR,$BinMGHR,$ar,$arMG, $totGenes,$assignedGenes,$nonMGgene) = countUpBin($MB{$bin},$isCanopy);
				# my %BinCnt = %{$BinHR};my %BinMGcnt = %{$BinMGHR};
				# my @CCs = @{$ar}; #list of genes.. no need to search again
				# my $LCAcomplL = calcLCAcompl(\@CCs,\%LCA);
				# if ($LCAcomplL < $LCAcompl[$T]) {
					# print "\n$uniqMBid LCA qual too low: $LCAcomplL"; next;
				# }
				
				# #list of marker genes in Bin
				# my @MGsG = @{$arMG};
				# my $MGgenes =  @MGsG ;
				# $genesInMAG{$uniqMBid} = \@CCs; # now $MAgene
				
				
				# #print "$totGenes  ".$BinCnt{"-1"}."\n";
				# #decide on which Bins are good...
				# my $BinNumL = $BinNumDeNovo; my $HiCnt=0; my $HiBin= "-1";
				
				# foreach my $bin (keys %BinCnt){ 
					# next if ($bin eq "-1");
					# if ($BinCnt{$bin} > $HiCnt){$HiCnt = $BinCnt{$bin}; $HiBin = $bin;}
				# }
				# my $unknwnFrac = 0;	if (exists $BinCnt{"-1"}){ $unknwnFrac = ($BinCnt{"-1"}/$nonMGgene) ;}
				# my $HiFrac = $HiCnt / $nonMGgene;#best match to other MG

				# #same numbers, but focus on marker genes (more important measure)
				# my $HiMGCnt=0;my $HiBinMG= "-1";
				# foreach my $bin (keys %BinMGcnt){ 
					# next if ($bin eq "-1");
					# if ($BinMGcnt{$bin} > $HiMGCnt){$HiMGCnt = $BinMGcnt{$bin}; $HiBinMG = $bin;}
				# }
				# my $MGunknwnFrac = 0;	if (exists $BinMGcnt{"-1"}){ $MGunknwnFrac = ($BinMGcnt{"-1"}/$MGgenes) ;}
				# my $MGHiFrac = $HiMGCnt / $MGgenes;#best match to other MG
				# my $MGHiFrac2= $MGHiFrac; my $MGHiFrac3 = $MGHiFrac;
				
				# if ($MGHiFrac > 0.6){#test if these are "core core" genes being hit
					# my $sumHit=0;my $sumHitOcc=0; my $Nhit=0; my $BinObs = $BinObs{$HiBinMG};
					# foreach my $gen (@MGsG){
						# my $locSc=0;
						# if (exists($gen2Bin{$gen})){#just catalogue how often these core MGs are occurring..
							# foreach my $sc (values%{$gen2Bin{$gen}}){ $locSc +=$sc;}
							# if (exists($gen2Bin{$gen}{$HiBinMG})){
								# my $HiBinCnt = $gen2Bin{$gen}{$HiBinMG};
								# $sumHit += ($HiBinCnt/$locSc); 
								# $sumHitOcc += $HiBinCnt;
								# $Nhit++; 
							# } 
						# }
							# #ALT: $sumHit += $gen2Bin{$gen}{$HiBinMG}; tests if this gene is often occurring in MGS
						
						
					# }
					# $MGHiFrac2 = ($sumHitOcc/($BinObs*$MGgenes)); #occurrence in core
					# $MGHiFrac3 = ($sumHit/$MGgenes); #not used by other MGS
					# #print "" . ($sumHit/$Nhit) . " $MGHiFrac $BinObs{$HiBinMG} 2:$MGHiFrac2 3:$MGHiFrac3\n$sumHitOcc/($BinObs*$MGgenes)\n$sumHit/$MGgenes\n";
				# }
				
				# #print "$uniqMBid $HiMGCnt $MGHiFrac $MGHiFrac2 H: $HiFrac unkn: $MGunknwnFrac $unknwnFrac\n";
				# if ($legacyV){
					# if ($unknwnFrac > 0.7){##too many unknown, make a new bin
						# if ($AllowNovelMGS[$T]==0){$unknwnBin++; next;}
						# $createPerTier[$T] ++;$binQualTier{$BinNumL} = $T;
						# $novelBin++;$regBins{$uniqMBid}=1;
						# $BinNumDeNovo ++ ; print "making new bin $BinNumL ($MGunknwnFrac)\n";
					# } elsif ( $HiFrac > 0.7 && $HiCnt > 30  ) { #matching to known bin..
						# $BinNumL = $HiBinMG; $matchBin++;$regBins{$uniqMBid}=1;
					# } else {
						# $ambBin{$uniqMBid}=1;next;
					# }
				# } else {
					# if ( $MGunknwnFrac > $unknwnTiers[$T]){##too many unknown, make a new bin
						# if ($AllowNovelMGS[$T]==0){$unknwnBin++; print " U ";next;}
						# $createPerTier[$T] ++;
						# $binQualTier{$BinNumL} = $T;
						# $novelBin++; $CanoNovel++ if ($isCanopy);
						# $regBins{$uniqMBid}=1;
						# $BinNumDeNovo ++ ; 
						# print "\nCreating from $uniqMBid new bin $BinNumL (". int(100*$MGHiFrac) . "%, ". int(100*$MGunknwnFrac) . "%)";
					# } elsif ( $MGHiFrac2 > $OverlapTiers[$T] && $MGHiFrac3 > $OverlapTiers[$T] ) { #matching to known bin..
						# $BinNumL = $HiBinMG; 
						# $matchBin++; $CanoAdd++ if ($isCanopy);
						# $regBins{$uniqMBid}=1;
						# print "\nMatch:  bin $BinNumL (". int(100*$MGHiFrac) . "%, ".  int(100*$MGHiFrac2) . "%, " . int(100*$MGHiFrac3) . "%)";
					# } else {
						# #print "Ambigous Bin\n";
						# $ambBin{$uniqMBid}=1;print " A ";next;
					# }
				# }
				# push(@{$BinGrps{$BinNumL}}, $uniqMBid);
				# $MAG2Bin{$uniqMBid} = $BinNumL;
				# push(@{$Bin2MAG{$BinNumL}}, $uniqMBid);
				# $BinObs{$BinNumL} ++;
				# #my $genesInseted=0;my $newGene=0;
				# foreach my $curGene (@CCs){
					# $gen2Bin{$curGene}{$BinNumL}++;
					# #$genesInseted++; 
					# $genesInBin{$curGene} ++;
					
				# }
				# #print "\n$genesInseted : $HiCnt : $newGene : $unknwnFrac\n" if ($unknwnFrac < 0.7);
				# foreach my $k (keys %genesInBin){ #save which genes were present in multi copy
					# if ($genesInBin{$k} > 1){
						# $GeneMultCopy{$k} ++;
					# }
				# }
			# }
			# $scnt++;
		# }
	#}