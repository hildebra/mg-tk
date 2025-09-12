#!/usr/bin/env perl
#uses a multi sample assembly to bin contigs with metabat
# perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec//MGS.pl /g/scb/bork/hildebra/SNP/GCs/alienGC2/ /local/bork/hildebra/MB2test/ /g/scb/bork/hildebra/SNP/GCs/alienGC2/Canopy2/clusters.txt
# perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec//helpers/MGS/compoundBinning.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ /local/bork/hildebra/MB2test/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Canopy4_AC/clusters.txt
#perl /hpc-home/hildebra/dev/Perl/MATAF3//MGS.pl -GCd /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/ -tmp /scratch/hildebra//GC/GC_Chicken//MAGs/ -nc 24 -canopies /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Canopy2/clusters.txt

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long qw( GetOptions );

#.1: used versioning, just updated Rhcl to (mem efficient) ward clustering, Rhcl checks, Rhcl log reading fixed
#.11: changed rhcl to Rcpphclust package: complete clustering, seems to give better core clusters than ward, also faster now
#.12: checkm2, mem optimizations
#.13: semiBin included
#.14: proGenomes3 added
#.15: GTDBmg added
#.17: process streamlining, HDD usage updates
#.18: updated to "jelly core" MGS clustering, advanced stats on MGS
#.19: set GTDBmg as default, strain resolution updated
#.20: added .LCA MG info to proritize MAGs ("tax clean" MAGs)
#.21: added option to skip Rhcl and DeepCanopies
#.22: -outD flag. complete rework of clusterMAGs.pl script (single processing)
#.23: clusterMAG binary
#.24: no-canopy fix
#.25: 9.12.23: added extraction of representative MGS genome in contigs (highest qual MAG)
#.26: 13.8.24: added -genomesPerFamily flag & function
#.27: 12.11.24: removed necessity for -canopies flag
my $MGSpipelineVersion = 0.26;

use Mods::IO_Tamoc_progs qw(getProgPaths jgi_depth_cmd);
use Mods::GenoMetaAss qw(readMap getDirsPerAssmblGrp readClstrRev unzipFileARezip getAssemblPath systemW gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive qsubSystemWaitMaxJobs);
use Mods::TamocFunc qw ( getFileStr checkMF);
use Mods::geneCat qw(readMG_LCA);
use Mods::Binning qw ( createBin2 createBinCtgs runMetaBat runCheckM runCheckM2 createBinFAA readMGS filterMGS_CM MB2assigns);

#sub MB2assigns;
sub getGoodMBstats;
sub printL;
#sub clusterMB2;
sub getConvergeWcanopy;sub CanopyPrep;
sub filterClustFile;
sub Rhclusts; sub evalRhcl_Bin; sub refine_Rhcl_MGS;
sub invertIndex;
sub createDeepCorrM;
sub replaceLowQualMGS4MAG;



#my $metab2Bin = getProgPaths("metabat2");

my $canBin = getProgPaths("canopy");
my $rareBin = getProgPaths("rare");
my $Rpath = getProgPaths("Rpath");
my $filtDeepCanPost = getProgPaths("filtDeepCan");
my $avx2Constr = getProgPaths("avx2_constraint",0); #"avx2"; #keyword that can be cluster specific


#add this? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4748697/figure/fig-1/
my $inD = "";#$ARGV[0];
my $outD = "";

#$inD.="/" unless($inD =~ m/\/$/);
my $doSubmit = 1;
my $numCore = 4;
my $canCore = 12;
my $memG = 150;#used only for binner
my $PearsCut = 0.3; my $SpearCut = 0.2;#correlation cuts for deep canopies
my $legacyV = 0;#legacy (pre Dec `22) parameters
my $rewrRHCL = 0; #redo RHCL analysis?
my $rewrTAX = 0;
my $rewrClusterMAGs = 0; #redo clusterMAGs analysis
my $rewrDeepCan =0;
my $doStrains = 0;
my $tmpD = ""; 
my $canopyF = "";
my $nodeTmpD = getProgPaths("nodeTmpDir");
my $useCheckM2 = 0; my $useCheckM1 = 1;
my $binSpeciesMG = 2; #defaults to semibin
my $ignoIncomplMAGs = 1; #if 1, will ignore if no Bin in MF3
my $useGTDBmg = "FMG";
my $checkMaxNumJobs = 800; #don't submit more than 800 jobs at once..
my $wait4stone = "";
my $useWeightedMGSscores = 1;
my $doBinCtgsPerFam =0 ; #extract bins for each family (or assembly_grp/sample if family missing)



my $stopAfterCluster = 0; #DEBUG flag!!

#contaminatation parameters for "good" MGS
my $contThre = 5; my $complThre = 90;
#these thresholds are used to include MGS in abundance calculations
my $complThreAbund = 80; my $contThreAbund=5;

my $useRHClust = 0; #outdated clustering with hierachical clusters in R..
my $useDeepCano = 0; #deep canopies? more confusion than worth it imho..

#options to pipeline..
GetOptions(
	"GCd=s"      => \$inD, 				#gene catalog dir
	"outD=s"     => \$outD,				#defaults to $inD/Bin_SB/
	"tmp=s" => \$tmpD,					#temp dir
	"submit=i" => \$doSubmit,			#1:submit jobs, 0: dry run. Default: 1
	"canopies=s" => \$canopyF,			#location of canopy clustering output file (clusters.txt)
	"smallCores=i" => \$numCore,		#cores used for normal jobs (not intensive)
	"bottleneckCores=i" => \$canCore,	#cores for compute intensive jobs
	"useRHClust=i" => \$useRHClust,		#1: do hierachical clustering of MGS genes. Default:0
	"redoRhcl=i" => \$rewrRHCL,			#rewrite R hierachical clusterings
	"redoCluster=i" => \$rewrClusterMAGs,
	"redoDeepCan=i" => \$rewrDeepCan,	#rewrite deep corraltions to Rhcl clusters
	"redoTax=i" => \$rewrTAX,			#rewrite tax annotations
	"MGset=s" => \$useGTDBmg,			#GTDB or FMG, which marker genes are used? Default: GTDB
	"wait4stone=s" => \$wait4stone,     #wait for these files to be created, refers currently exclusively to eggNOG annotations that are needed later
	"mem=i" => \$memG,					#memory used for intensive jobs
	"completeness=i" => \$complThre,	#what qual should final MGS have at least??
	"contamination=i" => \$contThre,	#contamination threshold for accepting MGS
	"strains=i" => \$doStrains,			#1: calc instra species strain phylogenies. Default: 0
	"useCheckM2=i" => \$useCheckM2,		#CheckM2 default qual checking of MAGs/MGS
	"useCheckM1=i" => \$useCheckM1,		#CheckM default qual checking of MAGs/MGS
	"binSpeciesMG=i" => \$binSpeciesMG,	#0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder
	"ignoreIncompleteMAGs=i" => \$ignoIncomplMAGs,	#1: assemblies without MAG calculations are ignored. Default: 1
	"legacy=i" => \$legacyV,			#1: use legacy code as pre Dec `22 (clustering is a bit more muddy, reported abundances slightly different, remember to use -MGset FMG). No longer supported. Default: 0
	"genomesPerFamily=i" => \$doBinCtgsPerFam,
);

#check all correct
checkMF();

#die "$useCheckM2 $useCheckM1\n";

print "\n-------------------------\nMGS v$MGSpipelineVersion\n-------------------------\n";
print "Running legacy version\n" if ($legacyV);
print "Using Rclust\n" if ($useRHClust);

die "-MGset option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");

die "Needs input dir arg (-GCd)!" if ($inD eq "");
#die "$doStrains\n";
#set up basic structures
$tmpD = $inD."/tmp/" if ($tmpD eq "");
my $QSBoptHR = emptyQsubOpt($doSubmit,"");
my %QSBopt = %{$QSBoptHR};

my $singleSample = 0;

die "No MGS set: -binSpeciesMG 0 \n" if ($binSpeciesMG == 0);
my $BinnerShrt = "MB2";
if ($binSpeciesMG == 2){$BinnerShrt = "SB";}#SemiBin
if ($binSpeciesMG == 3){$BinnerShrt = "MD";}
my $COGdir = "FMG";
if ($useGTDBmg eq "GTDB"){ 
	$COGdir = "GTDBmg";
}

$outD = $inD."/Bin_$BinnerShrt/" if ($outD eq "");
$outD .= "/" unless ($outD =~ m/\/$/);
my $logDir = $outD."LOGandSUB/";
my $annoDir = $outD."Annotation/";
my $RsubmDir = "$logDir/subm/";
#checkpoints
my $chkpDir = "$logDir/checkpoints/"; 
my $ABmgsSton = "$chkpDir/abund.mgs.stone";
my $ABmgsSton2 = "$chkpDir/abund.mgs_core.stone";
my $RHcCMstone = "$chkpDir/RhclCan.$BinnerShrt.cm.stone";
my $st1ston = "$chkpDir/Stage1.stone";
my $iniMB2sto = "$chkpDir/$BinnerShrt.cm.stone";
my $EXstone = "$chkpDir/RhclBinExtr.stone";
my $CMrefine = "$chkpDir/CMrefine.stone";
my $GTDBtaxSto = "$chkpDir/GTDBTK.stone";
my $BinExtrSto = "$chkpDir/BinExtr.stone";

#die "$ABmgsSton\n";
#main guide files for MGS
my $txtDoRhcl = ""; $txtDoRhcl = ".Rhcl" if ($useRHClust);
my $txtDoDeep = ""; $txtDoDeep = ".ext" if ($useDeepCano);
my $finalClusters2 = "$outD/$BinnerShrt.clusters";##$txtDoRhcl.can$txtDoRhcl"; 
my $finalClustersW = "$outD/$BinnerShrt.Wclusters"; #unweighted verssion..

my $finalClusters = $finalClusters2."pre";
my $finalClustersFilt = $finalClusters2.".filt";


if (-e "$inD/LOGandSUB/inmap.txt" || -e "$inD/LOGandSUB/GCmaps.inf"){ #this is the outdir of a whole MATAFILER run, or geneCat, doesn't matter
	$singleSample = 0;
	print "Compound Assembly MetaBatting..\n";
} 

system "mkdir -p $RsubmDir $tmpD $outD $logDir $annoDir $chkpDir 1> /dev/null";



my $GCd = "";my $mapF="";
#die Dumper($hrm);	
if (-e "$inD/LOGandSUB/GCmaps.inf"){
	my $tmp = `cat $inD/LOGandSUB/GCmaps.inf`; chomp $tmp;
	$mapF = $tmp;
	$GCd = $inD;
} else{
	$mapF = $inD."LOGandSUB/inmap.txt";
	#($hrm,$asGrpObj) = readMap($inD."LOGandSUB/inmap.txt");
}

#die "$mapF\n";
#figure out which compound assemblies there are..

#infer Assembly dirs & corrsponding bams with several Samples (compound assemblies)
my ($hrD,$hrM) = getDirsPerAssmblGrp($mapF);
my %map = %{$hrM};
my %DOs = %{$hrD};
my @DoosD = sort keys %DOs; #dirs of assembly groups



my $numSamples = scalar(  @{$map{opt}{smpl_order}}  );#@DoosD;
my $useCanopies=1;
if ($numSamples<10 || $canopyF eq ""){$useCanopies=0;}

if (-e "$inD/LOGandSUB/GCmaps.inf"){
	$inD = $map{outDir} if (exists($map{outDir} ));
}
my $ph1flag = 1;
$ph1flag = 0 if (-e "$outD/$BinnerShrt.clusters.obs");
my $MGSimproFlag = 1;
$MGSimproFlag = 0 if (-e $finalClusters2);

open LOG,">$logDir/pipeline.log";
printL "=====================================================\n";
printL "Running MGS v$MGSpipelineVersion pipeline in multi sample mode\n";
printL "GC dir: $inD\n";printL "Map: $mapF\n";
printL "Tmp dir: $tmpD\n";
printL "Ref Canopies: $canopyF\n" if ($canopyF ne "" && $useCanopies);
printL "MGS save dir: $outD\n";
if (!$useCanopies){
	printL "No Canopies used due to N<10 samples (N=$numSamples)\n"; 
}
printL "Cores used: $canCore\n";
printL "Memory: ${memG}G\n";
printL "Rewriting Rhcl\n" if ($rewrRHCL);
printL "Rewriting Taxonomy\n" if ($rewrTAX);
printL "Using CheckM2 quals\n" if ($useCheckM2);
printL "Using Binner $BinnerShrt\n";
printL "=====================================================\n";
my $cmSuffix = ".cm"; $cmSuffix = ".cm2" if ($useCheckM2); 

#clean up
system "rm -r $RsubmDir" if ($rewrRHCL);
system "rm -f $outD/$BinnerShrt.clusters*  $outD/$BinnerShrt.Wclusters*" if ($rewrClusterMAGs);
#my $FMGsubs = `wc -l $GCd/Matrix.$COGdir.mat | cut -f1 -d' '`; chomp $FMGsubs; $FMGsubs = int($FMGsubs);
# a whole lot faster.. but imprecise!
my $FMGsubs = `cat $GCd/$COGdir/*.LCA | wc -l`; chomp $FMGsubs; $FMGsubs = int($FMGsubs); 
if ($FMGsubs < 20){die "$GCd/$COGdir/*.LCA suspiciously small (N=$FMGsubs)\n/ Please ensure correctness\n";} #$GCd/Matrix.$COGdir.mat
#die;
my $CanoDir = $canopyF;$CanoDir=~s/\/[^\/]+$/\//; $CanoDir .= "Bins/";
#die "$CanoDir";
CanopyPrep($canopyF,$CanoDir);
$canopyF = $canopyF.".filt" if ($canopyF ne "");

#run metabat on each assembly group
my $cnt=0; my @jobs;
printL "Found ".scalar(@DoosD) ." assembly groups, ";
if ($ph1flag){
	printL "with Binnings\n"; 
} else {
	printL "calculating $BinnerShrt binnings\n";
}
foreach my $Doo (@DoosD){ #this loops ensures Binner predictions exist for each assembly
	#print "$Doo\n";
	last; #should be done in MATAFILER.. deactivate here..
	last if (!$ph1flag && -e $iniMB2sto);
	my $bef = "";
	my $tmpD2 = $tmpD."$Doo/";
	my $nodeTmpD2 = "$nodeTmpD/checkM/C$Doo/";
	#print "$nodeTmpD2\n";
	$bef .= "mkdir -p $tmpD2\n";# unless (-d $tmpD2);
	#my $allPaths = $DOs{$Doo}{wrdir};
	#my $smplIDtmp = $DOs{$Doo}{SmplID};
	my @smplIDs = @{$DOs{$Doo}{SmplID}};#split /,/,$smplIDtmp;
	my @paths = @{$DOs{$Doo}{wrdir}};#split /,/,$allPaths;
	#next if (@paths <=1);
	my $metaGD = getAssemblPath($paths[-1]);
	my $refFA = $metaGD."/scaffolds.fasta.filt";
	
	
	my $MBout = "$metaGD/Binning/$BinnerShrt/$smplIDs[-1]";
	my $postCmd = "";

	my $CM1done = 0; my $CM2done = 0; my $eBinAssStat=0;
	$CM1done = 1 if (-e "$MBout.cm" );$CM2done = 1 if (-e "$MBout.cm2" );
	$eBinAssStat =1 if (-e "$MBout.assStat");
	#print $MBout." $useCheckM1 $CM1done $useCheckM2 $CM2done $eBinAssStat\n";
	next if ($eBinAssStat && ( ( $useCheckM1 && $CM1done) || ($useCheckM2 && $CM2done ) ) );
	next if ($ignoIncomplMAGs);

	
	#die "$refFA\n";
	#my $refFA = "$inD/$Doo/metag/scaffolds.fasta.filt";
	my $MBcmd = "";
	if ($binSpeciesMG == 1){
		$bef .= jgi_depth_cmd(\@paths,$tmpD2."/depth",95,$numCore,$refFA);# unless (-e );
		$MBcmd = runMetaBat("$tmpD2/depth.jgi.depth.txt",$metaGD."/Binning/$BinnerShrt",$smplIDs[-1],$refFA);
	} elsif ($binSpeciesMG == 2){
		die "MGS.pl::SemiBin not implemented\n";
	} elsif ($binSpeciesMG == 3){
		die "MGS.pl::MetaDecoder not implemented\n";
	}
	#print $bef.$MBcmd;
	$bef = "" if ($MBcmd eq "");
	my $jobName = "Bin$cnt";
	my $mb2Qual = getProgPaths("mb2qualCheck_scr");
	$postCmd = "\n\nrm -rf $nodeTmpD2; mkdir -p $nodeTmpD2;\n$mb2Qual $refFA $MBout $nodeTmpD2 $numCore $useCheckM2\n" unless (-e "$MBout$cmSuffix"  && -e "$MBout.assStat");
	if ($MBcmd eq "" && $postCmd eq "") {next;}#print "next "; next;}
	#next;
	#die "$postCmd\n";
	$postCmd .= "rm -rf $tmpD2\n";
	#print "$MBout\n";
	#die "$bef$MBcmd$postCmd";
	print "$paths[-1]\n";
	my ($jobName2, $tmpCmd) = qsubSystem($paths[-1]."LOGandSUB/${BinnerShrt}_bin.sh",$bef.$MBcmd.$postCmd,$numCore,int($memG/$numCore)."G",$jobName,"","",1,[],\%QSBopt);
	$cnt++;
	push (@jobs, $jobName2);
	#die $paths[-1]."LOGandSUB/MB2_bin.sh";
}
qsubSystemJobAlive( \@jobs,\%QSBopt );


#check that really all cm 's are there
$cnt=0; my @missedMAGs=();
printL "Checking all MAGs are present\n";
foreach my $Doo (@DoosD){
	last if (-e "$iniMB2sto");
	printL "$BinnerShrt MAGs per sample group:\n" if ($cnt==0);
	my @paths = @{$DOs{$Doo}{wrdir}};#split /,/,$allPaths;
	my @smplIDs = @{$DOs{$Doo}{SmplID}};#split /,/,$smplIDtmp;
	my $metaGD = getAssemblPath($paths[-1]);
	my $MBout = "$metaGD/Binning/$BinnerShrt/$smplIDs[-1]";
	if ($ignoIncomplMAGs && (!-e $MBout || -s $MBout ==0 ) ){
		push (@missedMAGs, $smplIDs[-1]);
		#$missedMAGs++; 
		print "no MAG file: $MBout\n";
		next; 
	}
	#check if maybe emtpy Bin?
	if (!-e "$MBout$cmSuffix" ){
		my $ret = `cut -f2 $MBout | sort | uniq`;
		chomp $ret; if ($ret eq "0"){system "touch $MBout$cmSuffix ";}
	}
	my $ret = `wc -l $MBout$cmSuffix`; chomp $ret; $ret =~ s/.*\///;
	#printL "$smplIDs[-1]($Doo)\t$ret\n";
	die "Bin $MBout seems incomplete\n" unless (-e "$MBout$cmSuffix" && -e "$MBout.assStat");
	$cnt++;
}
print "Missed MAGs from " . scalar(@missedMAGs). " assembly groups\n:@missedMAGs\n" if (@missedMAGs);
system "touch $iniMB2sto" unless (-e $iniMB2sto);

#if ($cnt){	print "Waiting for jobs to finish.. restart when done\n";	exit(0);}

#----------------------------------------------------------------------------------------------------
#from here merging of MAGs into MGS

#just writes MAG.$BinnerShrt.assStat.summary, important for  reading in %valMBs
getGoodMBstats() if (!-e $finalClusters2 );#die;


#cluster MAGs based on shared genes between them
if ($ph1flag  || !-e "$outD/$BinnerShrt.clusters" ){
	#clusterMB2() if ($ph1flag  && !-e "$outD/$BinnerShrt.clusters" );
	my $clusscr = getProgPaths("clusterMGS_scr");
	my $canoIncl = ""; $canoIncl = "-canopies $canopyF" if ($useCanopies);
	my $cmd = "$clusscr -GCd $GCd -BinDir $outD -logDir $logDir -binSpeciesMG $binSpeciesMG -MGset $useGTDBmg -cores $numCore -useCheckM1 $useCheckM1 -useCheckM2 $useCheckM2 -legacy $legacyV $canoIncl 1>&2 > $logDir/clusterMGS_scr.log\n";
	print $cmd."\n"; 
	systemW $cmd;
}
die if ($stopAfterCluster); #DEBUGing only!!

#decide between weighted and unweighted scores for binning
if ($useWeightedMGSscores && !-e "${finalClusters2}UW" ){
	if (!-e $finalClustersW){die "Could not find required file $finalClustersW\n";}
	system "mv $finalClusters2 ${finalClusters2}UW";
	system "mv $finalClustersW $finalClusters2";
	#system "touch $finalClustersW.mov";
}

#alt motulizer?
#system "motulizer"
#$hr = readmotulizertable();

#create a core of metabat2 clusters, based on gene occurrence (between different samples)
my $RfilterMB2 = getProgPaths("filterMB2core");
my $postCmd = "$RfilterMB2 $finalClusters2\n" ; #creates $outD/MB2.clusters.core & $outD/MB2.clusters.ext
systemW $postCmd if (!-e $finalClusters2 || !-e "$finalClusters2.core");

if (0&&$useCanopies && $canopyF ne ""  && !-e $finalClusters2) {
	#merge and compare to canopy clusters
	my $iniGuideFile = "$finalClusters2.ext";#created by RfilterMB2 script
	#also get canopy MGS included (if available)
	my $mrgFi = "$iniGuideFile.can"; #contains metabat2, canopy extended + independent canopies 
	#print $mrgFi."\n\n";
	die;
	#outdated.. Canopies are now merged already in $clusscr
	getConvergeWcanopy($canopyF, $iniGuideFile,$mrgFi ) ;  
}


#summarize R clustering bins into gene assignments (and get stats of overlap etc)
if ($useRHClust && (!-e "$finalClusters" || $rewrRHCL)){
	#merge with deep canopies & canopies.. 
	#collect good cluster core... (really?) -> deactivated by default
	my $mrgFi = "$outD/$finalClusters2.ext.can";
	createDeepCorrM($mrgFi,"$outD/$finalClusters2.core") unless (-e $finalClusters);
	print "Creating Rhclusters\n"; 
	#this saves old + new Rhcl MGS into $finalClusters
	my ($href,$aref) = Rhclusts("$RsubmDir",$rewrRHCL,$mrgFi);
	#now evalutate how many genes are completely new, and how many genes are not included
	print "Evaluating R hclusters..\n";
	evalRhcl_Bin("$outD/$finalClusters2.ext",$finalClusters,$href,$aref);
	if (!-e $finalClusters){die "Somthing went wrong creating *.Rhclpre summary\nMissing $finalClusters\n";}
	refine_Rhcl_MGS($finalClusters);
	system "rm -f $RsubmDir/*.mat";
	filterClustFile($finalClusters2,$finalClustersFilt);
	$finalClusters2 = "$outD/$finalClusters2$txtDoRhcl.can$txtDoRhcl"; 
} else {
	print "Rhcl not requested, skipping postclustering in Rhcl\n";
	#is the default now::
	#$finalClusters2 = "$outD/$BinnerShrt.Wclusters";
	$finalClustersFilt = $finalClusters2.".core";
	#systemW "cp $mrgFi $finalClusters2";
}

#die;

#get checkM quality for new Bins
my $binD = "$outD/Genomes/MGS_GC/";system "mkdir -p  $binD" unless (-d $binD);
my $binDctg = "$outD/Genomes/MGS_ctg/";system "mkdir -p  $binDctg" unless (-d $binDctg);
my $binDctgFam = "$outD/Genomes/MGS_ctg_fam/";system "mkdir -p  $binDctgFam" unless (-d $binDctgFam);

#gget repr genomes for each MGS
if (!-e $BinExtrSto){
	print "\n\nCreating reference genome fasta's for all MGS based on contigs in\n$binDctg\n\n";

	if ($doBinCtgsPerFam){
		print "Also creating family-wise ref genomes\n";
		createBinCtgs($binDctgFam,$hrM,"$logDir/MAGvsGC.txt.gz",1);
		#die;
	}

	createBinCtgs($binDctg,$hrM,"$logDir/MAGvsGC.txt.gz",0);
	#do I really need per family genomes??
	
	#die;
	print "Creating reference genome fasta's for all MGS based on gene cat genes in\n$binD\n";
	createBin2($binD,"$finalClustersFilt","$GCd/compl.incompl.95.prot.faa","faa");
	createBin2($binD,"$finalClustersFilt","$GCd/compl.incompl.95.fna","fna");
	
	system "touch $BinExtrSto";
}

#die;

printL "---------------------------------------------------------------------\n";
printL "Stage I clustering done, MGS calculated.\nProgressing to Stage II: annotations, phylogenies and abundances\n";
printL "Using $finalClustersFilt as MGS rep\n";
printL "---------------------------------------------------------------------\n";


#clean up a bit..
system "rm -fr $annoDir/GTDB* $annoDir/kraken2* $annoDir/specI* $GCd//Anno/Tax/SpecI_MGS/" if ($rewrTAX);
system "touch $st1ston" unless (-e $st1ston);

my @jobs2wait=();

#basic quality checks are done at this point
#now get 1) taxonomy 2) phylogeny #) abundance matrix of MGS 4)abundance matrix MGS + specI (to capture unbinned species)

#GTDB tax & kraken2 tax
my $GTDBtaxF = "$annoDir/GTDBTK.tax";
if (!-e $GTDBtaxF || !-e"$annoDir/gtdbtk.summary.tsv" || !-e $GTDBtaxSto){
	my $GTDBtax = getProgPaths("taxPerMGSgtdb_scr");
	my $cmd = "$GTDBtax $binD $canCore $nodeTmpD/GTDBmgs/ $outD\n";
	$cmd .= "mv $outD/GTDBTK.tax $outD/gtdbtk.summary.tsv $annoDir\n\ntouch $GTDBtaxSto";
	#changed mem from 370 to 100 with GTDB-TK 2.1.0
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "150G"; 
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/GTDB.Rhcl.sh",$cmd,$canCore,int(200/$canCore)."G","GTDB_MGS","","",1,[],\%QSBopt);
	$QSBopt{tmpSpace} =$tmpSHDD;
	push(@jobs2wait,$jobName2);
	
}
#wait for checkm/GTDB
qsubSystemJobAlive( \@jobs2wait,\%QSBopt );
die "GTDBtax missing\n$GTDBtaxF\n" if ( !-e $GTDBtaxF);
#if (!-e "$finalClusters2$cmSuffix" && -e "$finalClusters3$cmSuffix"){system "mv $finalClusters3$cmSuffix $finalClusters2$cmSuffix";} #needs to be moved
die "CheckM2 scores missing\n$finalClusters2$cmSuffix\n" if (!-e "$finalClusters2$cmSuffix" );
#die "$cmSuffix\n";

#get only MGS at >$complThre compl, <5 contamination (middle qual), to be used in abundance, strains etc
#create filtered down final liste

#generate taxonmy from kraken2 assignments
if (1){
	my $kr2taxScr = getProgPaths("taxPerMGS_scr");
	my $cmd =  "$kr2taxScr $finalClustersFilt $GCd $annoDir/kraken2\n";# unless (-e "$finalClusters2.LCA");
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/krak2MGS.sh",$cmd,1,int(200/1)."G","KR2_MGS","","",1,[],\%QSBopt) ;
	$QSBopt{tmpSpace} =$tmpSHDD;
}


#redo tree and abundance:
#system "rm -r $outD/between_phylo/ $GCd//Anno/Tax/SpecI_MGS/ $outD/specI.tax";

#generate abundances per MGS
if (0 && !-e "$finalClusters2.matL0.txt"){ #deprecated, use specI based annotations instead..
	invertIndex($finalClusters2,"$finalClusters2.rev") unless (-e "$finalClusters2.rev");
	my $cmd = "$rareBin sumMat -i $GCd/Matrix.mat.gz -o $finalClusters2.mat -refD $finalClusters2.rev -t $numCore\n";
	$cmd .= "rm $finalClusters2.rev\n";
	printL $cmd;
#	systemW $cmd;
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/MGSabund.sh",$cmd,1,int(200/1)."G","AB1_MGS","","",1,[],\%QSBopt) ;
	$QSBopt{tmpSpace} =$tmpSHDD;
}

#die;
#once all tax annotations are done, infer consensus tax for MAGs
#annotate specI's with MAGs added..


unless (-e $ABmgsSton && -e "$outD/Annotation/Abundance/MGS.matL0.txt"){#-e "$GCd//Anno/Tax/${useGTDBmg}_MGS/MGS2speci.txt" ){
	my $specIabu = getProgPaths("specIGC_scr");
	my $cmdSI = "$specIabu -GCd $GCd -cores $canCore -MGS  $finalClustersFilt -MGStax $GTDBtaxF -outD  $outD -MGset $useGTDBmg\n";
	if ($legacyV){
		$specIabu = getProgPaths("specIGC_scr_v0");
		$cmdSI = "$specIabu $GCd $canCore $finalClustersFilt $GTDBtaxF\n";
	}
	$cmdSI .= "touch $ABmgsSton\n";
	$cmdSI .= "cp $GCd//Anno/Tax/SpecI_MGS/MGS2speci.txt $annoDir/specI.tax\n";
	printL "Merge with SpecI & get abundance \n";
	#print "$cmdSI\n";

	my @files = glob ("$GCd/FMG/tax/*tmp.m8");
	#die "@files\n$GCd/FMG/*tmp.m8\n";
	if (scalar(@files) >= 40){
		systemW $cmdSI ; 
	} else {
		print "Incomplete lambda FMG assignments.. submitting job\n";
		my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."/abundMGS.sh",$cmdSI,1,int(300/1)."G","AB2_MGS","","",1,[],\%QSBopt) ;
		$QSBopt{tmpSpace} =$tmpSHDD;
	}
}

unless (-e $ABmgsSton2){
	my $MMLscr = getProgPaths("MAGMGSLCA_scr");
	my $cmdSI2 = "$MMLscr -GCd $GCd -cores 4 -binD $outD;\ntouch $ABmgsSton2\n";
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/abundMGS_core.sh",$cmdSI2,4,int(100/4)."G","AB_MGS_core","","",1,[],\%QSBopt) ;
	$QSBopt{tmpSpace} =$tmpSHDD;
}

#die;
#/hpc-home/hildebra/geneCats/Chicken2/Cultured_genomes/99_ani_dRep/*.fasta
#and just call between MGS treebuild scr
my $treeMem = "100";
my $phyloBetween = getProgPaths("MGSPhyloBetween_scr");
my $baseTreeCmd = "$phyloBetween -GCd $GCd -MGS $finalClustersFilt -mem $treeMem -c $canCore -MSAprogram 4 -fast 0 ";
my $wait4tree = 2; $wait4tree = 2 if ($doStrains);
my $outDphylo = "$outD/between_phylo/";
my $ph1Cmd = "$baseTreeCmd -outD $outDphylo -wait2finish $wait4tree "; #superTree?
my $treedep="";


if (!-e "$outDphylo/phylo/IQtree_allsites.treefile" || !-e "$outD/between_phylo/phylo/IQtree_allsites.pdf"){
	print "Construct phylogeny of all MGS:\n$ph1Cmd\n";
	
	my $refTreeMsg = "\n################\n# If you want to include custom reference genomes, use \n# $baseTreeCmd-outD $outD/customRefs/ -refGenos [refs]\n################\n";
	print $refTreeMsg;
	$refTreeMsg = "#If you want to include custom reference genomes, use $baseTreeCmd -outD $outD/customRefs/ -refGenos [refs]";
	$ph1Cmd .= " -xtraMsg \"$refTreeMsg\";";

	#excute script
	my $ph1OUT = `$ph1Cmd`;
	#my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
	#my ($jobName2, $tmpCmd) = qsubSystem($logDir."/interMGSphylo.sh",$ph1Cmd,1,int(150/1)."G","MGSphylo","","",1,[],\%QSBopt) ;
	#	print "WAITID=$dep\n";
	print $ph1OUT;
	$ph1OUT =~ m/WAITID=(\d+)/;
	$treedep = $1;
}



#die;



#need to rewrite to new format, also check for passed MGS to include in intra-strain analysis
#needs MGSselection.txt
#from here on relies on files from R script, but this could be also auto generated (>80 compl, <5 cota) <- task done (may 20)
#reformat_4phylo($finalClustersFilt) unless (-e "$finalClustersFilt.mgs");


#process ends here unless strains need to be calculated
if ($doStrains == 0){
	printL "\n\nCompound Binning script finished.. no strain analysis set\n";
	close LOG;
	exit(0);
}
print "\n\n########################\nStarting strain delineation MGS\n########################\n";

if ($wait4stone ne ""){
	my $cntWaits=0;
	while ( !-e $wait4stone){
		if ($cntWaits == 0){print "\n\nWaiting for process generating $wait4stone\nIf the process aborted with error, please exit this routine as well\n\n";}
		$cntWaits++;
		sleep(10);
	}
}
#die "XX\n";

my $strain1scr = getProgPaths("MGS_strain1_scr");
my $iniTree = "$outD/between_phylo/phylo/IQtree_allsites.treefile";
#my $prunTree = "$outD/between_phylo/prunned.nwk";
#
#my $ph2Cmd = "$strain1scr $GCd $finalClustersFilt.mgs $canCore $iniTree 0 1\n";#$outD/between_phylo/phylo/IQtree.treefile\n";
my $ph2Cmd = "$strain1scr -GCd $GCd -MGS $finalClustersFilt -MGset $useGTDBmg -maxCores $canCore  -MGSphylo $iniTree -rmMSA 1 -onlySubmit 1 -submit 1 -reSubmit 0 -redoSubmissionData 0 \n#consider adapting options: -rmMSA 0 -presortGenes 1700 -maxGenes 500 -MGSminGenesPSmpl 5\n";#$outD/between_phylo/phylo/IQtree.treefile\n";
printL $ph2Cmd;
#systemW $ph2Cmd;
my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
my ($jobName2, $tmpCmd) = qsubSystem($logDir."/strainMGS.sh",$ph2Cmd,1,int(150/1)."G","strainKickoff",$treedep,"",1,[],\%QSBopt) ;
$QSBopt{tmpSpace} =$tmpSHDD;

#get phylogenies intra-species.. this requires a lot of power and best called from big cluster..
printL "Compound Binning script finished.. run above commands if you want to explore strains\n";
close LOG;
exit(0);




















#####################################################################
#####################################################################
#jsut reformated for better readibility.. 
sub reformat_4phylo{
	my ($FCF) =@_; #, $clusSelHR
	#my @allMGS=keys(%{$clusSelHR});
	#open I,"<$outD/MGSselection.txt" or die "Can't open $outD/MGSselection.txt\n"; 
	#while(<I>){chomp;push(@allMGS,$_);} close I;
	my $hr = readMGS($FCF); my %MGS = %{$hr};
	print "Selected ".scalar(keys %MGS) . " Bins for intrastrain phylo\n";
	open O,">$FCF.mgs" or die $!;
	my $MGcnt=0;
	foreach my $MG (keys %MGS){
		my $MG2 = $MG; $MG2 =~ s/_/:/;#for the MG2dram bins..
		$MG2 = $MG unless (exists($MGS{$MG2}));
		next unless (exists($MGS{$MG2}));
		$MG =~ s/:/_/; #make sure stupid : is completely gone...
		#print $MG."\n";
		$MGcnt++;
		print O $MG."\t".join(",",@{$MGS{$MG2}})."\n";
	}
	close O;
	print "Found $MGcnt MGs for withinstrain analysis\n";
}



sub replaceLowQualMGS4MAG{
	my ($finalClusters,$finalClusters2,$ChkMevalF) = @_;
	#die;
	return if (-e $finalClusters2 && -e $CMrefine);
	
	#good qual MGS, just leave as is
	my $hr = filterMGS_CM($ChkMevalF,$complThre,$contThre,0); #returns NON passing MGS #"$finalClusters$cmSuffix"
	my %MGSset = %{$hr}; my %totMGS;
	$hr = filterMGS_CM($ChkMevalF,$complThre,$contThre,1); #returns PASSING MGS #"$finalClusters$cmSuffix"
	my %goodMGSset = %{$hr};
	#1: read gene clusters and sort by good/bad qual
	#just create hashes of good ($C1)/bad ($MGSbad) MGS (and their genes)
	my %MGSbad; my %C1;
	open I,"<$finalClusters" or die "$finalClusters not present!\n";	
	while (<I>){
		chomp; my @spl = split /\t/; my $curMGS = $spl[0];
		$totMGS{$curMGS} = 1;
		#TODO create list of genes, to use at the end to create finalClusters2
		if (exists($MGSset{$curMGS})){#bad qual MGS
			push(@{$MGSbad{$curMGS}},$spl[1]);#print O2 join("\t",@spl)."\n"; 
		} else {
			#print O join("\t",@spl)."\n";
			push(@{$C1{$curMGS}},$spl[1]);
		}
	}
	close I; #close O; close O2;
	
	#2: link MGS to MAGs
	my %MGS2MAG; my %MAG2MGS;
	open I,"<$outD/$BinnerShrt.clusters.obs" or die "$BinnerShrt.clusters.obs not openable\n";
	while (<I>){
		chomp;my @spl = split /\t/;my $curMGS = $spl[0];
		next if (exists $C1{$curMGS}); #stop replacing these MGS..
#Bin     Observations    QualTier        Members
#MB2bin2 6       0       C2T1.4946,PD11T1.2870,PD12T1.2218,PD8T1.23552,C3T1.7496,C9T1.8945
		if (exists($MGSset{$curMGS})){
			$MGS2MAG{$curMGS} = $spl[3];
			my @spl2 = split /,/,$spl[3];
			foreach my $cc (@spl2){
				$MAG2MGS{$cc} = $curMGS;
			}
		}
	}
	close I;
	
	
	#3: collect cm stats for all MAGS (to compare later to their host MGS
	open I,"<$outD/MAG.$BinnerShrt.assStat.summary" or die "MGS.pl::MB2MAG sum not openable\n";
	my %bestMAG; my $MAGfnd=0;
	while (<I>){chomp;my @spl = split /\t/;
#Sample  MB2     totalL  meanL   ctgN    N20     N50     N80     G1k     G10k   G100k    G1M     CM_Compl        CM_Conta
#C10T0   10265   2590368 28156.1739130435        92      21345   44409   95817  92       67      3       0       94.63   0.00
		my $curMAG = $spl[0].".".$spl[1];
		if (exists($MAG2MGS{$curMAG})){
			my $curMGS = $MAG2MGS{$curMAG};
			#print "$curMGS\n";
			#print "XX $spl[12] $spl[13] XX\n";
			$MAGfnd++;
			if (exists($bestMAG{$curMGS})){#compare to find better MAG
				my $sco = $spl[11]-($spl[12]*2);
				my $sco2 = ( $bestMAG{$curMGS}{compl}-($bestMAG{$curMGS}{conta}*2) ) 
							* ($bestMAG{$curMGS}{N50} / $spl[6] );#takes assembly quality into account: will lower chance of mis-binnings
				if ($sco > $sco2 ){#&& $bestMAG{$curMGS}{N50} < $spl[6] * 0.8){#replace
					$bestMAG{$curMGS}{compl} = $spl[12];$bestMAG{$curMGS}{conta} = $spl[13];$bestMAG{$curMGS}{lengt} = $spl[2];$bestMAG{$curMGS}{MAG} = $curMAG;$bestMAG{$curMGS}{totL} = $spl[3];$bestMAG{$curMGS}{N50} = $spl[6];
				}
			} else {
				$bestMAG{$curMGS}{compl} = $spl[12];$bestMAG{$curMGS}{conta} = $spl[13];$bestMAG{$curMGS}{lengt} = $spl[2]; $bestMAG{$curMGS}{MAG} = $curMAG; $bestMAG{$curMGS}{totL} = $spl[3];$bestMAG{$curMGS}{N50} = $spl[6];
			}
		}
	}
	close I;
	print "Found $MAGfnd MAGs representing ".scalar(keys(%bestMAG)). " MGS, revaluating..\n";
	#now get list of MAGs in bad qual MGS
	my $logStr=""; my $logStr2="";
	my $MGSrepl=0;my$MGSnr=0;
	
	my %replaceMGS;
	foreach my $MGS (keys %MGSset){
		#print "  M  $MGS ";
		if (exists($bestMAG{$MGS}) ){ 
			my $curMAG = $bestMAG{$MGS}{MAG} ;
			my $scoMAG = $bestMAG{$MGS}{compl} - ($bestMAG{$MGS}{conta} *2);
			my $scoMGS = @{$MGSset{$MGS}}[0] - (@{$MGSset{$MGS}}[1]  *2);
			if ( $scoMGS < $scoMAG) { #( @{$MGSset{$MGS}}[1] > $contThre && @{$MGSset{$MGS}}[1] > $bestMAG{$MGS}{conta} ) ||
			#( @{$MGSset{$MGS}}[0] < $complThre && @{$MGSset{$MGS}}[0] < $bestMAG{$MGS}{compl} )){# only check conta, more limiting.. @{$MGSset{$MGS}}[0] > $bestMAG{$MGS}{compl} &&
				$MGSrepl++;
				$logStr .= "$MGS\t$curMAG\t@{$MGSset{$MGS}}[0]\t@{$MGSset{$MGS}}[1]\t$bestMAG{$MGS}{compl}\t$bestMAG{$MGS}{conta}\n";
				$replaceMGS{$curMAG} = $MGS; #this sets the flag to replace MGS with MAG instead..
			}
			#die "$logStr\n";
		} else {
			$MGSnr++;
			$logStr2 .= "$MGS\t$MGS\t@{$MGSset{$MGS}}[0]\t@{$MGSset{$MGS}}[1]\t\t\n";
		}
	}
	foreach my $MGS (keys %goodMGSset){
		$logStr .= "$MGS\t$MGS\t@{$goodMGSset{$MGS}}[0]\t@{$goodMGSset{$MGS}}[1]\t\t\n";
	}
	
	
	print "Replaced $MGSrepl, retained $MGSnr + ".scalar(keys(%C1)). " MGS from Rhcl.\n";
	
	#determined which ones to replace, log this 
	open O,">$logDir/MAGreplMGS.log" or die $!;
	print O "MGS\tMAG\tCompl_MGS\tCont_MGS\tCompl_MAG\tCont_MAG\n";
	print O $logStr;
	print O $logStr2;
	close O;
	
	#now get the genes for the new MAGs and rewrite the MGS asssociations for these..
	my ($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",2); 
	$hr2 = {}; my %ctg2gen = %{$hr1};

	#go through dirs to id the MAG
	my $MAGfound=0; my $totGenesAssigned=0;
	my %genesInBin; #stores the genecat genes in each MAG/MGS
	my %chmLineBin;

	foreach my $Doo (sort keys %DOs){
		my @smplIDs = @{$DOs{$Doo}{SmplID}}; my @paths = @{$DOs{$Doo}{wrdir}};
		my $metaGD = getAssemblPath($paths[-1]);
		my $MBf = $metaGD."/Binning/$BinnerShrt/$smplIDs[-1]";	my $MBfQual = $MBf.$cmSuffix;
		next unless (-e $MBfQual);
		($hr1,$hr2) = MB2assigns($MBf,$MBfQual);
		my %MB = %{$hr1};my %MBQ = %{$hr2};
		
		#go through each bin, and identify genes associated to it..
		
		foreach my $bin (sort keys %MB){
			my $uniqMBid = "$smplIDs[-1].$bin";
			next unless (exists($replaceMGS{$uniqMBid}));
			$chmLineBin{$replaceMGS{$uniqMBid}} = $MBQ{$bin}{line}; 
			$MAGfound++;
			#some stats on bin quality
			#die "$uniqMBid  $MBQ{$bin}{compl}  $MBQ{$bin}{conta}\n";
			#do extract the genes associated to bin and their reps in gene cat
			my @CCs = @{$MB{$bin}};
			foreach my $cc (@CCs){
				my $testKey = ">${cc}_$cnt";
				my $curGene = ""; my $cnt=0; my $miscnt=0;
				$curGene = $ctg2gen{$testKey} if exists($ctg2gen{$testKey});
				#count up how often each gene is found in already existing Bins (from other samples)..
				while ($cnt < 5 ||  $miscnt <10 || $curGene ne ""){
					if ($curGene ne ""){
						#$genesInBin{$curGene} = $testKey; 
						push(@{$genesInBin{$uniqMBid}},$curGene);
						$miscnt = 0;
					} else {$miscnt++;}
					$cnt++;
					$testKey = ">${cc}_$cnt";
					if (exists($ctg2gen{$testKey})){$curGene = $ctg2gen{$testKey};} else {$curGene = "";}
				}
			}
			#print "$uniqMBid Found ". scalar(keys(%genesInBin)) . " genes\n";
			$totGenesAssigned += @{$genesInBin{$uniqMBid}};#scalar(keys(%genesInBin));
		}
	}
	
	print "Found $MAGfound/". scalar(keys(%replaceMGS)) ."/". scalar(keys(%MGSset)) .", with $totGenesAssigned genes (potentially)\n";
	
	#write new genes out
	#final stats on real kept/discarded MGS
	my $MGSkept=0; my $MGSreplaced=0;my $badMGS=0;
	open O,">$finalClusters2" or die "Can't open finalClusters2 $finalClusters2 \n";
	my %collectedMGSset;
	foreach my $uMAG (keys %genesInBin){
		my $curMGS = $replaceMGS{$uMAG};$collectedMGSset{$curMGS} = 1;
		my @genes = @{$genesInBin{$uMAG}};
		foreach my $g (@genes){	print O "$curMGS\t$g\tuMAG\n";	}
		$MGSreplaced++;
	}
	#prior good qual.. just keep
	foreach my $MGS (keys %C1){
		next if (exists($collectedMGSset{$MGS} ));
		$collectedMGSset{$MGS} = 1;
		my @genes = @{$C1{$MGS}};
		foreach my $g (@genes){	print O "$MGS\t$g\tgMGS\n";	}
		$MGSkept++;
	}
	#bad qual that's just bad.
	foreach my $MGS (keys %MGSbad){
		next if (exists($collectedMGSset{$MGS} ));
		$collectedMGSset{$MGS} = 1;
		my @genes = @{$MGSbad{$MGS}};
		foreach my $g (@genes){	print O "$MGS\t$g\tbMGS\n";	}
	}
	close O;
	
	print "Kept $MGSkept, rewrote $MGSreplaced, $badMGS bad qual MGS into $finalClusters2\n"; #. scalar(keys(%collectedMGSset)) .
	
	#transfer checkM vals
	open I,"<$finalClusters$cmSuffix" or die "can't open checkM for rewrite: $finalClusters$cmSuffix\n";
	open O,">$finalClusters2$cmSuffix" or die "can't open rewrite checkM: $finalClusters2$cmSuffix\n";
	my $newC=0; my $oldC=0;
	while (my $lin = <I>){
		chomp $lin; my @spl = split /\t/,$lin;
		if (exists($chmLineBin{$spl[0]})){#replace with new val
			#die "$spl[0]\t$chmLineBin{$spl[0]}";
			print O "$spl[0]\t$chmLineBin{$spl[0]}\n";
			$newC++;
		} else {
			print O "$lin\n";
			$oldC++;
		}
	}
	close O; close I;
	print "Rewrote checkM file, $newC new and $oldC old entries\n";
	
	

	#all done?
	#system "cat $finalClusters2.t $finalClusters2.x > $finalClusters2";
	#system "rm $finalClusters2.t $finalClusters2.n $finalClusters2.x";
	#system "cp $binDpre/* $binD/";

	system "touch $CMrefine";
	return;
}

sub filterClustFile{
	my ($finalClusters,$finalClustersFilt) = @_;
	return if (-e $finalClustersFilt);
	my $hr = filterMGS_CM("$finalClusters2$cmSuffix",$complThreAbund,$contThreAbund); #returns passing MGS $contThre = 5; my $complThre
	my %subs = %{$hr};
	my $totalMGS = `wc -l $finalClusters2$cmSuffix | cut -f1 -d' '`; chomp $totalMGS;
	my $midQMGS = scalar keys %subs;
	print "Filtering $midQMGS MGS from $totalMGS MGS\n";
	open I,"<$finalClusters" or die $!; open O,">$finalClustersFilt" or die $!;
	my $curMGS=""; my $cont=0;my $cnt=0;my$tcnt=0;  
	my $ngenes = 0; my $lcnt=0;
	my %C1; my %C2;
	while (<I>){
		$lcnt++;
		chomp; my @spl = split /\t/; 
		if ($spl[0] ne $curMGS){
			$tcnt++;
			$curMGS = $spl[0];
			if (exists($subs{$curMGS})){
				$cnt++;
				$C1{$curMGS}=1;
				$cont=1;
			} else {
				$C2{$curMGS}=1;
				$cont=0;
			}
		}
		print O join("\t",@spl)."\n" if ($cont);
		$ngenes++;
	}
	close O; close I;
	print "Kept ".scalar (keys (%C1)) .", filtered ".scalar (keys (%C2)) ." MGS\nWrote $ngenes / $lcnt genes\n";
	#print "Found $midQMGS mid qual MGS (of $totalMGS)\n";
}

sub invertIndex{
	my ($in,$out) = @_;
	open I,"<$in" or die "can;t open $in\n";
	open O,">$out" or die "can;t open $out\n";
	while (<I>){
		chomp; my @spl=split /\t/;
		$spl[1] =~ s/:/_/g;
		print O "$spl[1]\t$spl[0]\n";
	}
	close I; close O;
}
sub evalRhcl_Bin{
	my ($oriF,$newF,$rhr,$aar) = @_;
	#print "$oriF\n$newF\n";
	my %ret = %{$rhr};
	my @allB = @{$aar};
	my $hr=readMGS($oriF); my %ori = %{$hr};
	$hr=readMGS($newF); my %new = %{$hr};
	#just go over all bins..
	for my $c (@allB){	
		if (!exists($new{$c}) || !exists($ori{$c})){
			#print "can't find $c\n";
			$ret{inter}{$c} = -1;$ret{newCano}{$c} = -1;$ret{missAll}{$c} = -1;
			next;
		}
		my %union; my %isect;
		my @o = @{$ori{$c}}; my @n=@{$new{$c}};
		foreach my $e (@o, @n) {
			$union{$e}++ && $isect{$e}++
		}
		my $isect = scalar (keys %isect);
		my $union = scalar(keys %union);
		$ret{inter}{$c} = $isect;$ret{newCano}{$c} = @n-$isect;$ret{missAll}{$c} = @o-$isect;
	}
	writeSummaryRhclusts(\%ret,$aar,$newF.".txt");
}

sub refine_Rhcl_MGS(){
	my ($finalClusters) = @_;
	#extract prots/nts per MGS - pre final, used to refine Rhcl clusters that might have gone wrong..
	my $binDpre = "$outD/RhclClust.pre/";system "mkdir -p $binDpre " unless (-d $binDpre);
	if (!-e $EXstone){ 
		#system "rm -r $binD\n" if (-d $binD);
		createBinFAA("$binDpre",$finalClusters,"$GCd/compl.incompl.95.prot.faa","faa");
		createBinFAA("$binDpre",$finalClusters,"$GCd/compl.incompl.95.fna","fna");
		system "touch $EXstone\n";
	}

	#quality check via checkM or checkM2
	my $ChkMevalF = $finalClusters.$cmSuffix;#".cm"; 

	my @jobs2wait = ();
	if ( $useCheckM1 && (!-e $ChkMevalF || !-e $RHcCMstone)){
		printL "running checkM on new Bins..\n";
		#checkm on each Bin proteins
		my $req_CMmem = 200;	my $cmC = "";
		$cmC .= runCheckM($binDpre,$ChkMevalF,"$nodeTmpD/cmMGS/",$numCore,0) ;	
		$cmC .= "\ntouch $RHcCMstone\n";
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."/checkM.Rhclst.sh",$cmC,$numCore,int($req_CMmem/$numCore)."G","ChMrhcl","","",1,[],\%QSBopt);
		push(@jobs2wait,$jobName2);
	}

	#checkM2
	if ( $useCheckM2 && (!-e $ChkMevalF || !-e $RHcCMstone)){
		printL "running checkM2 on new Bins..\n";
		my $req_CMmem = 50;	my $cmC = "";
		$cmC .= runCheckM2($binDpre,$ChkMevalF,"$nodeTmpD/cmMGS/",$canCore,0) ;	
		$cmC .= "\ntouch $RHcCMstone\n";
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."/checkM2.Rhclst.sh",$cmC,$canCore,int($req_CMmem/$canCore)."G","ChMrhcl","","",1,[],\%QSBopt);
		push(@jobs2wait,$jobName2);
	}

	#wait for checkM & checkM2
	qsubSystemJobAlive( \@jobs2wait,\%QSBopt );@jobs2wait = ();

	print "\nReplacing low qual MGS with higher qual MAGs\n";
	#search MGS at 80,5 qual -> replace these with high scoring MAGs
	replaceLowQualMGS4MAG($finalClusters,$finalClusters2,$ChkMevalF);
	system "rm -r $binDpre";

}

sub writeSummaryRhclusts{
	my ($href,$aref,$ofile) = @_;
	my %ret = %{$href};
	my @allB = @{$aref};
	my $cnt=0;
	my $fail=0;my$succ=0;
	open O,">$ofile" or die "Can't open output file $ofile\n";
	print O "Name\tStatus\tQuality\tMaxGenes\tMB2Genes\tClusteredGenes\tMarkGenes\tInternalCorr\tk\tMedianOccu\tCorrStat\tOverlapOri\tNewGenesCano\tGenesOnlyMB2\n";
	for my $c (@allB){
		if (!exists($ret{qual}{$c})){ #failed binning
			print O "$c\tF\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\n";
			$fail++;
		} else {
			print O "$c\t$ret{status}{$c}\t$ret{qual}{$c}\t$ret{MaxG}{$c}\t$ret{SizMB2}{$c}\t$ret{SizNow}{$c}\t$ret{MGs}{$c}\t$ret{cop}{$c}\t$ret{k}{$c}\t$ret{medOcc}{$c}\t$ret{CorMode}{$c}\t$ret{inter}{$c}\t$ret{newCano}{$c}\t$ret{missAll}{$c} \n";
			$succ++;
		}
	}
	close O;
	print "Rhclustering step finished: $succ Successful, $fail Failed\n";
}

sub submitRhcl{
	my ($filesAref,$memFac,$rewr,$logDir2,$inD) = @_;
	my $RpostBin = getProgPaths("PostBinning");
	my @files = @{$filesAref};
	my @jobs;
	#my %ret;
	my $cnt =-1;
	my $ncore = 1;
	my $submissions=0;my $nonSubs=0;
	foreach my $file (@files) {
		#print "$inD/$file\n";
		next unless ($file =~ m/\.mat$/ && -e "$inD/$file" ); 
		qsubSystemWaitMaxJobs($checkMaxNumJobs);

		my $memFac2=1;
		$cnt++;
		#next unless ($file =~ m/^MGS/);
		my $cbin = $file; $cbin =~ s/\.mat//;
		my $submQsub=1;
		my $LOGerrF = "$logDir2/${cbin}.mat.Rhclst.sh.etxt"; 
		my $LOGstdF = "$logDir2/${cbin}.mat.Rhclst.sh.otxt"; 
		if (!$rewr && -e $LOGstdF){#read summary
			my $LOGstr = getFileStr($LOGstdF); 
			#don't resubmit if job is likely still running..
			# resubmission useless, already done or failed once!..
			if ((-e "$inD/${cbin}core.genes" && $LOGstr =~ m/:::Correct:/) || $LOGstr =~ m/:::FAIL::/ || $LOGstr =~ m/Bin not good/|| $LOGstr =~ m/Too few non-zeros/){
				$submQsub=0;
			} #everything else is prob an error..
			if (-e $LOGerrF){
				$LOGstr = getFileStr($LOGerrF); 
				if ($LOGstr =~ m/oom-kill event/){
					$memFac2 = 3;
				}
			}
		}
		#print  "$submQsub $LOGstdF\n";
		#die;
		if ( $submQsub){ 
			my $tmpL = `wc -l $inD/$file | cut -f1 -d' '`; chomp $tmpL;
			my $mem = 10;	#$mem = 20 if ($tmpL > 15e3);	
			my @preCons = @{$QSBoptHR->{constraint}};
			if ($tmpL > 15e3){$mem = 20;$ncore=2;}	
			if ($tmpL > 2e4){$mem = 20;$ncore=3;}	#push(@{$QSBoptHR->{constraint}}, $avx2Constr);
			if ($tmpL > 4e4){$mem = 40;$ncore=4;}	#$mem = 100 if ($tmpL > 5e4);
			if ($tmpL > 5e4){$mem = 50;$ncore=6;}	#$mem = 100 if ($tmpL > 5e4);
			if ($tmpL > 4e4){push(@{$QSBoptHR->{constraint}}, $avx2Constr);}
			#$mem = 150 if ($tmpL > 6e4);	$mem = 200 if ($tmpL > 7e4);
			$mem *= $memFac * $memFac2; $mem = int($mem/$ncore); $mem .= "G";
			#print "M $mem  ";
			my $cmd = "$RpostBin $outD/$finalClusters2 $inD/$file $GCd/$COGdir.subset.cats $Rpath core $ncore\n";
			my $jobName = "Rh$cbin";
			#$QSBoptHR->{useLongQueue} = 1;
			if (0){
				system $cmd;
			}elsif (1){
					my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
					my ($jobName2, $tmpCmd) = qsubSystem($logDir2."/${cbin}.mat.Rhclst.sh",$cmd,$ncore,$mem,$jobName,"","",1,[],\%QSBopt);
					$QSBopt{tmpSpace} =$tmpSHDD;
				$submissions++;
				push (@jobs, $jobName2);
			}
			#reset conditions
			@{$QSBoptHR->{constraint}} = @preCons;
		} else { $nonSubs++;}
	}
	
	qsubSystemJobAlive(\@jobs,\%QSBopt);

	return ($submissions);
}
sub Rhclusts{
#script will either submit R script, or read in summary stat for each post R clustering..
	my ($inD, $rewr,$mrgFi) = @_;
	my $logDir2 = $logDir."/Rhcl/";
	system "mkdir -p $logDir2" unless (-d $logDir2);
	#die;
	opendir DIR,$inD or die $!;
	my %ret;
	my @files = readdir(DIR);
	closedir (DIR);
	my @allBins;
	@files = sort(@files);
	foreach my $file (@files) {
		my $cbin = $file; $cbin =~ s/\.mat//;
		push(@allBins,$cbin);
		$ret{status}{$cbin} = "N";
	}
	#submit Rhcl jobs..
	my $cntRsubs = 1; my $submissions = 1;
	while ($submissions){
		print "Rhcl submission round $cntRsubs\n";
		$submissions = submitRhcl(\@files,1,$rewr,$logDir2,$inD);
		$cntRsubs++;
		if ($cntRsubs > 6){last;}
	} 
	#if ($nonSubs > 0){		#writeSummaryRhclusts(\%ret,\@allBins,$mrgFi.".Rhcl.txt");	}
	###### step 2: read clustering output
	
	my $OX=""; #collects all genes & bins used in total..
	my $coreGenesFnd=0; my $coreGenesMissed=0; my $RHCLsuccess=0;
	foreach my $file (@files) {
		next unless (-e "$inD/$file" && $file =~ m/\.mat$/);
		my $cbin = $file; $cbin =~ s/\.mat//;
		push(@allBins,$cbin);
		$ret{status}{$cbin} = "N";
		if (!-e "$inD/${cbin}core.genes"){#normal not to have this for some cases
			#print "Failed to detect core genes for $inD/${cbin}\n";
			$coreGenesMissed++;
			#next;
		}else {
			$coreGenesFnd ++;
		}
		

		if (!$rewr && -e "$logDir2/${cbin}.mat.Rhclst.sh.otxt"){#read summary
			my $clogRchl = "$logDir2/${cbin}.mat.Rhclst.sh.otxt";
			my $LOGstr = getFileStr($clogRchl); 
			#print "$LOGstr\n";
	#Correct: MB2bin181, G2, MaxG:72503 [S=1193/1536, MG=27, cop=0, k=32199, medOcc=20, spearman]
			if ($LOGstr =~ m/\":::Correct:/){
				#catalog these core genes..
				open IX,"<$inD/${cbin}core.genes" or die "cant open $inD/${cbin}core.genes\n";
				while(my $l=<IX>){chomp $l; $OX .= "$cbin\t$l\tRhcl\n";} close IX;
				$ret{status}{$cbin} = "C";
				$RHCLsuccess++;
			} elsif ($LOGstr =~ m/\":::FAIL::/ || $LOGstr =~ m/Too few non-zeros/){# (\S+), (G\d), MaxG:(\d+) \[S=(\d+)\/(\d+), MG=(\d+), cop=([0-9.]+), k=(\d+), medOcc=([0-9.]+), (.*)\]/){
				#print "Failed RHCL for $cbin"; 
				$ret{status}{$cbin} = "F";
				next;
			} else {
				print "Inconclusive Rhcl report: $clogRchl\n";
				next;
			}
				#collect stats..
			#print "$LOGstr\n";
			#:::Correct: MGS0372, G3, MaxG:1190 [S=1187 (O:0), MG=19 (O:0), cop=0.396, k=2, medOcc=21, complete]
			#:::Correct: MB2bin1, G4, MaxG:4298 [S=1701 (O:1867), MG=83 (O:24), cop=0.102, k=128, medOcc=67, complete]
			$LOGstr =~ m/: (\S+), (G\d), MaxG:(\d*) \[S=(\d+) \(O:(\d+)\), MG=(\d+) \(O:\d+\), cop=([0-9.]+), k=(\d+), medOcc=([0-9.]+), (.*)\]/;
			die "wrong cluster name: \"$1\" != \"$cbin\"\n$LOGstr\n$clogRchl" unless ($1 eq $cbin);
			$ret{qual}{$cbin} = $2; $ret{MaxG}{$cbin} = $3;
			$ret{SizMB2}{$cbin} = $5; $ret{SizNow}{$cbin} = $4;
			$ret{MGs}{$cbin} = $6; $ret{cop}{$cbin} = $7;
			$ret{CorMode}{$cbin} = $10; 
			$ret{k}{$cbin}=$8; $ret{medOcc}{$cbin}=$9;
		
		}
	}
	
	
	###### step 3: create merge clustering file of Rhcl & MB2 & Canopies
	
	#write new clusters (if found any)
	if ($OX eq ""){
		return (\%ret,\@allBins);
	}
	my $hr = readMGS($mrgFi); my %defBins = %{$hr};
	my $xtrWrts=0; my $failedRhcl =0 ;
	open OO,">$finalClusters" or die $!; 
	#write clusters that were good in R clustering
	print OO $OX; 
	#write every cluster that was not good for Rhcl just based on canopies / metabat2..
	foreach my $cl (keys %defBins){
		my $cl2 = $cl; $cl2 =~ s/:/_/;
		if (!exists($ret{status}{$cl2}) || $ret{status}{$cl2} ne "C"){
			$failedRhcl ++ if (exists($ret{status}{$cl2}) && $ret{status}{$cl2} eq "F");
			my @wrt = @{$defBins{$cl}};
			foreach my $wr (@wrt){
				print OO "$cl\t$wr\tcan.ext\n";
			}
			$xtrWrts++;
		}
		
	}
	close OO;
	print "Wrote $xtrWrts / " . scalar(keys %defBins) ." default bins, $failedRhcl failed Rhcl, $RHCLsuccess successes\n";
	$OX="";
	#die;
	return (\%ret,\@allBins);
	
}

sub CanopyPrep{
	my ($inFc,$binCanDir) = @_;
	if ($canopyF eq ""){print"Canopy not requested, will skip step\n";return ;}
	#read genes in canopies..
	my $ChkMevalF = "$inFc.filt$cmSuffix";
	return if (-s $ChkMevalF);
	my %cans; my %canCnts; my %can2gene;
	printL "Prepping Canopy MGS (format, Bin quality)..\n";
	printL "$inFc\n";
	open I,"<$inFc" or die "can't open canopy file $inFc\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#$cans{$spl[1]} = $spl[0];
		$canCnts{$spl[0]} ++;
		push(@{$can2gene{$spl[0]}},$spl[1]);
	}
	close I;
	open O,">$inFc.filt" or die $!;
	my $canCNT=0;
	foreach my $can (keys %canCnts){
		next if ($canCnts{$can} < 700);
		my @loc = @{$can2gene{$can}};
		foreach my $cg (@loc){
			print O "$can\t$cg\n";
		}
		$canCNT++;
	}
	close O;
	printL "Kept $canCNT/". scalar(keys(%canCnts)) . " Canopy MGS\n";
	createBinFAA($binCanDir,"$inFc.filt","$GCd/compl.incompl.95.prot.faa","faa");
	#checkM2
	if ( $useCheckM2 && (!-e $ChkMevalF )){
		printL "running checkM2 on new Canopy MGS..\n";
		my $req_CMmem = 50;	my $cmC = "";
		$cmC .= runCheckM2($binCanDir,$ChkMevalF,"$nodeTmpD/cmCANO/",$canCore,0) ;	
		#$cmC .= "\ntouch $RHcCMstone\n";
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."/checkM2.cano0.sh",$cmC,$canCore,int($req_CMmem/$canCore)."G","ChMrhcl","","",1,[],\%QSBopt);
		push(@jobs2wait,$jobName2);
	}

	#wait for checkM & checkM2
	qsubSystemJobAlive( \@jobs2wait,\%QSBopt );@jobs2wait = ();

}


sub getConvergeWcanopy( $ $ $){
	my ($inFc,$inFm,$mrgF) = @_;
	#read genes in canopies..
	my %cans; my %canCnts; my %can2gene;
	print "$inFc\n";
	open I,"<$inFc" or die "can't open canopy file $inFc\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		$cans{$spl[1]} = $spl[0];
		$canCnts{$spl[0]} ++;
		push(@{$can2gene{$spl[0]}},$spl[1]);
	}
	close I;
	#foreach my $cc(keys %canCnts){		if ($canCnts{$cc} < 200){			delete $canCnts{$cc};			delete $cans{$cc};		}	}
	#read genes in MB2
	
	my %mbs; my %mbCnts; #my %finalMBS;
	open I,"<$inFm" or die "can't open matabat2 filtered file $inFm\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		$mbs{$spl[1]} = $spl[0];
		#$finalMBS{$spl[0]}}{$spl[1]} = 1;
		$mbCnts{$spl[0]} ++;
	}
	close I;
	#find 90% overlaps
	my %merg; my %c2m; my %m2c;
	my %canCnt2;
	foreach my $ck (keys %cans){
		if (exists($mbs{$ck})){
			$c2m{ $cans{$ck} } { $mbs{$ck} } ++;
			$m2c{ $mbs{$ck} } { $cans{$ck} } ++;
			$canCnt2{$cans{$ck}} ++; 
		}
	}
	my $reciCnt=0; my $nreciCnt=0; my $ambiCnt=0; my $indeoenCnt=0; 
	my $defThr=0.7;
	my %indeList;
	foreach my $cc (keys %c2m){
		next if ($canCnts{$cc} < 700);
		my @srtC2M = sort {$c2m{$cc}{$b} <=> $c2m{$cc}{$a}} keys %{$c2m{$cc}};
		my $bestM = $srtC2M[0];
		my $isInde=0;
		my $cF = $c2m{$cc}{$bestM} / ($canCnt2{$cc});
		#test for independence (no overlap) or ambiguity
		if ( $cF < $defThr && $cF > 0.1 && $canCnt2{$cc} > 100){
			$ambiCnt++;
			next;
			#print "Possible ill defined overlap $cF $cc\n";
		} elsif ( $cF < 0.1 || $canCnt2{$cc} < 100 ){#definetly independent canopy, save for later
			$isInde=1;
			$indeoenCnt++;
			$indeList{$cc} = 1;
			#print "$canCnts{$cc} ";
			next;
		}
		#from here it's associated to something
		next unless ($cF >= $defThr && $canCnt2{$cc} > 300); #pretty strong signal
		my @srtM2C = sort {$m2c{$bestM}{$b} <=> $m2c{$bestM}{$a}} keys %{$m2c{$bestM}};
		if ($srtM2C[0] eq $cc){#reciprocal
			#print "reci";
			$reciCnt++;
			#add genes to this MB2 bin to print together later.. actually stupid, that's what the deepclustering is for
			#foreach my $cgene (keys %can2gene{$cc}){
			#	unless (exists($finalMBS{$bestM}}{$cgene})){
			#		$finalMBS{$bestM}}{$cgene} = 1;
			#	}
			#}
			
		} else { 
			#print "nonReci $cc $bestM $c2m{$cc}{$bestM}: @srtM2C\n";
			$nreciCnt++;
			#print "$canCnt2{$cc}\n";
		}
		push(@{$merg{$cc}} , $bestM);
	}
	print "\nreci hit: $reciCnt non-reci hit: $nreciCnt ambigous: $ambiCnt New Bins: $indeoenCnt\n";
	system "rm -f $inFm.can; cp $inFm $inFm.can";
	open O,">>$mrgF";
	foreach my $ic (keys %indeList){
		foreach my $icg ( @{$can2gene{$ic}} ){
			print O "$ic\t$icg\t-1\t-1\t-1\n";
		}
	}
	close O;
	return "$mrgF";
}

sub getGoodMBstats{
	my $logfile = "$outD/MAG.$BinnerShrt.assStat.summary";
	printL "\n\nPer sample MAG stats are in $logfile\n\n";
	return if (-e $logfile && -s $logfile>0);# && -e $finalClusters);
	printL "recording per MAG assembly stats..";
	open L,">$logfile.1" or die $!;
	my $headWritten=0;
	foreach my $Doo (sort keys %DOs){
		my @smplIDs = @{$DOs{$Doo}{SmplID}}; my @paths = @{$DOs{$Doo}{wrdir}};
		my $metaGD = getAssemblPath($paths[-1]);
		my $MBf = $metaGD."/Binning/$BinnerShrt/$smplIDs[-1]";	my $MBfQual = $MBf.$cmSuffix;
		next unless (-e $MBfQual);
		my ($hr1,$hr2) = MB2assigns ($MBf,$MBfQual);
		my %MB = %{$hr1};my %MBQ = %{$hr2}; my %valMBs;
		next if (scalar(keys %MB) == 0); #empty MAG file..
		foreach my $bin (keys %MB){
			if ($MBQ{$bin}{compl}< 80 || $MBQ{$bin}{conta}> 5 ){next;}
			$valMBs{$bin} = "$MBQ{$bin}{compl}\t$MBQ{$bin}{conta}";
		}
		open I,"<$MBf.assStat" or die "Can't open assStat $MBf.assStat\n";
		if (!$headWritten){
			$headWritten=1;
			my $hd = <I>; chomp $hd;
			print L "Sample\t$hd\tCM_Compl\tCM_Conta\n";
		}
		while (my $li = <I>){
			chomp $li;
			$li =~ m/^(\S+)\t/; next unless (exists($valMBs{$1}));
			print L "$smplIDs[-1]\t$li\t$valMBs{$1}\n";
		}
		close I;
	}
	close L;
	system "mv $logfile.1 $logfile;\n";
	print " Done\n";
}




sub printL{
	my ($msg) = @_;
	print $msg;
	print LOG $msg;
}






sub createDeepCorrM{
	my ($guid,$MB2core) = @_;
	if (!$useDeepCano){print "Skipping Deep canopies\n";return;}
	my $corrPrefix = "$BinnerShrt.ext.can";
	my $CorrE = "$outD/$corrPrefix.pear.corr";
	my $CorrS = "$outD/$corrPrefix.spear.corr";
	my $CorrM = "$outD/$corrPrefix.mrg.corr";
	system "rm -rf $CorrM $RsubmDir" if ($rewrDeepCan || $rewrRHCL);
	system "mkdir -p $RsubmDir" unless (-d "$RsubmDir");
	if (!$useCanopies){
		system "touch $CorrM";
		return;
	}
	return if (-e $CorrM);
	system "rm -f $CorrE $CorrS" if ($rewrDeepCan );
	my $canDefOpt = "$canBin -i $GCd/Matrix.mat.scaled.gz   --referenceMB2 $guid ";
	$canDefOpt .= "--maxMB2genes 1000 --dont_use_mmap -n $canCore -b --stop_criteria 0 --filter_max_top3_sample_contribution 1 ";
	$canDefOpt .= " --cag_filter_min_sample_obs 2 --dont_create_progress_stat_file "; #--redundant_guides
	my $postCmd = ""; 
	my $Esto = "$chkpDir/deepE.sto"; my $Ssto = "$chkpDir/deepS.sto"; 
	my @jobs;
	if (!-e $Esto ){
		$postCmd = "$canDefOpt -o $CorrE --max_canopy_dist $PearsCut ; touch $Esto \n";
		my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
		my ($jobName2, $tmpCmd2) = qsubSystem($logDir."/Can.deep1.sh",$postCmd,$canCore,"8G","DEEP1","","",1,[],\%QSBopt);
		$QSBopt{tmpSpace} =$tmpSHDD;
		push(@jobs, $jobName2);
	}
	if (!-e $Ssto){
		$postCmd = "$canDefOpt -o $CorrS --use_spearman --max_canopy_dist $SpearCut; touch $Ssto \n";
		my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
		my ($jobName21, $tmpCmd21) = qsubSystem($logDir."/Can.deep2.sh",$postCmd,$canCore,"8G","DEEP2","","",1,[],\%QSBopt);
		$QSBopt{tmpSpace} =$tmpSHDD;
		push(@jobs, $jobName21);
		#die "\n\nrestart after qsubs have finished \n";
		#print "Waiting for qsubs to finish (abort if they died)..\n";
	}
	qsubSystemJobAlive( \@jobs,\%QSBopt );
	if (-e $CorrE && -e $CorrS){
		$postCmd = "";
		print "Extracting correlating genes from deep canopy ..\n";
		$postCmd .= "$filtDeepCanPost $CorrE $CorrS $MB2core $PearsCut $SpearCut $CorrM\n"; #make sure MB2 core is imputed here
		#$postCmd= "cat $CorrE $CorrS > $CorrM\nrm -rf $RsubmDir/\nmkdir -p $RsubmDir/\n";
		$postCmd .= "$rareBin submatrices -i $GCd/Matrix.mat.scaled.gz -o $RsubmDir -reference $CorrM -t $numCore\n";
		systemW $postCmd;
	}  else {
		die "Can't find $CorrE and $CorrS\n";
	}
	print "Created deep canopies..\n";

}

