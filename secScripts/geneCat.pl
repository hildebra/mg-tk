#!/usr/bin/env perl
#script to create secondary summary stats over MATAFILER runs
#1) creates a gene catalog from predicted genes using cd-hit-est
#2) sums up diamond tables -> moved to helper script "combine_DIA.pl"
#3) makes marker genes to proteins
#4) assigns diamond(FuncAssign) / FOAM to gene finished gene catalog
#usage: ./geneCat.pl [mappingFile] [in/outdir] [numCore] [cdhit-ID]
#ex ./geneCat.pl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/soil_ex.map,/g/scb/bork/hildebra/data2/refData/Soil/PNAS_refs/other_soil.map,/g/scb/bork/hildebra/data2/refData/Soil/howe2014/iowa.map,/g/scb/bork/hildebra/data2/Soil_finland/soil_map.txt /g/scb/bork/hildebra/SNP/GCs/SoilCatv3 1 95 /g/bork3/home/hildebra/data/TAMOC/FinSoil/GlbMap/extraGenes_all.fna /g/bork3/home/hildebra/data/TAMOC/FinSoil/GlbMap/extraGenes_all.faa
#faa generation from gene  catalog 
#ex ./geneCat.pl ? /g/bork1/hildebra/SNP/GC/T2_GNM3_ABR protExtract 
#./geneCat.pl ? /g/scb/bork/hildebra/SNP/GNMass2_singl/ 1 95
#./geneCat.pl MGS /g/scb/bork/hildebra/SNP/GCs/  #canopy clustering
#./geneCat.pl FuncAssign /g/scb/bork/hildebra/SNP/GCs/ #FAOM eggNOG assignment of genes
#./geneCat.pl FMG_extr /g/scb/bork/hildebra/SNP/GCs/SimuGC/  #gets FMGs in separate folder
use warnings;
use strict;
use File::Basename;
use Getopt::Long qw( GetOptions );

use Cwd; use English;
use Mods::GenoMetaAss qw( fileGZe splitFastas readMapS systemW readGFF getAssemblPath resolve_path);  
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive);
use Mods::IO_Tamoc_progs qw(getProgPaths buildMapperIdx);
use Mods::TamocFunc qw(getSpecificDBpaths readTabbed3 checkMF);
use Mods::FuncTools qw(assignFuncPerGene calc_modules);
use Mods::geneCat qw(readGeneIdx  readGeneIdxSpl sortFNA attachProteins attachProteins2 attachProteins3 );

sub geneCatFlow;

sub readCDHITCls;sub readFasta; 
sub writeBucket;  sub announceGeneCat; sub printL;
sub mergeClsSam; sub secondaryCls;
sub cleanUpGC; sub nt2aa; sub clusterFNA;
#sub systemW; #system execution (stops this program upon error)
sub protExtract; #extracts proteins seqs for each cluster
sub combineClstr;#rewrites cd-hit cluster names to my gene Idx numbers
sub FOAMassign;
sub geneCatFunc; sub geneCatFunc_emapper;
sub readSam;
sub prepCDhit; #gets genes from single dirs, sorts into new files, prepares library layout
sub getCanopyDir; sub canopyCluster; #MGS creation
sub krakenTax; #assign tax to each gene via kraken
sub kaijuTax; #assign tax to each gene via kraken
sub specITax;
sub writeMG_COGs; sub ntMatchGC;
sub clusterSingleStep;sub clusterMultiStep;

#.1: : FH
#.27: added external genes support
#.28: gene length filter
#.34: switched to explicit arg parsing and minimap2
#.35: submitLocal added
#.36: bugs, bugs and 
#.37: mmseqs DNA clust
#.38: single step clustering 
#.39: MMseq2 updates
#.40: checkM2, semiBin, proGenomes3, FMG & GTDB MG support
#.41: MG framework, LCA for MGs, decluter check reworked, GTDBmg extraction in GC
#.42: split preprocessing in multiple jobs
#.43: deactivate decluter by default, extrAllE100GC.pl parallelized, more LOGandSUB subfolders for greater clarity, highmem cluster submissions
#.44: eggNOGmapper integration
#.45: 1.5.24: -calcSupplCovSmpls and functionality added
#.46: 23.5.24: added foldseek first trial implementation
#.47: small debugs to fix foldseek integration
#.48: 28.7.24: bugfix for resuming checkpointed GC
#.49: 13.11.24: auto remove canopies if <10 samples
my $version = 0.49;
$| = 1;

my $justCDhit = 1; #always set default to 0, to dangerous otherwise..
my $doSubmit = 1; my $qsubNow = 1;
my $doStrains = 0; #flag if SNP calling was done on metag (and should be processed here)
my $doMags = 1; #flag whether to start canopy clustering, metabat2 & subsequent merging..
my $allinClust = 1; #cluster 5P, 3P, inc and compl in the same step
my $oldNameFolders= -1;
my $doFMGseparation = 1; #cluster FMGs separately?
my $doGeneMatrix =1; #in case of SOIL I really don't need to have a gene abudance matrix / sample
my $numCor = 20; my $totMem = 200; #in G
my $modeStone = ""; #if a stone needs to be touched in that mode..
my $numCor3 = -1; my $numCor0 = -1;
my $totMem3=-1;
my $totMem5=-1;
my $useGTDBmg = "GTDB";#FMG, GTDB
my $toLclustering=1;#just write out, no sorting etc
my $bactGenesOnly = 0; #set to zero if no double euk/bac predication was made
my $canopyAutoCorr=0.15;
my $ignoreIncompleteMAGs = 1;
my $batchNum = -1;
my $smplSep = "__"; #separator of samples and gene id in all MF fastas.. pretty static by now, do not modify!
my $selfScript=Cwd::abs_path($PROGRAM_NAME);#dirname($0)."/".basename($0);
my $rtkFunDelims = "-funcHieraSep \";\" -funcHAnnoAND \",\" -funcAnnoOR \"|\" "; #used to define for rtk how tax annotation strings treat hierachies, and annotations that should be summed (AND) or should be treated as equally likely (OR)

#check all correct
checkMF();

#
#--------------------------------------------------------------program Paths--------------------------------------------------------------
my $magPi = getProgPaths("MAGpipe");
my $mmseqs2Bin = getProgPaths("mmseqs2");
my $cdhitBin = getProgPaths("cdhit");#/g/bork5/hildebra/bin/cd-hit-v4.6.1-2012-08-27/cd-hit
my $vsearchBin = "";#"/g/bork5/hildebra/bin/vsearch1.0/bin/vsearch-1.0.0-linux-x86_64";
my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $hmmBin3 = getProgPaths("hmmsearch");#  hmmer3  <- might need this explicit version if problem
my $hmmBestHitScr = getProgPaths("hmmBestHit_scr");#needed for HMMer based func assignments
#my $tabixBin = "/g/bork5/hildebra/bin/samtools-1.2/tabix-0.2.6/./tabix";
#my $bgzipBin = "/g/bork5/hildebra/bin/samtools-1.2/tabix-0.2.6/./bgzip";
my $pigzBin = getProgPaths("pigz");
my $specIphyloScript = getProgPaths("specIphylo_scr");
my $decluterGC = getProgPaths("decluterGC_scr");
my $kmerScr = getProgPaths("kmerPerGene_scr");
my $rareBin = getProgPaths("rare");#"/g/bork3/home/hildebra/dev/C++/rare/rare";
my $GCcalc = getProgPaths("calcGC_scr");#"perl $thisDir/secScripts/calcGC.pl";
my $sortSepScr = getProgPaths("sortSepReLen_scr");#"perl $thisDir/secScripts/sepReadLength.pl";
my $extre100Scr = getProgPaths("extre100_scr");#"perl $thisDir/helpers/extrAllE100GC.pl";
my $genelengthScript = getProgPaths("genelength_scr");#= "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/geneLengthFasta.pl";
my $GCscr = getProgPaths("geneCat_scr");
my $mini2Bin = getProgPaths("minimap2"); #a lot faster than bowtie2 and better suited..

my $avx2Constr = getProgPaths("avx2_constraint",0); #"avx2"; #keyword that can be cluster specific

# ---------------------------  general constants --------------------------------
my $mini2IdxFileSuffix = ".mmi";
my $countMatrixP = "Matrix";
my $countMatrixF = "$countMatrixP.mat";

my $NodeTmpDir = getProgPaths("nodeTmpDir",0);



#die "$extre100Scr\n";
#my $thisDir = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/";

#databases
my $FOAMhmm = getProgPaths("FOAMhmm_DB");#"/g/bork3/home/hildebra/DB/FUNCT/FOAM/FOAM-hmm_rel1.hmm";
my $ABresHMM = getProgPaths("ABresHMM_DB");#"/g/bork/forslund/morehmms/Resfams.hmm";

#constant path to important files in metagenomic assemblies
my $path2nt = "/genePred/genes.shrtHD.fna";
my $path2GPdir = "/genePred/";
my $path2CS = "/ContigStats/";
my $path2aa = "/genePred/proteins.shrtHD.faa";
my $path2gff = "/genePred/genes.gff";
my $path2FMGids = "ContigStats/FMG/FMGids.txt";
my $COGdir = "COG";
 

if ($bactGenesOnly){
	$path2nt = "/genePred/genes.bac.shrtHD.fna";
	$path2aa = "/genePred/proteins.bac.shrtHD.faa";
	$path2gff = "/genePred/genes.bac.gff";
}
#--------------------------------------------------------------program Paths--------------------------------------------------------------

#global vars
my %FMGfileList;



my $mapF = "?";#$ARGV[0];#"/g/bork5/hildebra/data/metaGgutEMBL/MM.txt";
my $GCdir = "";
my $cdhID = 95; my $minGeneL = 100;
my $extraRdsFNA = "";
#only used in MGS.pl script:
my $useCheckM1= 0; my $useCheckM2 =1;
#only in MGS.pl:
my $binSpeciesMG = 2;#0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder

#$justCDhit = $ARGV[4] if (@ARGV > 4);
#$extraRdsFNA = $ARGV[5] if (@ARGV > 5);#FNA with (predicted) genes, that are to be artificially added to the new gene catalog (and clutered with new genes)
my $extraRdsFAA = "";
#$extraRdsFAA = $ARGV[6] if (@ARGV > 6);#FAA with proteins corresponding to  $extraRdsFNA
#die $cdhID."\n";
my $clustMMseq = 1;#mmseqs2Bin
#$oldNameFolders = $ARGV[4] if (@ARGV>=4 && $ARGV[2] eq "protExtract" );
my $GLBtmp = getProgPaths("globalTmpDir") . "/GC/";
my $tmpDir = $GLBtmp;#getProgPaths("nodeTmpDir") .  "/GC/";
my $tmpDirDef= $tmpDir;
my $mode = "geneCat";
my $fastaSplits="500M";
my $funcAligner = "diamond"; #diamond or foldseek
my $curDB_o = "KGM,TCDB,CZy,ABRc";#NOG,#,ACL"; #"mp3,PTV,KGM,TCDB,CZy,NOG,ABRc,ACL" #default databases to use in functional assignments
my $submitLocal = 1; #default mode now
my $submSys = "";
my $out = ""; my $refDB = ""; #for nt matches via minimap2
my $requireAllAssemblies = 1; #break if some assemblies not present?
my $SmplStart= -1; my $SmplStop= -1; my $SmplBatch= -1; #related to subSmplPrep
my $skipPreCheck = 1;
my $doDecluter = 0;
my $CalcgGneMatSuppl = 1;#calc supplementary coverage samples?



die "geneCat.pl: please provide more arguments\n" if (@ARGV < 2);
GetOptions(
#Directories/files
	"o|GCd=s"  => \$GCdir, #main save location for gene catalog and supporting files
	"tmp=s" => \$tmpDir, #tmp dir, global availalbe
	"glbTmp=s" => \$GLBtmp, #global tmp dir, same as -tmp usually
	"map=s" => \$mapF, #mapping file, can be a combination of several .map files to combine different datasets (e.g. -map file1.map,f2.map)
#run modes
	"m|mode=s" => \$mode, #possible modes: mergeCLs CANOPY specI kraken kaiju FMG_extr FOAM ABR FuncAssign protExtract ntMatchGC geneCat subprepSmpls
#cluster options
	"clusterID=i" => \$cdhID, #identity at which to cluster gene catalog, default: 0.95
	"minGeneL=i" => \$minGeneL, #minimal gene length for gene to be included in gene catalog, default: 100
	"extraGenesNT=s" => \$extraRdsFNA, #add genes (nt) from external sources, e.g. from complete genomes
	"extraGenesAA=s" => \$extraRdsFAA, #add genes (AA) from external sources, e.g. from complete genomes
	"mmseqC=i" => \$clustMMseq, #1: use mmseqs2 instead of CD-HIT for gene clustering
	"decluterMatrix=i" => \$doDecluter, #declutering of gene matrix? by default deactivated, was more useful for canopy based MGS, can intro unwanted biases as long as gene/prots don't get removed
#flow control
	"1stepClust=i" => \$allinClust, #cluster incomplete genes separate?
	"submitLocal=i" => \$submitLocal, #pretty important run mode switch, to submit jobs while geneCat is runnning single core 
	"submSystem=s" => \$submSys,
	"continue|justCDhit=i" => \$justCDhit, #flow control, 1: continue with found files 0: delete existing (partial) gene cat, start again
	"c|cores=i" => \$numCor,
	"c0|cores0=i" => \$numCor0, #specifcally cores only for the big main clustering job..
	"c3|cores3=i" => \$numCor3, #for small jobs that really don't require that much power.. 
	"mem=i" => \$totMem, #max mem
	"stone=s" => \$modeStone,
	"mem3=i" => \$totMem3, #max mem for smaller jobs
#sample processing related
	"calcSupplCovSmpls=i" => \$CalcgGneMatSuppl, #if suppl reads were mapped, report gene abundances in these as separate samples (columns?)
	"oldStyleFolders=i" => \$oldNameFolders, #deprecated. only used for results calculated with an older MATAFILER version
	"requireAllAssemblies=i" => \$requireAllAssemblies, #normally not exposed, continues even if some assemblies not present..
	"sampleBatches=i" => \$batchNum, #how many batches to use for initial accumulation of genes? (200-500 samples per batch recommended)
#Binning/MGS related
	"binSpeciesMG=i" => \$binSpeciesMG, #use MAGs to create MGS? 1= metaBat2, 2=SemiBin, 3=metaDecoder
	"useCheckM2=i" => \$useCheckM2, #1: use checkM2 completeness predictions, Default: 1
	"useCheckM1=i" => \$useCheckM1, #1: use checkM completeness predictions, Default: 0
	"doStrains=i" => \$doStrains, #1: calculate intraSpecific phylogenies on each MGS
	"doMags=i" => \$doMags, #1: start canopy clustering, metabat2 & subsequent merging into MGS
	"canopyAutoCorr=f" => \$canopyAutoCorr, #canopy clustering parameter to filter autocorrelated genes prior to canopy clustering
#Marker Genes/ taxonomy
	"MGset=s" => \$useGTDBmg, #use either FMG or GTDB marker genes to compare and merge MAGs and calculate their abundance
#flags for specific modes
	"out=s" => \$out, #output dir, only used in modes protExtract ntMatchGC 
	"functDB=s" => \$curDB_o, #for FuncAssign mode: functional DBs to annotate gene cat to 
	"refDB=s" => \$refDB, #for ntMatchGC mode: reference fasta DB 
	"fastaSplit=i" => \$fastaSplits, #for FuncAssign mode: split geneCat into chunks to parallelize jobs. Default: 500M 
	"functAligner=s" => \$funcAligner, #either "diamond" or "foldseek"
	"SmplStart=i" => \$SmplStart, #for subprepSmpls
	"SmplStop=i" => \$SmplStop, #for subprepSmpls
	"SmplBatch=i" => \$SmplBatch, #for subprepSmpls
);


checkMF(2);


if ($numCor3 == -1){
	$numCor3 = $numCor;
	if ($numCor > 10){
		$numCor3 = int($numCor/3);
	}
	if ($numCor3 < 1){$numCor3 = 1;}
}
if ($numCor0 == -1){
	$numCor0 = $numCor;
}
if ($totMem3 == -1){
	$totMem3 = $totMem;
	if ($totMem3 > 150){
		$totMem3 = int($totMem/5);
	}
	if ($totMem3 < 40){$totMem3 = 40;}
}
if ($totMem5 == -1){
	$totMem5 = $totMem3;
	if ($totMem5 > 40){
		$totMem5 = int($totMem5/2);
	}
	if ($totMem5 < 10){$totMem5 = 10;}
}




#basic outdir handling..
$GCdir.="/" unless ($GCdir =~ m/\/$/);
if ($GCdir eq "" ){die"no valid output (gene cat) dir specified\n";}
#print "\n\n$GCdir\n";

$GCdir = resolve_path($GCdir);

#die "\n\n$GCdir\n";

if ($tmpDirDef eq $tmpDir){
	$GCdir =~ m/\/([^\/]+)\/?$/;
	$tmpDir .= $1."/"; $GLBtmp.=$1."/";
}

#my $GCdir = $GCdir;#."GeneCatalog/";
my $qsubDir = $GCdir."LOGandSUB/";

my $primaryClusterFNA= "compl.incompl.$cdhID.fna";
my $primaryClusterCLS= "compl.incompl.$cdhID.fna.clstr";

#FMG or GTDB marker genes?
die "-MGset option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");
my $speciesCutoff = "specI_cutoff";
#my $speciesGTDB = "specI_GTDB"; my $speciesDir = "specIPath";my $speciesLink = "specI_lnks"; 
if ($useGTDBmg eq "GTDB"){ 
	print "Using GTDB marker genes\n";
	$path2FMGids = "ContigStats/GTDBmg/marker_genes_meta.tsv";
	$COGdir = "GTDBmg";$speciesCutoff = "GTDB_cutoff";
	#$speciesLink = "GTDB_lnks"; $speciesGTDB = "GTDB_GTDB"; $speciesDir = "GTDBPath";
} else {
	print "Using FMG marker genes\n";
}

#my $logdir = $GCdir."LOGandSUB";
if ($mapF eq ""){
	die "empty mapfile given\n";
}
$mapF =~ s/\/\//\//g;
if ($mapF =~ m/^\??$/){
	if (-e "$qsubDir/GCmaps.inf"){
		$mapF = `cat $qsubDir/GCmaps.inf`;
		die "extracted mapf from $qsubDir/GCmaps.inf\n does not exist:\n$mapF\n" if (!-e $mapF&& $mapF !~ m/,/);
	} else {
		die "Can't find expected copy of inmap in GC outdir: $GCdir\n";
	}
} elsif (-e "$qsubDir/GCmaps.inf" && -e "$qsubDir/GCmaps.inf"){
	
	my $mapFOri = `cat $qsubDir/GCmaps.ori`;$mapFOri =~ s/\/\//\//g;
	if ($mapFOri = $mapF){ #same as in input arg.. great, replace with local copies!
		$mapF = `cat $qsubDir/GCmaps.inf`;
	} else {
		die "Continuing run, but inmap does not seem to match!\nOriginal map: $mapFOri\n-map arg: $mapF\nExiting.. delete gene cat folder before proceeding (or use original map)\n";
	}
}

my %map; my %AsGrps; my @samples; my $numSmpls=0;
if (-e $mapF || $mapF =~ m/,/){ #read mapping file(s)
	print "MAP=".$mapF."\n";
	if (!-e $mapF && $mapF !~ m/,/){die"Could not find map file (first arg): $mapF\n";}
	#die $mapF."\n";
	my ($hr,$hr2) = readMapS($mapF,$oldNameFolders);
	%map = %{$hr};
	$oldNameFolders = $map{opt}{folderStruct} ;
	$GCdir = $map{opt}{outDir} if ($GCdir eq "" && exists($map{opt}{outDir} ));
	@samples = @{$map{opt}{smpl_order}};
	%AsGrps = %{$hr2};
	$numSmpls=scalar(@samples);
	#die "@samples\n";
	#die $map{outDir}."XX\n";
}

#my $defaultsCDH =""; 
my $bucketCnt = 0; my $cnt = 0; my @bucketDirs = ();
my $bdir = $GCdir."B$bucketCnt/";#dir where to write the output files..
system "mkdir -p $qsubDir" unless (-d $qsubDir);
my $QSBoptHR = emptyQsubOpt($doSubmit,"",$submSys);#,"bash"
$QSBoptHR->{qsubDir} = $qsubDir;
#my %QSBopt = %{$QSBoptHR};

#$defaultsCDH = "-d 0 -c 0.$cdhID -g 0 -T $numCor -M ".int(($totMem+30)*1024) if (@ARGV>3);

if ($mode eq "mergeCLs"){#was previously mergeCls.pl
	mergeClsSam($tmpDir,$cdhID,$GCdir);
	exit(0);
} elsif ($mode eq "subprepSmpls"){ #subpart sample gene extractions..
	addingSmpls($SmplStart,$SmplStop,$SmplBatch);
	exit(0);
} elsif($mode eq "CANOPY"){ #create MGS/MGU with canopy clustering
	my $numCor2 = $numCor;
	#if (@ARGV > 2){$numCor2 = $ARGV[2];}
	#die "$numCor2\n";
	canopyCluster($GCdir,$tmpDir,$numCor2);
	exit(0);
}elsif ($mode eq "specI" || $mode eq "kraken" | $mode eq "kaiju" ){ #tax assigns
	my $numCor2 = $numCor;
	#if (@ARGV > 2){$numCor2 = $ARGV[2];}
	if ($mode eq "kraken"){
		krakenTax($GCdir,$tmpDir,$numCor2);
	}elsif($mode eq "specI"){
		specITax($GCdir,$numCor2);
	} else {
		kaijuTax($GCdir,$tmpDir,$numCor2);
	}
	exit(0);
}elsif ($mode eq "FMG_extr"){
	die "deprecated, use helpers/extrAllE100GC.pl\n";
	writeMG_COGs($GCdir);
	exit(0);
} elsif($mode eq "FOAM" || $mode eq "ABR"){ #FOAM functional assignment
	FOAMassign($GCdir,$tmpDir,$mode);
	exit(0);
} elsif($mode eq "FuncAssign"){ #eggNOG/KEGG/CAZy functional assignment
	$QSBoptHR->{doSubmit} = 1;
	#die "$fastaSplits\n";
	my @DBs = split /,/,$curDB_o;
	my $clean=1;
	@DBs = do { my %seen; grep { !$seen{$_}++ } @DBs };
	foreach my $curDB (@DBs){
		my $numCor2 = $numCor;
		if ($curDB eq "mp3"){$numCor2=1;}
		geneCatFunc($GCdir,$tmpDir."/GCanno/",$curDB,$numCor2,$clean);
		$clean=0;
	}
	print "All function assigns done\n";
	exit(0);
} elsif($mode eq "FuncEMAP"){ #eggNOG functional assignment
	$QSBoptHR->{doSubmit} = 1;
	my $clean=1;
	if ($fastaSplits eq "500M"){$fastaSplits="150M";}
	geneCatFunc_emapper($GCdir,$NodeTmpDir."/GCannoEMAP/",$numCor,$clean,$fastaSplits,$modeStone);
	exit(0);
} elsif ($mode eq "protExtract"){#
	protExtract($GCdir,$extraRdsFAA,$out); #GCdir
	exit(0);
} elsif ($mode eq "ntMatchGC"){
	ntMatchGC($GCdir,$refDB,$out);
	exit(0);
} elsif ($mode ne "geneCat" ){
	die "wrong mode called $mode\n";
}


#start logging big genecat run..
system "mkdir -p $qsubDir" unless (-d $qsubDir);
open LOG,">$qsubDir/GeneCat.log";
system "echo \'$version\' > $GCdir/version.txt";
announceGeneCat();

#autoset batchnum
if ($batchNum <1){
	$batchNum = int(scalar(@samples) / 200 )+1;
}
if ($batchNum <1){$batchNum = 1;}


#test if base folders even exist..
my $stoneDir = "$GCdir/checkpoints/";#store checkpoints
my $prepStone = "$stoneDir/GenesCollated.stone";
if (!-e $prepStone && (!-e "$GCdir/B0//compl.fna" || !-e "$GCdir/B0//compl.fna.gz" || !-e "$GCdir/B0/compl.srt.fna.gz") ){
	printL "Genes were not collated, therefore switching to creation mode\n" if ($justCDhit == 1 );
	$justCDhit = 0;
}
#prep base dirs..
if ($justCDhit==0){
	if (-d $GCdir && -d $qsubDir){printL "Warning: outdir $GCdir exists.. delete and recreate? (7s wait)\n"; sleep 7;}
	system ("rm -rf $GCdir\n");#mkdir -p $GCdir/globalLOGs");
} 
system "mkdir -p $stoneDir" unless (-d $stoneDir); 


#copying genes from dirs into prep files..
prepCDhit();

#big main step
geneCatFlow($bdir,$bucketCnt,$GCdir,$map{opt}{outDir});#,$GCdir);

print "FInished GC step\n";
close LOG;
exit(0);










#####################################################################
#####################################################################
sub clusterMultiStep{
	my ($complStone,$incomplStone,$clnLnStone,$cogStone,$bdir,$OutD,$cmd,$dep1) = @_;
	#cd-hit way
	my $copycat=0;
	my $REF = "$tmpDir/compl.$cdhID.fna";

	if( -e $complStone || 
		(-s "$bdir/$primaryClusterFNA" && -s "$bdir/$primaryClusterCLS") ){ #even further along..
	#-s 0.8
		$cmd .= "cp $bdir/compl.$cdhID.fna* $tmpDir\n"; $copycat=1;
		if (-s "$bdir/compl.$cdhID.fna.gz" ){#|| (-s "$bdir/compl.$cdhID.fna.clstr.gz" && -s "$bdir/compl.$cdhID.fna.ctsv.gz")){
			$cmd .= "gunzip $tmpDir/compl.$cdhID.fna*\n" ;
		} 
	}else {
		#$cmd .= sortFNA($bdir,"compl",$toLclustering,$tmpDir,$numCor);
		#$cmd .= $cdhitBin."-est -i $bdir/compl.fna -o $tmpDir/compl.$cdhID.fna  -n 9 -G 1 -aS 0.95 -aL 0.6 $defaultsCDH \n" ;
		$cmd .= clusterFNA( "$bdir/compl.fna", $REF ,0.9,0.6,"$cdhID",$numCor,0,$tmpDir."/mainCL/",$clustMMseq,$totMem);
		$cmd .= "touch $complStone\n";
		$cmd .= "cp $tmpDir/compl.$cdhID.fna* $bdir\n" unless ($copycat);
	}
	if ($submitLocal){
		if (!$copycat){
			$QSBoptHR->{useLongQueue} = 1;
			my @preCons = @{$QSBoptHR->{constraint}};
			push(@{$QSBoptHR->{constraint}}, $avx2Constr) if ($clustMMseq);#--constraint=sse4
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep,$qcmd) = qsubSystem($qsubDir."mainCDhit.sh",$cmd,$numCor,int($totMem/$numCor)."G","MCLGC","","",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$QSBoptHR->{useLongQueue} = 0;
			@{$QSBoptHR->{constraint}} = @preCons;

			qsubSystemJobAlive( [$dep],$QSBoptHR ); 
			die "Can't find Stone $complStone\n" unless (-e $complStone);
		} else { systemW $cmd;}
		$cmd = "";
	}
	
	#5', 3' complete genes and incompletes
	
	my $bwtIdx = $REF.".bw2";
	my $bwtIncLog = "$qsubDir/bowTie_incompl.log"; my $bwt35Log = "$qsubDir/bowTie_35.log";
	my $bwtIdxB0 = "$bdir/SAM/compl.$cdhID.fna.bw2";
	system "mkdir -p $bdir/SAM";
	
	my $bwtCore = $numCor;
	
	#$bwtCore=16;# if ($bwtCore > 25);
	
	my @miniParJobs = (); my $mapIdxJob="";
	#build bowtie index?
	if ( !-s "$bdir/SAM/35compl.$cdhID.align.sam"  || !-s "$bdir/SAM/P35compl.NAl.pre.$cdhID.fna" || 
			!-s "$bdir/SAM/incompl.$cdhID.align.sam" || !-s "$bdir/SAM/incompl.NAl.pre.$cdhID.fna" ){
		my ($tmpCmd,$bwtIdxT) =  buildMapperIdx($REF,$numCor3,1,3);
		$cmd .=$tmpCmd;
		$bwtIdx = $bwtIdxT;
		if ($submitLocal){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep,$qcmd) = qsubSystem($qsubDir."mapperIdx.sh",$cmd,$numCor3,int($totMem3/$numCor)."G","IDXmGC","","",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$cmd = ""; $mapIdxJob = $dep;
		}

		#unless (-e $bwtIdxB0.".rev.2.bt2" && -e $bwtIdxB0.".4.bt2"&& -e $bwtIdxB0.".2.bt2"){
		#	$cmd .= $bwt2Bin."-build --threads $bwtCore -q $REF $bwtIdx\n" ;
#			} else {
#				$cmd .= "cp $bwtIdxB0* $tmpDir\n";
#			}
	}
	#die "$bwtIdx \n $REF\n";
	my $longbwt2opt = " -L 15 -i S,0,2.50 --local --norc  --no-hd --no-sq "; #--no-unal
	my $stdBowtie2opt = " -sensitive --local --norc  --no-hd --no-sq "; #--no-unal
	my $mini2Base = "$mini2Bin -2 -a -t $numCor3 --secondary=no -x asm20 "; #--sam-hit-only
	my $xxtra = "";
	#p35 incomplete genes
	if ( !-s "$bdir/SAM/35compl.$cdhID.align.sam"  ){
		if (-e "$bdir/5Pcompl.fna.gz" && !-e "$bdir/5Pcompl.fna"){
			$cmd .= "zcat $bdir/5Pcompl.fna.gz > $tmpDir/35Pcompl.fna\n" 
		} else {
			$cmd .= "cat $bdir/5Pcompl.fna > $tmpDir/35Pcompl.fna\n";#>> $bdir/incompl.fna \n";
		}
		if (-e "$bdir/3Pcompl.fna.gz" && !-e "$bdir/3Pcompl.fna"){
			$cmd .= "zcat $bdir/3Pcompl.fna.gz >> $tmpDir/35Pcompl.fna\n\n";
		} else {
			$cmd .= "cat $bdir/3Pcompl.fna >> $tmpDir/35Pcompl.fna\n\n";
		}
		$cmd .= "touch $tmpDir/35compl.$cdhID.align.sam \n";
		if (0){
			#take care of long reads - only needed for bowtie2
			$cmd .= "$sortSepScr 8000 $tmpDir/35Pcompl.fna\n";
			$cmd .= $bwt2Bin.$longbwt2opt;
			$cmd .= "--un $tmpDir/P35compl.NAl.pre.$cdhID.fna.long -p 4 -x $bwtIdx -f -U  $tmpDir/35Pcompl.fna.long > $tmpDir/35compl.$cdhID.align.sam 2> $bwt35Log\n";
			#and bulk of reads
			$cmd .= $bwt2Bin.$stdBowtie2opt." -p $bwtCore ";
			$cmd .= "--un $tmpDir/P35compl.NAl.pre.$cdhID.fna -x $bwtIdx -f -U  $tmpDir/35Pcompl.fna >> $tmpDir/35compl.$cdhID.align.sam 2>> $bwt35Log\n";
			
			#fix missing newlines
			$cmd .= "cat $tmpDir/P35compl.NAl.pre.$cdhID.fna.long >> $tmpDir/P35compl.NAl.pre.$cdhID.fna\n rm $tmpDir/P35compl.NAl.pre.$cdhID.fna.long\n";
			$cmd .= "\nsed -i -r 's/([ACGT])>/\\1\\n>/g' $tmpDir/P35compl.NAl.pre.$cdhID.fna\n";
			$xxtra = "$tmpDir/P35compl.NAl.pre.$cdhID.fna";
		} else {#mini2 way
#				foreach my $cog ( @COGlst){#do COGs separate .. too complicated, just ignore for now..
#					$cmd .= "$mini2Base $FMGFL2{$cog}$mini2IdxFileSuffix $tmpDir/35Pcompl.fna | grep -v '^\@' > $tmpDir/35compl.$cdhID.align.$cog.sam \n"; #$samBin view 
#				}
			$cmd .= "if [ -s $tmpDir/35Pcompl.fna ] ; then $mini2Base $bwtIdx $tmpDir/35Pcompl.fna | grep -v '^\@' > $tmpDir/35compl.$cdhID.align.sam ; fi\n"; #$samBin view 
		}
		$cmd .= "cp $xxtra $tmpDir/35compl.$cdhID.align.sam $bdir/SAM/\n" ;
		
		#$bwtCore=50 if ($bwtCore > 50);
	 } else {
		$cmd .= "cp $bdir/SAM/35compl*.sam  $tmpDir/\n"; # $bdir/SAM/P35compl.NAl.pre.$cdhID.fna
		if ($submitLocal){systemW $cmd;$cmd="";}
	 }
	if ($submitLocal && $cmd ne ""){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		my ($dep,$qcmd) = qsubSystem($qsubDir."mini35.sh",$cmd,$numCor3,int($totMem3/$numCor3)."G","m35GC",$mapIdxJob,"",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$cmd = ""; push(@miniParJobs,$dep); #qsubSystemJobAlive( [$dep],$QSBoptHR ); 
	}

	$xxtra = "";
	if (!-s "$bdir/SAM/incompl.$cdhID.align.sam" ){
		my $queryFNA = "$bdir/incompl.fna";
		if (-e "$bdir/incompl.fna.gz" && !-e $queryFNA){
			#$cmd .= "gunzip $bdir/incompl.fna.gz\n";
			$queryFNA = "$bdir/incompl.fna.gz";
		}
		$cmd .= "touch $tmpDir/incompl.$cdhID.align.sam \n";
		#take care of long reads
		if (0){
			$cmd .= "$sortSepScr 8000 $bdir/incompl.fna\n";
			$cmd .= $bwt2Bin. $longbwt2opt;
			$cmd .= "--un $tmpDir/incompl.NAl.pre.$cdhID.fna.long -x $bwtIdx -f -U $bdir/incompl.fna.long > $tmpDir/incompl.$cdhID.align.sam 2> $bwtIncLog\n";
			$cmd .= $bwt2Bin.$stdBowtie2opt." -p $bwtCore ";
			$cmd .= "--un $tmpDir/incompl.NAl.pre.$cdhID.fna -x $bwtIdx -f -U $bdir/incompl.fna >> $tmpDir/incompl.$cdhID.align.sam 2>> $bwtIncLog\n";
			$cmd .= "cat $tmpDir/incompl.NAl.pre.$cdhID.fna.long >> $tmpDir/incompl.NAl.pre.$cdhID.fna\n rm $tmpDir/incompl.NAl.pre.$cdhID.fna.long\n";
			$cmd .= "sed -i -r 's/([ACGT])>/\\1\\n>/g' $tmpDir/incompl.NAl.pre.$cdhID.fna\n";
			$xxtra = "$tmpDir/incompl.NAl.pre.$cdhID.fna";
		} else {#mini2 now
			$cmd .= "if [ -s $queryFNA ]; then $mini2Base $bwtIdx $queryFNA | grep -v '^\@' > $tmpDir/incompl.$cdhID.align.sam ; fi\n"; # $samBin view 2> $bwt35Log
		}
		$cmd .= "cp $xxtra $tmpDir/incompl.$cdhID.align.sam $bdir/SAM/\n" ;
	 } else {
		$cmd .= "cp $bdir/SAM/incompl*.sam  $tmpDir/\n"; #$bdir/SAM/incompl.NAl.pre.$cdhID.fna 
		if ($submitLocal){systemW $cmd;$cmd="";}

	 }

	if ($submitLocal && $cmd ne ""){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		my ($dep,$qcmd) = qsubSystem($qsubDir."miniIncom.sh",$cmd,$numCor3,int($totMem3/$numCor3)."G","mInGC",$mapIdxJob,"",1,[],$QSBoptHR); $cmd = "";
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		push(@miniParJobs,$dep); #
	}
	
	#merge cluster track files (as they all match to same ref DB set)
	#$cmd .= "rm -f $bwtIdx"."*\n"; #<- leave in
	#if (!-s "$bdir/SAM/incompl.NAl.pre.$cdhID.fna" || !-s "$bdir/SAM/P35compl.NAl.pre.$cdhID.fna" ){
	$cmd .= "perl $selfScript -mode mergeCLs -o $OutD -MGset $useGTDBmg -tmp $tmpDir -1stepClust $allinClust -clusterID $cdhID -c $numCor -map $mapF\n"; #1 was BIG before.. but kinda useless paramenter now
	#}
	#die $cmd."\n";
	$cmd .= "\nsed '/^\$/d' $tmpDir/$primaryClusterFNA > $tmpDir/$primaryClusterFNA.tmp;rm $tmpDir/$primaryClusterFNA;mv $tmpDir/$primaryClusterFNA.tmp $tmpDir/$primaryClusterFNA\n";
	$cmd .= "touch $clnLnStone\n";
	$cmd .= "$pigzBin -f -c -p $numCor $tmpDir/$primaryClusterFNA > $bdir/$primaryClusterFNA.gz\n"; #just make sure this is backed up..
	$cmd .= "$pigzBin -f -c -p $numCor $tmpDir/$primaryClusterCLS > $bdir/$primaryClusterCLS.gz\n"; #just make sure this is backed up..
#		$cmd .= "rm -r $bdir/SAM/\n"; #not really needed any longer, delete..
	$cmd .= "mv $tmpDir/log/Cluster.log $qsubDir\n"; #this is from merge script..
	$cmd .= "touch $incomplStone\n";
	#still multi Core..
	if ($submitLocal){
		$QSBoptHR->{useLongQueue} = 1;
		my @preCons = @{$QSBoptHR->{constraint}};
		push(@{$QSBoptHR->{constraint}}, $avx2Constr) if ($clustMMseq);
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		my ($dep,$qcmd) = qsubSystem($qsubDir."mergeIncom.sh",$cmd,$numCor,int($totMem3/$numCor)."G","mEInGC",join(";",@miniParJobs,$dep1),"",1,[],$QSBoptHR); $cmd = "";
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
		@{$QSBoptHR->{constraint}} = @preCons;
		qsubSystemJobAlive( [$dep],$QSBoptHR ) ; 
		die "Can't find Stone $incomplStone\n" unless (-e $incomplStone);
	}
	return $cmd;
}

sub gzifelscat{
	my ($inF) = @_;
	if (-e "$inF.gz" && !-e "$inF"){
		return "zcat $inF.gz";
	} else {
		return "cat $inF";
	}

}
sub clusterSingleStep{
	my ($complStone,$incomplStone,$clnLnStone,$cogStone,$bdir,$OutD,$cmd,$dep1) = @_;
	my $DB = "$tmpDir/compl.35inc.fna.gz";
	return "" if (-e $incomplStone);
	
	#these steps dont use local SSD, tmpsapce is on scratch..
	my $preHDDspace = ${$QSBoptHR}{tmpSpace};
	${$QSBoptHR}{tmpSpace} = "0";#"${totMem}G"; #$totMem #doesn't need much, stores on scrach
	my $dep =""; my $qcmd = "";
	
	my $OFcompl = "$bdir/compl.fna.gz";my $OFincompl = "$bdir/incompl.fna.gz";my $OF5in = "$bdir/5Pcompl.fna.gz";my $OF3in = "$bdir/3Pcompl.fna.gz";
	
	$cmd .= "touch $tmpDir/compl.$cdhID.fna\n";
	$cmd .= "echo \"Concatenating input files of complete & incomplete genes..\"\n";
	$cmd .= "rm -f $DB; cp $OFcompl $DB; cat $OFincompl $OF5in $OF3in >> $DB\n";
	#previously worked on unzipped files..
	#$cmd .= gzifelscat("$bdir/compl.fna"). " > $DB;\n";$cmd .= gzifelscat("$bdir/5Pcompl.fna"). " >> $DB;\n";
	#$cmd .= gzifelscat("$bdir/3Pcompl.fna"). " >> $DB;\n";$cmd .= gzifelscat("$bdir/incompl.fna"). " >> $DB;\n";
	#now on zipped ones
	$cmd .= "touch $clnLnStone\n";
	
	
	#already done? skip this..
	if (-e $complStone || (-e $clnLnStone && -e $DB) ){
		$cmd = "" ;
	} else {#submit as single core job..
		#my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
		#($dep,$qcmd) = qsubSystem($qsubDir."prepFiles.sh",$cmd,1,"10G","prepFiles","","",1,[],$QSBoptHR); $cmd = "";
		#$QSBopt{tmpSpace} =$tmpSHDD;
	}
	#cat stuff can run on single core..
	if ($submitLocal && $cmd ne ""){
		#die "$cmd\n";
		print "Will merg pre files..\n";
		systemW $cmd;$cmd="";
	}

	$cmd .= clusterFNA( "$DB", "$tmpDir/$primaryClusterFNA" ,0.0,0.0,"$cdhID",$numCor0,0,$tmpDir."/fullCL/",$clustMMseq,$totMem);
	
	$cmd .= "$pigzBin -f -c -p $numCor0 $tmpDir/$primaryClusterFNA > $bdir/$primaryClusterFNA.gz\n"; #just make sure this is backed up..
	$cmd .= "$pigzBin -f -c -p $numCor0 $tmpDir/$primaryClusterCLS > $bdir/$primaryClusterCLS.gz\n"; #just make sure this is backed up..
	$cmd .= "touch $complStone $qsubDir/Cluster.log\n";
	if (-e $complStone){
		$cmd = "mkdir -p $tmpDir\n" ;
		$cmd .= gzifelscat("$bdir/$primaryClusterCLS"). " > $tmpDir/$primaryClusterCLS;\n";
		$cmd .= gzifelscat("$bdir/$primaryClusterFNA"). " > $tmpDir/$primaryClusterFNA;\n";
		if ($submitLocal){systemW $cmd;$cmd="";}
	}
	
	if ($submitLocal && $cmd ne ""){
		my @preCons = @{$QSBoptHR->{constraint}};
		#die "@preCons\n";
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		push(@{$QSBoptHR->{constraint}}, $avx2Constr) if ($clustMMseq);
		$QSBoptHR->{useHiMemQueue} = 1 if ($totMem > 1000);
		my ($dep,$qcmd) = qsubSystem($qsubDir."fullClust.sh",$cmd,$numCor0,int($totMem/$numCor0)."G","fullGC",$dep,"",1,[],$QSBoptHR); $cmd = "";
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$QSBoptHR->{useHiMemQueue} = 0;
		@{$QSBoptHR->{constraint}} = @preCons;
		qsubSystemJobAlive( [$dep,$dep1],$QSBoptHR ) ; 
		die "Can't find Stone $complStone\n" unless (-e $complStone);
	} elsif ($dep1 ne "") {
		qsubSystemJobAlive( [$dep1],$QSBoptHR ) ; 
	}
	die "clustering failed\n" unless (-e $complStone);
	${$QSBoptHR}{tmpSpace} = $preHDDspace;
	
	$cmd .= "\n\nperl $selfScript -mode mergeCLs -MGset $useGTDBmg -1stepClust $allinClust -o $OutD -tmp $tmpDir -clusterID $cdhID -c $numCor0 -map $mapF\n"; #still need to add the COG genes..
	$cmd .= "touch $incomplStone\n";

	if ($submitLocal){
		die "Can't find $cogStone\n" unless (-e $cogStone);
		systemW $cmd;$cmd="";
	}

	#die $cmd;
	
	return $cmd;

}


sub geneCatFlow($ $ $ $ ){
	my ($bdir,$nm,$OutD,$assDirs) = @_;
	my $cmd = "";
	my $totMemL = $totMem; if ($totMem > 50){$totMemL=50;}
	my $totMemS = int($totMem/4); if ($totMemS <50){$totMemS = 50;}
	#my %FMGfileList = %{$hr}; <- is globally defined 
	
	
		#stones to track progress..
	my $complStone = "$stoneDir/2.complCl.stone";
	my $cogStone = "$stoneDir/2.COGclust.stone";
	my $incomplStone="$stoneDir/3.incomplCl.stone";
	my $clnLnStone = "$stoneDir/1.clnLines.stone";
	my $matrixSton = "$stoneDir/5.matrix1.stone";
	my $moveStone = "$stoneDir/6.moveRel.stone";
	my $protStone = "$stoneDir/6.protExtr.stone";
	my $declStone = "$stoneDir/7.decluter.stone";
	my $FMGstone = "$stoneDir/8.FMG.stone";
	my $cogLCAStone = "$stoneDir/8.COG.LCA.stone";
	my $geneStatsStone = "$stoneDir/9.geneStats.stone";
	my $krakStone = "$stoneDir/10.krak.stone";my $SIstone="$stoneDir/10.specI.stone"; 
	my $canopyStone = "$stoneDir/10.canopy.stone";
	my $funcStone = "$stoneDir/10.func.stone";
	my $emapStone = "$stoneDir/10.emap.stone";
	my $preMGSstone="$stoneDir/11.preMGS.stone";
	
	
	#get list of potential COGs
	if (-d "$bdir/$COGdir") { #COGs were created
		opendir(DIR, "$bdir/$COGdir/") or die $!;
		my @cogfiles = grep {/^preclus\..*\.fna$/ && -f "$bdir/$COGdir/$_" } readdir(DIR); close DIR;
		for (my $i=0; $i<@cogfiles; $i++){
			$cogfiles[$i] =~ m/^preclus\.(.*)\.fna$/;
			$FMGfileList{$1} = "$bdir/$COGdir/".$cogfiles[$i] unless (-e "$bdir/$1.fna.clstr" || -e "$bdir/$1.fna.ctsv");
		}
	}

	#load marker gene specific cutoffs
	my $dir2cutoffs = getProgPaths($speciesCutoff,0);
	#DEBUG: print "$dir2cutoffs\n\n";open I,"<$dir2cutoffs" or die "Cant open $dir2cutoffs\n";;
	
	my %FMGcutoffs = %{readTabbed3($dir2cutoffs,1)};
	#die;
	if (0){ #DEBUG
		my @kk = keys %FMGcutoffs;
		die "$kk[0]: $FMGcutoffs{$kk[0]} ". @kk ."\n @kk\n";
	}
	my $relaxFMG = 0;#0.005;


	
	$cmd .= "rm -rf $tmpDir\n" unless (-e $clnLnStone && !-e $complStone);
	$cmd .="mkdir -p $tmpDir\n";
	#die "$clnLnStone\n$complStone\n$cmd\n";
	if ($submitLocal && !-e $moveStone){systemW $cmd;$cmd="";}
	
	#cluster FMGs
	my %FMGFL2 ; my $dirflag=0; my $cpFromP = 0;
	my @COGlst = keys %FMGfileList;
	#die "@COGlst\n";
	$cmd .= "mkdir -p $bdir/$COGdir/\n" ;
	$cmd .= "mkdir -p $tmpDir/$COGdir/\n" ;
	my $preclustMMseq = $clustMMseq;
	my $useMMSEQs4COG = 1;
	print "Using mmseqs2 for COG clustering: $useMMSEQs4COG\n";
	foreach my $cog ( @COGlst){
		
		die "can't find $cog in FMGcutoffs list\n" unless (exists $FMGcutoffs{$cog});
		$FMGFL2{$cog} = "$bdir/$COGdir/$cog.$cdhID.fna";
		#print "$FMGFL2{$cog} \n";
		if (!-e $FMGFL2{$cog} || !-s  $FMGFL2{$cog}){#"$bdir/COG/$cog.$cdhID.fna"){
			#$cmd .= $cdhitBin."-est -i $FMGfileList{$cog} -o $tmpDir/COG/$cog.$cdhID.fna -n 9 -G 1 -aS 0.95 -aL 0.6 -d 0 -c ". $FMGcutoffs{$cog}/100 ." -g 0 -T $numCor\n";
			# $clustMMseq = 0; #use mmseq, and use it's slow mode instead..  
			$cmd .= clusterFNA($FMGfileList{$cog},$FMGFL2{$cog},0.9,0.0,($FMGcutoffs{$cog}/100)-$relaxFMG,$numCor3,1,"$NodeTmpDir/$cog/",$useMMSEQs4COG,$totMemL);
		} else {
			$cpFromP = 1;
			#$cmd .= "cp $bdir/COG/$cog.$cdhID.fna* $tmpDir/COG/;"; 
		}
	}
	$clustMMseq = $preclustMMseq;
	foreach my $cog ( @COGlst){
		last; #deactivated.. will be detected in 
		#my ($tmpCmd,$bwtIdxT) =  buildMapperIdx($FMGFL2{$cog},$numCor,1,3);
		#if (!-e $bwtIdxT){$cmd .= $tmpCmd."\n";}
	}
	
	$cmd .= "cp $bdir/$COGdir/*.$cdhID.fna* $tmpDir/$COGdir/\n"; 
	$cmd .= "touch $cogStone\n";
	#die $cmd."\n";
#		$cmd .=  "cp $tmpDir/COG/COG*.$cdhID.fna* $bdir/COG\n";
	$cmd .= "\n";
	my $COGdep="";
	if ($submitLocal ){ #this loop submits marker gene clusterings
		if (-e $FMGstone){
			$cmd = "";
		}elsif (!$cpFromP){
			my @preCons = @{$QSBoptHR->{constraint}};
			push(@{$QSBoptHR->{constraint}}, $avx2Constr) if ($clustMMseq);
			my $preHDDspace = ${$QSBoptHR}{tmpSpace};
			${$QSBoptHR}{tmpSpace} = "50G";#"${totMem}G"; #$totMem #doesn't need much, stores on scrach
			my $memCOG = int($totMemL/2);if ($memCOG < 50){$memCOG=50;}
			my ($dep,$qcmd) = qsubSystem($qsubDir."cogCluster.sh",$cmd,int($numCor/2),($memCOG/$numCor/2)."G","cCLGC","","",1,[],$QSBoptHR);
			@{$QSBoptHR->{constraint}} = @preCons;
			${$QSBoptHR}{tmpSpace} = $preHDDspace;
			$COGdep = $dep;
		} else { systemW $cmd;}
		$cmd = "";
	}

	#die "$stoneDir/complCl.stone";

	
	
	if (!-e $incomplStone) { #map incompletes on complete clusters & merge sams to cdhit format
		if ($allinClust){#default.. just cluster all genes (without marker genes) at 95% nt id
			$cmd .= clusterSingleStep($complStone,$incomplStone,$clnLnStone,$cogStone,$bdir,$OutD,$cmd,$COGdep);
		} else {
			$cmd .= clusterMultiStep($complStone,$incomplStone,$clnLnStone,$cogStone,$bdir,$OutD,$cmd,$COGdep);
		}
	} elsif (!-e $declStone) { #restore files.. for post processing steps
		$cmd .= "zcat $bdir/$primaryClusterFNA.gz > $tmpDir/$primaryClusterFNA;\n" unless (-e "$tmpDir/$primaryClusterFNA");
		$cmd .= "zcat $bdir/$primaryClusterCLS.gz > $tmpDir/$primaryClusterCLS;\n" unless (-e "$tmpDir/$primaryClusterCLS");
		if ($submitLocal){systemW $cmd;$cmd="";}
	}
	qsubSystemJobAlive( [$COGdep],$QSBoptHR ); 
	#from now on all single core jobs..
	#die;
	
	$cmd .= "cp $tmpDir/cluster.ids* $OutD 2>/dev/null || : \n" unless (-e "$OutD/cluster.ids.primary" || -e $matrixSton);
	if ($submitLocal){systemW $cmd; $cmd ="";}
	#remove blank lines..
#	if (!-e $clnLnStone){
#		$cmd .= "\nsed '/^\$/d' $tmpDir/compl.incompl.$cdhID.fna > $tmpDir/compl.incompl.$cdhID.fna.tmp;rm $tmpDir/compl.incompl.$cdhID.fna;mv $tmpDir/compl.incompl.$cdhID.fna.tmp $tmpDir/compl.incompl.$cdhID.fna\n";
#	}

	if (!-e "$GCdir/LOGandSUB/${COGdir}.clusN.log"){
		die "Can;t find $GCdir/LOGandSUB/${COGdir}.clusN.log\nDid merge happen?";
	}
	
		#move to final location
	if (!-e $moveStone){
		$cmd .= "#cp relevant files to outdir and zip the rest\n";
		$cmd .= "mv $tmpDir/${primaryClusterFNA}* $OutD\n" unless (-e "$OutD/$primaryClusterFNA");
		$cmd .= "mkdir -p $qsubDir/${COGdir}\ncp $tmpDir/$COGdir/*.fna.c* $tmpDir/${COGdir}.clusN.log $qsubDir/${COGdir} 2>/dev/null || :\n" if (@COGlst > 0);
		$cmd .= "touch $moveStone\n";
		if ($submitLocal){systemW $cmd;		$cmd = "";}
		die "Can't find $moveStone\n" if ($submitLocal && !-e $moveStone);
		my $cmdL = "$pigzBin -p $numCor $bdir/*\n";
		if ($submitLocal){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep,$qcmd) = qsubSystem($qsubDir."pigz1_GC.sh",$cmdL,$numCor3,int(20/$numCor3)."G","pigzGC","","",1,[],$QSBoptHR); 
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
		} else {$cmd .= $cmdL;}
	}
	if ( !-e $moveStone ){
		die "GC Files not found at final location $OutD; missing $moveStone\n";
	}
	

	
	my @matDeps = ();
	#calc gene matrix.. can as well be run later on finished file..
	#requires genes2row file
	if (!-e $matrixSton && $doGeneMatrix){
		my $numCorL = 4; 
		my $geneMatSupplFlag = " -calcSupplCov "; if (!$CalcgGneMatSuppl){$geneMatSupplFlag="";}
		if ($submitLocal){die"Can;t find $OutD/$primaryClusterCLS" unless (-e "$OutD/$primaryClusterCLS");}
		my $newMapp = ""; $newMapp = "-oldMapStyle" if ($oldNameFolders > 0);
		my $cmd1 = "$rareBin geneMat -i $OutD/$primaryClusterCLS -t $numCorL -o $OutD/$countMatrixP $newMapp -map $mapF -refD $assDirs $geneMatSupplFlag -gz\n"; #add flag -useCoverage to get coverage estimates instead
		#print "$cmd1\n";
		my $cmd2 = "$rareBin geneMat -i $OutD/$primaryClusterCLS -t $numCorL -o $OutD/Mat.cov $newMapp -map $mapF -refD $assDirs $geneMatSupplFlag -useCoverage -gz\n\nrm $OutD/Mat.cov.genes2rows.txt"; #coverage mat
		my $cmd3 = "$rareBin geneMat -i $OutD/$primaryClusterCLS -t $numCorL -o $OutD/Mat.med $newMapp -map $mapF -refD $assDirs $geneMatSupplFlag -useCovMedian -gz\nrm $OutD/Mat.med.genes2rows.txt"; #coverage mat
		$cmd1 .= "\ntouch $matrixSton\n";
		
		if ($submitLocal){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep1,$qcmd1) = qsubSystem($qsubDir."genemat1.sh",$cmd1,$numCorL,int($totMem3/$numCorL)."G","GM1","","",1,[],$QSBoptHR);
			my ($dep2,$qcmd2) = qsubSystem($qsubDir."genemat2.sh",$cmd2,$numCorL,int($totMem3/$numCorL)."G","GM2","","",1,[],$QSBoptHR);
			my ($dep3,$qcmd3) = qsubSystem($qsubDir."genemat3.sh",$cmd3,$numCorL,int($totMem3/$numCorL)."G","GM3","","",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			@matDeps = ($dep1,$dep2,$dep3);
			#qsubSystemJobAlive( [$dep1,$dep2,$dep3],$QSBoptHR ); 
			#systemW $cmd;$cmd = "";
		} else {
			$cmd .= $cmd1 . $cmd2 . $cmd3;
		}

		#die $cmd;
		#$cmd .= "wait\n";
		#$cmd .= "$pigzBin -p $numCor $OutD/Matrix.mat\n";
		#$cmd .= "$pigzBin -p $numCor $OutD/Mat.cov* \n";
		#$cmd .= "$pigzBin -p $numCor $OutD/Mat.med* \n";
		
		#die "Can't find $matrixSton\n" if ($submitLocal && !-e $matrixSton);
	}


	qsubSystemJobAlive( \@matDeps,$QSBoptHR ); 
	die "Matrix were not created, $matrixSton\n" if ($submitLocal && !-e $matrixSton);




	#get protein sequences for each gene & rewrite seq names to numbers
	unless (-e "$OutD/compl.incompl.$cdhID.prot.faa" && -e $protStone){
		$cmd .= "perl $selfScript -mode protExtract -tmp $tmpDir -MGset $useGTDBmg -c $numCor3 -clusterID $cdhID -map \"?\" -o $OutD -extraGenesAA \"$extraRdsFAA\" -oldStyleFolders $oldNameFolders\n";
		$cmd .= "touch $protStone\n";
		#die "$cmd\n\n";
		if ($submitLocal){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep,$qcmd) = qsubSystem($qsubDir."protExtr.sh",$cmd,1,$totMem3."G","protExt","","",1,[],$QSBoptHR); @matDeps = ($dep);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			#systemW $cmd;		
			$cmd = "";
			
		}
	}
	print "Waiting for protein extraction..\n";
	qsubSystemJobAlive( \@matDeps,$QSBoptHR ); 
	die "Prot extraction unsuccessful\n" unless (-e $protStone);
	#die;
	
	#now decluter based on proteins.
	if (@samples > 2 && !-e $declStone && $doDecluter){
		my $localExe=1;$localExe=0 if ($submitLocal);
		$cmd .= "$decluterGC $OutD $tmpDir $numCor $localExe $totMem3 $declStone\n";
		#$cmd .= "touch $declStone\n";
		if ($submitLocal){systemW $cmd;		$cmd = "";}
		#$GCd/decluter/declut.stone
		die "Can't find $declStone\n" if ($submitLocal && !-e $declStone);
	}elsif (!$doDecluter){
		system "touch $declStone" unless (-e $declStone);
	}

	#get 100 marker genes
	my $depExtr = "";
	#single cores...
	if ($submitLocal && !-e $FMGstone){
		$cmd .= "#get marker genes and create matrices for these\n";
		$cmd .= "$extre100Scr $OutD $tmpDir/FMG1/\n";
		$cmd .= "touch $FMGstone\n";
		#systemW $cmd;		$cmd = "";
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		my ($dep,$qcmd) = qsubSystem($qsubDir."GC_ExtrE100.sh",$cmd,1,int($totMem3*2)."G","GCe100","","",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$depExtr = $dep; $cmd = "";
	}
	
	#start LCA for marker genes..
	
	unless (-e $cogLCAStone){
		$cmd .= getProgPaths("MG_LCA_scr") . " -GCd $OutD -tmp $tmpDir/LCA/ -MGset $useGTDBmg -c $numCor ;\n" ;
		$cmd .= "touch $cogLCAStone\n";
		if ($submitLocal && !-e $cogLCAStone){
			#this is now a control script, that submits further jobs..
			print "submitting MG LCA\n";
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep1,$qcmd1) = qsubSystem($qsubDir."MG_LCA.sh",$cmd,1,"10G","MG_LCA",$depExtr,"",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
#			my ($dep1,$qcmd1) = qsubSystem($qsubDir."MG_LCA.sh",$cmd,$numCor,int($totMem3/$numCor)."G","MG_LCA",$depExtr,"",1,[],$QSBoptHR);
			$depExtr = $dep1;$cmd="";
		}
	}
	#die;
	#and get specI's.. the road to SNP genes then
	#some dependencies..
	my $SIdep=""; 
	#deactivated for now, this part is no longer needed..
	if (0 && !-e $SIstone){
		my $siScr = getProgPaths("specIGC_scr");
		$cmd .= "$siScr -GCd $OutD -cores $numCor3 -tmp $tmpDir/SI/ -MGset $useGTDBmg \n\n";

		#$cmd .= "$selfScript -GCd $OutD -m specI -MGset $useGTDBmg -c $numCor3\n"; #-o $OutD 
		$cmd .= "touch $SIstone\n";
		if ($submitLocal && $cmd ne ""){
			print "submitting specI tax abundance..\n";
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep,$qcmd) = qsubSystem($qsubDir."specI_GC.sh",$cmd,1,int(1/$numCor3)."G","SIGC","$depExtr","",1,[],$QSBoptHR);
			$cmd="";$SIdep = $dep;
			$QSBoptHR->{tmpSpace} = $tmpSHDD;
		}
	}
	#and calculate kmer per gene
	if (!-e $geneStatsStone){
		my $cmd1="";my $cmd3="";my $cmd2="";
		$cmd1 = "$kmerScr $OutD $numCor\ngzip $OutD/$primaryClusterFNA.kmer\n" unless (-e "$OutD/$primaryClusterFNA.kmer.gz");
		$cmd2 = "$GCcalc $OutD/$primaryClusterFNA $OutD/compl.incompl.$cdhID.fna.GC \n" unless (-e "$OutD/compl.incompl.$cdhID.fna.GC");
		$cmd3 = "$genelengthScript $OutD/$primaryClusterFNA $OutD/$primaryClusterFNA.length \n" unless (-e "$OutD/$primaryClusterFNA.length");
		if ($submitLocal){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep,$qcmd) = qsubSystem($qsubDir."GC_GC.sh",$cmd2,1,int(30/1)."G","GCGC","","",1,[],$QSBoptHR);
			my $dep2="";
			#my ($dep2,$qcmd2) = qsubSystem($qsubDir."kmer_GC.sh",$cmd1,1,int($totMem/1)."G","kmGC","","",1,[],$QSBoptHR);
			my ($dep3,$qcmd3) = qsubSystem($qsubDir."len_GC.sh",$cmd3,1,int(30/1)."G","leGC","","",1,[],$QSBoptHR);
			($dep,$qcmd) = qsubSystem($qsubDir."CheckGstats.sh","touch $geneStatsStone\n",1,"1G","GeStats","$dep;$dep2;$dep3","",1,[],$QSBoptHR); 
			$QSBoptHR->{tmpSpace} =$tmpSHDD;

		} else {
			$cmd .= $cmd1.$cmd2.$cmd3;
			$cmd .= "touch $geneStatsStone\n";
		}
		#do not check since this is just running in the background till exhaustion..
		#die "Can't find $geneStatsStone\n" if ($submitLocal && !-e $geneStatsStone);
	}
	



	#MAG related..
	$cmd .= "#taxonomic assignments of all genes via kraken\n";
	$cmd .= "$selfScript -mode kraken -MGset $useGTDBmg -o $OutD -c $numCor3\n";
	$cmd .= "touch $krakStone\n";
	if (-e $krakStone){$cmd="";}
	if ($submitLocal && $cmd ne ""){
		print "submitting kraken tax abundance..\n";
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		my ($dep,$qcmd) = qsubSystem($qsubDir."krak_GC.sh",$cmd,$numCor3,int($totMem5/$numCor3)."G","krakGC","","",1,[],$QSBoptHR); $cmd="";
		$QSBoptHR->{tmpSpace} =$tmpSHDD;

	}
	
	#functional annotations.. just run some by default
	$cmd .= "#functional assignments of all genes via diamond\n";
	$cmd .= "$selfScript -mode FuncAssign -MGset $useGTDBmg -o $OutD -c $numCor3 -functDB $curDB_o -functAligner $funcAligner \n";
	$cmd .= "touch $funcStone\n";
	if (-e $funcStone){$cmd="";}
	if ($submitLocal && $cmd ne ""){
		print "submitting diamond func abundance..\n";
		my ($dep,$qcmd) = qsubSystem($qsubDir."func_GC.sh",$cmd,1,int($totMem3)."G","funcGC","","",1,[],$QSBoptHR); $cmd="";
		#systemW $cmd; $cmd = "";
	}


	$cmd .= "#functional assignments via eggNOGmapper\n";
	#-c $numCor3 .. use max 6 cores for this due to single core emapper final step
	$cmd .= "$selfScript -mode FuncEMAP -MGset $useGTDBmg -o $OutD -c 6 -stone $emapStone \n";
	#$cmd .= "touch $emapStone\n";
	if (-e $emapStone){$cmd="";}
	if ($submitLocal && $cmd ne ""){
		print "submitting eggNOGmapper func abundance..\n";
		my ($dep,$qcmd) = qsubSystem($qsubDir."emap_GC.sh",$cmd,1,int($totMem3)."G","emapGC","","",1,[],$QSBoptHR); $cmd="";
	}
	
	
	
	#die;

	#$cmd .= "$selfScript -mode kaiju -MGset $useGTDBmg -o $OutD\n"; #kaiju is too instable.. don't use
	
	#$cmd .= cleanUpGC($bdir,$OutD,$cdhID);
	$cmd .= "#Canopy clustering\n";
	$cmd .= "#" unless ($doMags);
	$cmd .= "$GCscr -mode CANOPY -o $OutD -c $numCor -tmp $tmpDir/MGS\n";
	$cmd .= "touch $canopyStone\n";
	if (-e $canopyStone || $numSmpls < 10){$cmd="";}
	my $CANdep="";
	if ($submitLocal && $cmd ne "" ){
		print "submitting canopy clustering..\n";
		$CANdep = canopyCluster($GCdir,"$tmpDir/cano/",$numCor,$canopyStone);
		$cmd="";
		#this needs to call the canopy function here, in order to get the correct dependency for the canopy call..
		#my ($dep,$qcmd) = qsubSystem($qsubDir."canopy_GC.sh",$cmd,$numCor,int($totMem3/$numCor)."G","canGC","","",1,[],$QSBoptHR); $cmd="";
		#$CANdep = $dep;
	}
	$cmd .= "\ntouch $preMGSstone\n";
	my $canopyExpectedDir = getCanopyDir($OutD);
	#$cmd .= "#MetaBat2 single sample (sample group) MAGs\n";
	#cross compare MB2, extract more genes via canopy, fix via correlation stats
	#is dependent on SNPs being called..
	if ($doStrains){
		$cmd .= "#run from job distributing node:\n";
		$cmd .= "#" if (!$submitLocal);
		#$cmd .= "$specIphyloScript $OutD $numCor"; #I deactivated this script, as same as strainWithin.pl, just without bins.. so useless
	}
	$cmd .= "\n\n";
	#$cmd .= "rm -f $bdir/SAM/compl.$cdhID.fna.b* $bdir/SAM/compl.$cdhID.fna$mini2IdxFileSuffix\n";
	$cmd .= "echo \"cleaning up B0 tmp dir\"\nrm -rf $bdir\n"; #not really needed any longer
	#die $cmd;
	#single cores...
	if ($submitLocal){systemW $cmd;		$cmd = "";}
	$cmd .= "\n\n#next step - MAG creation\n#MGS and MAG creation - can take awhile...\n";
	$cmd .= "#" unless ($doMags);
	#mem before was $totMem, but this is more a limiatation of checkM that should be used..
	
	#estimate ram size..
	my $MGSoutD = "$OutD/Bin_SB/";
	#if ($submitLocal){systemW $cmd;		$cmd = "";}
	if (-e $qsubDir."MGS.sh"){
		my $tmpS = `grep -v '^#' $qsubDir/MGS.sh`;
		if ($tmpS =~ m/-outD\s+(\S+)/){
			$MGSoutD = $1;
		}
	}
	my $canoStr = "-canopies $canopyExpectedDir/clusters.txt ";
	$canoStr = "" if ($numSmpls < 10);
	#$cmd .= "#wait for eggnogmapper to finish\nuntil [ -f $emapStone ];do sleep 5; done\n" unless (-e $emapStone);
	$cmd .= "$magPi -mem 150 -GCd $OutD -tmp $tmpDir/MAGs/ -bottleneckCores $numCor $canoStr -strains $doStrains -useCheckM2 $useCheckM2 -useCheckM1 $useCheckM1 -wait4stone $emapStone -binSpeciesMG $binSpeciesMG -ignoreIncompleteMAGs $ignoreIncompleteMAGs -MGset $useGTDBmg -outD $MGSoutD \n";

	print $cmd."\n\n";
	

	

	if ($submitLocal && $cmd ne ""){
		print "submitting MGS script..\n";
		my $idxFileSize = -s "$OutD/compl.incompl.95.fna.clstr.idx"; $idxFileSize /= (1024*1024 * 1024); #size in kB->MB->GB
		$QSBoptHR->{useLongQueue} = 1;
		#is a single core script, should be treated and submitted differently from gene cat script 
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		
		my ($dep,$qcmd) = qsubSystem($qsubDir."MGS.sh",$cmd,1,int($idxFileSize*12+15)."G","MGSgc",$CANdep.";".$depExtr,"",1,[],$QSBoptHR); $cmd="";
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;

	}

	#$cmd .= "$magPi $GCd $tmpDir/MAGs $numCor\n";
	#add gc calcs
	#these incompletes have to be added to big gene catalog in last step by themselves (avoid mis-center clustering

	#tabix prep of matrix
	#tabix doesn't work properly.. use sed idx'ing instead
	#$cmd.= "$bgzipBin $OutD/Matrix.mat\n" ;
	#$cmd.= "$tabixBin -S 1 -s 1 $OutD/Matrix.mat.gz\n";
	#die $cmd."\n";
	$cmd .= "rm -f -r $tmpDir\n";
	
	$cmd .= "\n\necho \"=======================================\"\necho \"Finished Gene Catalog script\"\necho \"=======================================\"\n\n";
	my $jobName = "CD_$nm";
	
	#die "$cmd\n";


	if (!$submitLocal){
		print $qsubDir."CDHITexe.sh\n";
		my $jobd = "";
		$QSBoptHR->{useLongQueue} = 1;
		my ($dep,$qcmd) = qsubSystem($qsubDir."CDHITexe.sh",$cmd,$numCor,int($totMem/$numCor)."G",$jobName,$jobd,"",$qsubNow,[],$QSBoptHR);
		if ($qsubNow==0){
			print "$qcmd\n";
		}
	}
	#die ($cmd);
	
}

sub ntMatchGC{
	my ($GCdir,$refDB,$out) = @_;
	die "No -refDB given\n" if ($refDB eq "");
	die "No -GCd given\n" if ($GCdir eq "");
	die "No -out given\n" if ($out eq "");
	my $pctID = 0.1; my $pctCov=0.75; my $mapQual=15;
	my $smtBin = getProgPaths("samtools");
	my $bamfilter = getProgPaths("bamFilter_scr");
	my $geneFNA = "$GCdir/$primaryClusterFNA";
	$out =~ m/^(.*\/)[^\/]*$/;my $outD=$1;
	system "mkdir -p $outD" unless (-d $outD);
	my $iTO = "$out.bam";
	my $iTO2 = "$out.txt";
	my ($tmpCmd,$bwtIdxT) =  buildMapperIdx($refDB,$numCor,1,3);
	my $cmd = "";
	$cmd .= "$tmpCmd\n";
	my $mini2Base = "$mini2Bin -2 -a -t $numCor --secondary=no -Y -x asm20 "; #--sam-hit-only
	$cmd .= "$mini2Base $bwtIdxT $geneFNA | $bamfilter $pctID $pctCov $mapQual | $smtBin view -b1 -@ $numCor -F 4 - > $iTO\n";
	#die "$cmd\n\n";
	#convert bam to txt file
	$cmd .= "$smtBin view $iTO | cut -f1,3 > $iTO2 \n";

	$cmd .= "$rareBin sumMat -i $GCdir/$countMatrixF.gz -o $out.mat -t $numCor -refD $iTO2 \n"; #$rtkFunDelims
	my ($dep,$qcmd) = qsubSystem($qsubDir."NucMap.sh",$cmd,$numCor,int($totMem5/$numCor)."G","IDXmGC","","",1,[],$QSBoptHR);
	#die $cmd;
}

sub addingSmpls{
	my ($startSmpl,$stopSmpl,$batch) = @_;
	print "Adding up Smpls $startSmpl:$stopSmpl\n" ;#if ($justCDhit==0);
	my @skippedSmpls;my $wrongSmplNms = ""; my @rmSrcDirs;
	my %uniqueSampleNames; 	my $doubleSmplWarnString = ""; 
	my @OCOMPL = (); my @O3P=(); my @O5P = (); my @OINC = (); #these arrays store complete & incomplete fasta seqs

	open QLOG,">$qsubDir/GeneCompleteness.txt.$batch";
	#print QLOG "Smpl\tComplete\t3'_compl\t5'_compl\tIncomplete\tTotalGenes\n";

	my %allFMGs; #stores FMGs, seperated by COG subsets 
	my $OFcompl = "$tmpDir/compl.fna.gz";my $OFincompl = "$tmpDir/incompl.fna.gz";my $OF5in = "$tmpDir/5Pcompl.fna.gz";my $OF3in = "$tmpDir/3Pcompl.fna.gz";
	my $OFcomplF = "$bdir/compl.fna.gz";my $OFincomplF = "$bdir/incompl.fna.gz";my $OF5inF = "$bdir/5Pcompl.fna.gz";my $OF3inF = "$bdir/3Pcompl.fna.gz";
	#prep system
	system "mkdir -p $tmpDir" unless (-d $tmpDir);
	my $OC; my $O5; my $O3; my $OI;
	#die "shouldn't \n";
	open $OC,"| gzip -f -3 -c  >$OFcompl.$batch" or die "Can't open $OFcompl.$batch\n"; #close $OC;
	open $O5,"| gzip -f -3 -c > $OF5in.$batch" or die "Can't open $OF5in.$batch\n"; #close $O5;
	open $O3,"| gzip -f -3 -c > $OF3in.$batch" or die "Can't open $OF3in.$batch\n"; #close $O3;
	open $OI,"| gzip -f -3 -c > $OFincompl.$batch" or die "Can't open $OFincompl.$batch\n"; #close $OI;
	print "Direct output to file\n";

	#now really add all files together
	my $JNUM= -1; my $stopAndRedo=0;
	my @missedSmpls=();
	#foreach my $smpl(@samples){
	for (my $JNUM=$startSmpl; $JNUM<$stopSmpl ; $JNUM++){
		next if ($JNUM >= @samples);
		my $smpl = $samples[$JNUM];
		#last if ($justCDhit==1);
		my $dir2rd = $map{$smpl}{wrdir};
		$dir2rd = $map{$smpl}{prefix} if ($dir2rd eq "");
		#die "\n".$smpl."X\n";
		if ($map{$smpl}{ExcludeAssem} eq "1"){next;}
		if ($dir2rd eq "" ){#very specific read dir..
			if ($map{$smpl}{SupportReads} ne ""){
				$dir2rd = "$GCdir$smpl/";	
			}elsif ($requireAllAssemblies){
				print "Can't find valid path for $smpl, continue without\n";push(@missedSmpls, $smpl);next;
			} else {
				die "Can;t find valid path for $smpl\n";
			}
		} 
		
		#check if mult assembly and adapt
		my $assGo = 0;
		my $cAssGrp = $map{$smpl}{AssGroup};
		$AsGrps{$cAssGrp}{CntAss} ++;	
		print "$JNUM - $smpl - $cAssGrp - $AsGrps{$cAssGrp}{CntAss}:$AsGrps{$cAssGrp}{CntAimAss}\n";
		if ($AsGrps{$cAssGrp}{CntAss}  >= $AsGrps{$cAssGrp}{CntAimAss} || 
				$map{$smpl}{assFinSmpl} eq $smpl){ 
			$assGo = 1;
		}
		unless ($assGo){print "Not last in comb assembly: ".$map{$smpl}{dir}."\n"; next;}
		my $SmplName = $map{$smpl}{SmplID};
		#$dir2rd = "/g/scb/bork/hildebra/SNP/SimuL/sample-0/";
		my $metaGD = "$dir2rd/assemblies/metag/";
		if (!-e "$metaGD/longReads.fasta.filt.sto" && !-e "$metaGD/scaffolds.fasta.filt"){$metaGD = `cat $dir2rd/assemblies/metag/assembly.txt`; chomp $metaGD;}
		my $inFMGd = "$metaGD/ContigStats/FMG/";
		#print "\n$metaGD\n";
		#print "$dir2rd/assemblies/metag/scaffolds.fasta.filt\n";
		if ((! fileGZe("$metaGD/scaffolds.fasta.filt") || ! fileGZe("$metaGD/longReads.fasta.filt")) && ! fileGZe("$metaGD/$path2nt")){# && -d $inFMGd){
			print "Skipping $dir2rd\n";
			die "no ass1\n $metaGD\n$dir2rd/assemblies/metag/assembly.txt\n" unless (-e "$metaGD/scaffolds.fasta.filt" && !-e "$metaGD/longReads.fasta.filt" ); 
			die "no ass2\n $$metaGD/longReads.fasta.filt\n" unless (-e "$metaGD/longReads.fasta.filt" ); 
			die "no NT\n" unless (-e "$metaGD/$path2nt");
			push(@skippedSmpls,$map{$smpl}{dir});
			next;
		}
		#next;
		print "==== ".$dir2rd." ====\n";
		#print LOG "==== ".$dir2rd." ====\n";
		my $inGenesF = "$metaGD/$path2nt";
		#my $inGenesFs = $inGenesF; $inGenesFs =~ s/\.fna$//;
		my $fnaHref = readFasta($inGenesF);
		my %fnas = %{$fnaHref}; my @scnts = (0,0,0,0,0);
		my $gffHref= readGFF("$metaGD/$path2gff");
		my %gff = %{$gffHref};
		my %curFMGs; #FMGs and their ID
		my %curFMGsTag; #checks that all FMGs were present in fasta.. can point to corrupted files
		if ($doFMGseparation){
			open I, "<$metaGD/$path2FMGids" or die "cant open FMGids:\n$metaGD/$path2FMGids\n";
			my $cnt = 0;
			while (my $line = <I>){
				chomp $line; my @spl = split(/\s+/,$line); #MM1__C104459_L=563;_1 COG0552
				if ($cnt ==0 ){ 
					die "wrong FMGid: $line\n" unless ($spl[0] =~ m/(^.*)__C/);
					#check that samples was only used once
					if (!$oldNameFolders && $SmplName ne $1 && !exists($map{altNms}{$1})){
						my $errStr = "expected Sample name $SmplName, found $1\n";
						$wrongSmplNms .= $errStr;
						push (@rmSrcDirs, $dir2rd,$metaGD) ;
					}
					if (exists($uniqueSampleNames{$1})){
						my $errStr = "$dir2rd: Can't use sample names twice: $1\n$spl[0]\ndblTo:$uniqueSampleNames{$1}\n";
						print STDERR $errStr;
						$doubleSmplWarnString .= $errStr;
					} else { $uniqueSampleNames{$1} = $dir2rd;}
					
				}
				die "fasta headers are not in required format:$spl[0]\n$metaGD/$path2FMGids\n" unless ($spl[0] =~ m/.*__.*L=\d+=_\d/);
				$curFMGs{">".$spl[0]} = $spl[1];
				$curFMGsTag{">".$spl[0]} = 0;
				#die "$spl[0]\n";
				$cnt++;
			}	
			close I;
		}
		#split into buckets
		my $tooShrtCnt=0; my $prevSmpID = "";
		foreach my $hd (keys %fnas){
			if (length($fnas{$hd}) <= $minGeneL && !exists $curFMGs{$hd}){$tooShrtCnt++;next;}
			my $shrtHd = $hd ;	#$shrtHd =~ m/(\S+)\s/; $shrtHd = $1;
			#print $shrtHd."\n$hd\n";			print "$fnas{$hd}\n";
			die "fasta headers are not in required format:$shrtHd\n$inGenesF\n" unless ($shrtHd =~ m/.*__.*L=\d+=_\d/);
			my @spl = split /__/,$shrtHd;
			if ($spl[0] ne $prevSmpID){if ($prevSmpID eq "") {$prevSmpID = $spl[0];} else {die "Mix of several samples?? $shrtHd, $spl[0] detected, expected sample $prevSmpID !! Aborting\n\n";} }
			unless (exists $gff{$shrtHd}){
				system "rm -r $metaGD/$path2GPdir $metaGD/$path2CS"; $stopAndRedo=1; 
				print "can't find gff entry for $shrtHd\nDeleting entire gene prediction and coverage..";next;
			}
			unless ($gff{$shrtHd} =~ m/;partial=(\d)(\d);/){ die "Incorrect gene format for gene $hd \n in file $inGenesF\n";}
			if (exists $curFMGs{$hd} && $doFMGseparation){
				$allFMGs{$curFMGs{$hd}}{$hd} = $fnas{$hd}; $scnts[4] ++;
				$curFMGsTag{$hd} = 1;
			} elsif ($1==0 && $2==0){ #complete genes
				print $OC $shrtHd."\n".$fnas{$hd}."\n";
				$scnts[0] ++;
			} elsif ($1==0 && $2==1){ #3' complete
				print $O3 $shrtHd."\n".$fnas{$hd}."\n";
				$scnts[1] ++;
			} elsif ($1==1 && $2==0){ #5' complete
				print $O5 $shrtHd."\n".$fnas{$hd}."\n";
				$scnts[2] ++;
			} else { #gene fragments, just map
				print $OI $shrtHd."\n".$fnas{$hd}."\n";
				$scnts[3] ++;
			}
		}
		
		
		my $totCnt = $scnts[0] + $scnts[1] + $scnts[2] + $scnts[3] ;
		my $ostr= "SmplID::$SmplName:$prevSmpID : ";
		if ($totCnt>0){
			$ostr .= $totCnt."+$scnts[4] genes: ".sprintf("%.3f",($scnts[0]/$totCnt*100)). "% complete, ";
			$ostr.=sprintf("%.1f",$scnts[1]/$totCnt*100)."% 3' compl, ";
			$ostr.=sprintf("%.1f",$scnts[2]/$totCnt*100). "% 5' compl, ";
			$ostr.=sprintf("%.1f",$scnts[3]/$totCnt*100). "% incompl, ";
			$ostr.=sprintf("%.1f",$tooShrtCnt). " < $minGeneL";
			
		} else {
			$ostr = "0 genes found.\n";
		}
		
		if ($doFMGseparation){
			foreach my $k (keys %curFMGsTag){
				if ($curFMGsTag{$k} == 0){
					print "FATAL:: missing gene $k\n";
					$stopAndRedo=1;
				}
			}
		}
		print $ostr."\n"; 
		print QLOG "$SmplName\t$prevSmpID\t$scnts[0]\t$scnts[1]\t$scnts[2]\t$scnts[3]\t$totCnt\t$tooShrtCnt\n";
		$cnt++;
		if ( 0 && $cnt % 10 == 0) {#write out & submit cdhit job #not used any longer!
			my $bdir = $GCdir."B$bucketCnt/";
			writeBucket(\@OCOMPL,\@O3P,\@O5P,\@OINC,$bdir,$bucketCnt);
			$bucketCnt++;push(@bucketDirs,$bdir);
			@OCOMPL=();@O3P=();@O5P=();@OINC=();#clean old seqs
		}
	}
	
	if (@missedSmpls){
		print "The following samples were without assembly/gene predictions:\n@missedSmpls\n";
		open SMR,">$qsubDir/Missed_samples.txt";
		print SMR join("\n",@missedSmpls);
		close SMR;

	}
	if ($doubleSmplWarnString ne "" || $wrongSmplNms ne ""){
		print "Recommended to remove: \n'rm -r @rmSrcDirs'\n\n" if (@rmSrcDirs > 0);
		print $doubleSmplWarnString."\n\n\n$wrongSmplNms\n\n\n";
		print "incomplete; aborting process\n";
		exit(21);
	}
	if ($stopAndRedo){
		print "Something wrong while extracting metagenomic genes.. you will need to rerun MATAFILER, fatal assemblies have already been deleted and will now need to be re-assembled.\n";
		exit(31);
	}

	#any extra reads (e.g. from ref genomes?)


	if ($justCDhit==0 && $extraRdsFNA ne ""){
		my $fnaHref = readFasta($extraRdsFNA);
		my %fnas = %{$fnaHref}; 
		my $xcnts = 0; my $tooShrtCnt=0;
		foreach my $hd (keys %fnas){
			if (length($fnas{$hd}) <= $minGeneL){$tooShrtCnt++;next;}
			my $shrtHd = $hd ;	$shrtHd =~ m/(\S+)\s/; $shrtHd = $1;
			#
			#just assume that every gene is complete
			print $OC $shrtHd."\n".$fnas{$hd}."\n";
			$xcnts ++;
		}
		print "Added $xcnts genes from external source\nSkipped $tooShrtCnt Genes (too short $minGeneL)\n";
		print QLOG "Added $xcnts genes from external source\n";
	}
	close QLOG;

	print "\n\n--skipped: ".join(",",@skippedSmpls)."\n" if (@skippedSmpls > 0);

	#die();
	#writeBucket(\@OCOMPL,\@O3P,\@O5P,\@OINC,$bdir,$bucketCnt);
	#write marker genes separate
	foreach my $cog (keys (%allFMGs)){
		system "mkdir -p $bdir/$COGdir/" unless (-d "$bdir/$COGdir/");
		my %cogFMG = %{$allFMGs{$cog}};
		my $ccogf = "$tmpDir/$COGdir/preclus.$cog.fna";
		open Ox,">$ccogf.$batch" or die "Can't open COG output file $ccogf.$batch\n";
		#$FMGfileList{$cog} =  "$ccogf";
		foreach my $geK (sort { length($cogFMG{$b}) <=> length($cogFMG{$a}) } keys %cogFMG) {
			print Ox $geK."\n".$cogFMG{$geK}."\n";
		}
		close Ox;
	}

#	die;
	close $O3;close $OC;close $O5; close $OI;
	sleep (2); #give IO enough time
	
	#already in this process start appending files.. faster than waiting..
	print "###\n###\n###\nData collection finished.. adding to main files\n###\n";
	my @transferFiles = ($OFcompl,$OFincompl,$OF5in,$OF3in);
	my @destFiles = ($OFcomplF,$OFincomplF,$OF5inF,$OF3inF);
	
	for (my $i=0;$i<@transferFiles;$i++){
		my $curTransfer = $transferFiles[$i];
		my $dest = $destFiles[$i];
		my $lockFile = $curTransfer.".lock";
		while (-e $lockFile){sleep(10);}
		if (!-e $lockFile){systemW "touch $lockFile;" ;
		sleep(2);
		systemW "cat $curTransfer.$batch >> $dest;";
		systemW "rm -f  $lockFile $curTransfer.$batch;";
		sleep(2);
		} else {die "lock already existed while attempting to write!!\n\n";}
	}
}

sub prepCDhit(){
	
	return if (-e $prepStone && $justCDhit );
	system "mkdir -p $GCdir" unless (-d $GCdir);
	system "mkdir -p $qsubDir" unless (-d $qsubDir);

	#copy maps
	my @maps = split(/,/,$mapF);
	my @newMaps;my $cntMaps=0;
	foreach my $mm (@maps){
		#system "cp $mm $qsubDir/map.$cntMaps.txt"; 
		system "envsubst < $mm  > $qsubDir/map.$cntMaps.txt";
		push (@newMaps,"$qsubDir/map.$cntMaps.txt"); $cntMaps++;
	}
	open O,">$qsubDir/GCmaps.inf"; print O join ",",@newMaps; close O;
	open O,">$qsubDir/GCmaps.ori"; print O join ",",@maps; close O;
	$mapF = join ",",@newMaps; #always work with copied versions of maps
	#die();

	 my @probSample=();
	#detecting and removing wrong / corrupted samples..
	#first check if all input is present
	my $JNUM= -1;
	foreach my $smpl(@samples){
		$JNUM++;
		last if ($justCDhit==1 || $skipPreCheck);
		print "Checking if all requires input files are present..\n" if ($JNUM==0);
		my $dir2rd = $map{$smpl}{wrdir};
		$dir2rd = $map{$smpl}{prefix} if ($dir2rd eq "");
		if ($map{$smpl}{ExcludeAssem} eq "1"){next;}
		if ($dir2rd eq "" ){#very specific read dir..
			if ($map{$smpl}{SupportReads} ne ""){
				$dir2rd = "$GCdir$smpl/";	
			} else {
				die "Can;t find valid path for $smpl\n";
			}
		} 
		my $assGo = 0;
		my $cAssGrp = $map{$smpl}{AssGroup};
		$AsGrps{$cAssGrp}{CntAss} ++;	
		my $metaGD  = getAssemblPath($map{$smpl}{wrdir},"",0);#"$dir2rd/assemblies/metag/";
		if ($metaGD eq ""){ #something wrong..
			print "No assmebly file: $map{$smpl}{wrdir}\n" ;
			push @probSample, $smpl;
			next;
		}
		#my $inFMGd = "$metaGD/ContigStats/FMG/";
		if ($AsGrps{$cAssGrp}{CntAss}  >= $AsGrps{$cAssGrp}{CntAimAss} ){ $assGo = 1; $AsGrps{$cAssGrp}{CntAss}=0;}
		unless ($assGo){ next;}#print "Not last in comb assembly: ".$map{$smpl}{dir}."\n";
		if ((!-e "$metaGD/scaffolds.fasta.filt" || !-e "$metaGD/longReads.fasta.filt") && !-e "$metaGD/$path2nt"){# && -d $inFMGd){
			print "Skipping $dir2rd\n";
			die "no ass1\n $metaGD\n$dir2rd/assemblies/metag/assembly.txt\n" unless (-e "$metaGD/scaffolds.fasta.filt" && !-e "$metaGD/longReads.fasta.filt" ); 
			die "no ass2\n $metaGD/longReads.fasta.filt\n" unless (-e "$metaGD/longReads.fasta.filt" ); 
			die "no NT\n" unless (-e "$metaGD/$path2nt");
			#push(@skippedSmpls,$map{$smpl}{dir});
			next;
		}
		my $inGenesF = "$metaGD/$path2nt";
		my $inGenesGFF = "$metaGD/$path2gff";
		my $problem = 0;
		if (!-e $inGenesF){
			print "Gene predictions not present: $inGenesF\n" ;
			$problem = 1;
		}
		if (!-e $inGenesGFF){
			print "Gene annotations not present: $inGenesGFF\n" ;
			$problem = 1;
		}
		my $FMGf = "$metaGD/$path2FMGids";#"$inFMGd/FMGids.txt";
		if ($doFMGseparation && !-e $FMGf ){
			print "No FMG ids: $FMGf\n$JNUM : $dir2rd\n" ;
			$problem = 1;
		}
		if ($problem ){
			push @probSample, $smpl;
		}

	}
	if (@probSample){
		print "Found problematic samples:\n@probSample\nPlease check mapping & contig stats\n";
		exit (23) if ($requireAllAssemblies);
	} else {	print "All required input files seem to be presents.\n"; }

	if ($justCDhit == 0 || !-e $prepStone){ #recreate input files..
	
		my $maxSmpls = scalar(@samples);
		print "Preparing splitting preprocessing of $maxSmpls metagenomes in $batchNum batches.\n";
		my $batch = 0; my @jobs;
		system "mkdir -p $tmpDir/$COGdir" unless (-d "$tmpDir/$COGdir");
		system "mkdir -p $qsubDir/preprocess/" unless (-d "$qsubDir/preprocess/");
		system "rm -f touch $prepStone.1";
		for ( $batch = 0; $batch < $batchNum;$batch ++){
			my $locTo = int($maxSmpls/$batchNum*(1+$batch));
			my $locFrom = int($maxSmpls/$batchNum*($batch));
			$locFrom ++ if (($batch +1) == $batchNum); ## last sample.. just add one extra to be safe
			$locFrom = 0 if ($batch == 0);
			#print "$locFrom,$locTo\n";
			  
			my $cmd = "$selfScript -mode subprepSmpls -GCd $GCdir -map $mapF -tmp $tmpDir -SmplStart $locFrom -SmplStop $locTo -SmplBatch $batch";
			#systemW $cmd."\n";
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my $numCor = 3;
			my ($jdep,$txtBSUB) = qsubSystem($qsubDir."/preprocess/Preprocess.$batch.sh",$cmd,$numCor,int(30/$numCor)."G","PrPr$batch","","",1,[],$QSBoptHR);
			push(@jobs,$jdep);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;

			# addingSmpls($locFrom,$locTo,$batch);  subprepSmpls
		}
		my $cmd2 = "touch $prepStone.1";
		my ($jdep,$txtBSUB) = qsubSystem($qsubDir."Preprocess.check.sh",$cmd2,1,"1G","CheckPrPr",join(";",@jobs),"",1,[],$QSBoptHR);
		qsubSystemJobAlive( [@jobs,$jdep],$QSBoptHR ); 
		die "Something went wrong in the sample prep..\nCheck $qsubDir/Preprocess.*.sh\n" unless (-e "$prepStone.1");

		#die;
		
		my $OFcompl = "$bdir/compl.fna.gz";my $OFincompl = "$bdir/incompl.fna.gz";my $OF5in = "$bdir/5Pcompl.fna.gz";my $OF3in = "$bdir/3Pcompl.fna.gz";
		if (1){
			#this is now down within the "addingSmpls()" function..
		}elsif ($batch == 0 ){ #only 1 round, simple move
			#system "mv $OFcompl.0 $OFcompl; mv $OFincompl.0 OFincompl; mv $OF5in.0 $OF5in; mv $OF3in.0 $OF3in;";
		} else {
			#print "Concatenating main gene files\n";systemW "cat $OFcompl.* > $OFcompl; rm -f $OFcompl.*;";systemW "cat $OFincompl.* > $OFincompl;rm -f $OFincompl.*;";systemW "cat $OF5in.* > $OF5in;rm -f $OF5in.*";systemW" cat $OF3in.* > $OF3in;rm -f $OF3in.*";
		}
				
		if (-d "$tmpDir/$COGdir") { #COGs were created
			sleep(5);
			print "Concatenating FMG/GTDB gene files\n";
			opendir(DIR, "$tmpDir/$COGdir/") or die $!;
			my @cogfiles = grep {/^preclus\..*\.fna\.\d+/ && -f "$tmpDir/$COGdir/$_" } readdir(DIR); close DIR;
			for (my $i=0; $i<@cogfiles; $i++){
				my $outfile = $cogfiles[$i];
				$outfile =~ s/\.\d+$//;
				next if (-e "$bdir/$COGdir/$outfile");
				#print "$bdir/$COGdir/$outfile\n";
				systemW "cat $tmpDir/$COGdir/$outfile.* > $bdir/$COGdir/$outfile; rm -f $tmpDir/$COGdir/$outfile.* ";
			}
		}
		
		open QLOG,">$qsubDir/GeneCompleteness.txt";
		print QLOG "Smpl\tinternalID\tComplete\t3'_compl\t5'_compl\tIncomplete\tTotalGenes\tGenesTooShort\n";
		close QLOG;
		systemW "cat $qsubDir/GeneCompleteness.txt.* >> $qsubDir/GeneCompleteness.txt; rm $qsubDir/GeneCompleteness.txt.*";
		print "Done concatenating\n";
		
		system "touch $prepStone" unless (-e $prepStone);
	}

}

sub krakenTax{
	my ($GCd,$tmpD,$NC) = @_;
	#my $krkBin = getProgPaths("kraken2");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
	my $krk2Bin = getProgPaths("kraken2");
	my $oriKrakDir = getProgPaths("Kraken2_path_DB");
	my $ete3taxid = getProgPaths("taxid2tax_scr");
	my $miniDB = getProgPaths("Kraken2_mini");
	my $outD = $GCd."/Anno/Tax/";
	system "mkdir -p $outD" unless (-d $outD);
	system "mkdir -p $tmpD" unless (-d $tmpD);
	my @thrs = (0.01,0.02,0.04,0.06,0.1,0.2,0.3);
	my $geneFNA = "$GCd/$primaryClusterFNA";
	#my $curDB = "$oriKrakDir/minikraken_2015";
	my $curDB = "$oriKrakDir/$miniDB";
	#die $curDB."\n";
	#kraken1
	#my $cmd .= "$krkBin --preload --threads $NC --fasta-input  --db $curDB  $geneFNA >$tmpD/rawKrak.out\n";
	#kraken2
	#--use-names --use-mpa-style --report $tmpD/rawKrak.rep
	my $cmd = "";
	my $krak1= 1;	if ($krk2Bin ne ""){ $krak1=0; }
	#only configured for kraken2
	$cmd .= "$krk2Bin --threads $NC --db $curDB $geneFNA  --confidence 0.0  | grep '^C' | cut -f2,3 > $outD/krak2.out\n";

	#die "$cmd\n";
	#for (my $j=0;$j< @thrs;$j++){
	#	$cmd .= "$krkBin-filter --db $curDB  --threshold $thrs[$j] $tmpD/rawKrak.out | $krkBin-translate --mpa-format --db $curDB > $outD/krak_$thrs[$j]".".out\n";
	#}
	print "Starting kraken assignments of the gene catalog\n";
	systemW $cmd unless (-e "$outD/krak2.out");
	my $krakTaxFile = "$outD/krak2.txt";
	if (1 || !-e $krakTaxFile){
		my %allTaxs; my %gene2tax;
		open I,"<$outD/krak2.out" or die "Can't open kraken2 output\n$outD/krak2.out\n";
		while (<I>){chomp;my @spl=split/\t/;$gene2tax{$spl[0]} = $spl[1]; $allTaxs{$spl[1]}=1;}
		close I;
		my @allTs = keys %allTaxs;
		my @taxLvl = qw( d__ p__ c__ o__ f__ g__ s__);
		$cmd = " $ete3taxid ". join(" ",@allTs);#source activate Python27;  ."; source deactivate Python27;";
		#die "$ete3taxid\n";
		my $strings = `$cmd`; @allTs = split /\n/,$strings; $strings="";
		print "All kraken assignments are done\n";
		for (my $x=0;$x<@allTs;$x++){
			my @spl = split /\t/,$allTs[$x]; my $tID = shift @spl;
			for (my $i=0;$i<@spl;$i++){$spl[$i] = $taxLvl[$i].$spl[$i];}
			$allTaxs{$tID} = join(";",@spl);
		}
		#print rewriting with extended tax
		open O,">$krakTaxFile" or die $!;
		foreach my $k (keys %gene2tax){
			print O "$k\t$allTaxs{$gene2tax{$k}}\n";
		}
		close O;
	}
	#die "$strings\n";
	if (0){ #kraken1
		for (my $j=0;$j< @thrs;$j++){
			open I,"<$outD/krak_$thrs[$j].out" or die "could not find file $outD/krak_$thrs[$j].out";
			open O,">$outD/krak_$thrs[$j].txt" or die "could not open file $outD/krak_$thrs[$j].txt";
			while(<I>){
				chomp;
				s/\|/;/g;
				print O $_."\n";
			}
			close I;close O;
		}	
	}
	$cmd = "$rareBin sumMat -i $GCd/$countMatrixF.gz -o $outD/krak.mat -t $NC -refD $krakTaxFile $rtkFunDelims \n";
	systemW $cmd;
}


sub printL{
	my ($msg) = @_;
	print $msg;
	print LOG $msg;
}

sub announceGeneCat{
	printL "===================================\n";
	printL "GeneCat v$version\n";
	printL "===================================\n";
	printL "Clustering at $cdhID nt id, using $numCor cores\n";
	printL "outdir: $GCdir\n";
	printL "map:    $mapF\n";
	printL "tmpdir: $tmpDir\n";
	if ($submitLocal){printL "Running in local submission mode (keep alive long time on single core)\n"} 
	else {printL "Creating fire and forget bash submission script (a bit slower since single large core machine is used)\n";}
	if ($justCDhit==0) {printL "Collating Genes from mapping file\n" ;
	} else { printL "Continuing started run in $GCdir\n";}
	printL "===================================\n";
	
}

sub specITax{
	my ($GCd,$nc) = @_;
	my $siScr = getProgPaths("specIGC_scr");
	my $cmd = "$siScr -GCd $GCd -cores $nc -tmp $tmpDir -MGset $useGTDBmg \n";
	#die "$cmd";
	print "Calculating tax abundance via SpecI's\n$cmd\n";
	systemW $cmd;
}
sub kaijuTax{#different tax assignment for gene catalog
	my ($GCd,$tmpD,$NC) = @_;
	my $kaijD = getProgPaths("kaijuDir");
	my $kaijBin = "$kaijD/./kaiju";
	my $KaDir = getProgPaths("Kaiju_path_DB");
	my $outD = $GCd."/Anno/Tax/";
	system "mkdir -p $outD" unless (-d $outD);
	system "mkdir -p $tmpD" unless (-d $tmpD);
	my @thrs = (0.01,0.02,0.04,0.06,0.1,0.2,0.3);
	my $geneFNA = "$GCd/$primaryClusterFNA";
	#die $curDB."\n";
	#paired read tax assign
	my $kaDB = "-t $KaDir/nodes.dmp -f $KaDir/kaiju_db.fmi";
	my $cmd .= "$kaijBin $kaDB -z $NC -i  $geneFNA -o $tmpD/rawKaiju.out\n";
	$kaDB = "-t $KaDir/nodes.dmp -n $KaDir/names.dmp";
	$cmd .= "$kaijD/./addTaxonNames $kaDB -i $tmpD/rawKaiju.out -o $tmpD/Kaiju1.anno -u -p \n";
	$cmd .= "sort $tmpD/Kaiju1.anno > $outD/Kaiju.anno\n";
	print "Starting kaiju assignments of the gene catalog\n";
	systemW $cmd;
	print "All kaiju assignments are done\n";
}
sub getCanopyDir{
	my ($GCd) = @_;
	my $oD = "$GCd/Canopy";
	my $xtraLstr = "";
	if ($canopyAutoCorr > 0){
		$xtraLstr = "_AC";
	}
	$oD .= "$xtraLstr/";
	return $oD;
}
sub canopyCluster{
	my ($GCd,$tmpD,$NC) = @_;
	my $canopyStone = "";
	if (@_ > 3){$canopyStone = $_[3];}
	my $canBin = getProgPaths("canopy");#"/g/bork3/home/hildebra/bin/canclus/cc_x64.bin";
	my $canGuide = getProgPaths("makeCanoGuides_scr");#"/g/bork3/home/hildebra/bin/canclus/cc_x64.bin";
	my $matF_pre = "$GCd/$countMatrixF.gz";
	my $matF = "$GCd/$countMatrixP.mat.scaled";
	my $cmd = "";
	my $doDeepCanopy = 0;
	unless (-e "$matF.gz"){
		print "Normalizing GC matrix.\n";
		$cmd .= "$rareBin normalize -i $matF_pre -t $NC -o $matF -gz\n";
		#print "Done\n";
	}
	my $oD = getCanopyDir($GCd);

	system "mkdir -p $oD" unless (-d $oD); #rm -r $oD
	my $jdeps = "";
	#die "Can't find matrix infile $matF\n" unless (-e $matF);
	#./cc.bin -i Matrix.mat.scaled.txt -o Matrix.mat.scaled.clusters -c Matrix.mat.scaled.profiles -p MGS:drama -n 32  --profile_measure 75Q --stop_criteria 0 --filter_max_top3_sample_contribution 0.7 --max_canopy_dist 0.1 --max_merge_dist 0.1
	if (!-e "$oD/clusters.txt" && !-e "$oD/profiles.txt "){
		print "Initial Canopy clustering\n";
		my $xtra = "";
		if ($canopyAutoCorr > 0){
			$xtra .= " --sampleDistMatFile $oD/smpl_dist.mat --sampleDistLog $oD/autocorr.log --sampleMinDist $canopyAutoCorr ";
		}
		$cmd .= "$canBin -i $matF.gz -o $oD/clusters.txt -c $oD/profiles.txt -p MGS $xtra --dont_use_mmap -n $NC --progress_stat_file $oD/progress.txt --profile_measure 75Q -b --stop_criteria 100000 --filter_max_top3_sample_contribution 0.7 --max_canopy_dist 0.1 --max_merge_dist 0.1\n\n";
	}
	$cmd .= "touch $oD/guideMGS.txt $oD/clusters.txt\n";
	$cmd .= "$canGuide $oD/clusters.txt $oD/profiles.txt $oD/guideMGS.txt 100\n" unless (-e "$oD/guideMGS.txt");
	#deep canopy prep
	if ($doDeepCanopy && !-e "$oD/deepClus_spea.txt"){
		$cmd .= "echo \"Deep Canopy clustering Phase II\"\n";
		$cmd .= "$canBin -i $matF.gz -o $oD/deepClus_spea.txt -r $oD/deepClus_spea_partial.txt --dont_use_mmap -n $NC -g $oD/guideMGS.txt -b --stop_criteria 0 --filter_max_top3_sample_contribution 1 --max_canopy_dist 0.5 --max_canopy_dist_part 0.2 --cag_filter_min_sample_obs 2 --dont_create_progress_stat_file --use_spearman \n";
	}
	if ($doDeepCanopy && !-e "$oD/deepClus_pear.txt"){
		$cmd .= "echo \"Deep Canopy clustering Phase III\"\n";
		$cmd .= "$canBin -i $matF.gz -o $oD/deepClus_pear.txt -r $oD/deepClus_pear_partial.txt --dont_use_mmap -n $NC -g $oD/guideMGS.txt -b --stop_criteria 0 --filter_max_top3_sample_contribution 1 --cag_filter_min_sample_obs 2 --max_canopy_dist 0.5 --max_canopy_dist_part 0.2 --dont_create_progress_stat_file \n";
	}
	#my own old options..
	#--die_on_kill --stop_criteria 250000 --cag_filter_min_sample_obs 5 --cag_filter_max_top3_sample_contribution 1 --filter_max_top3_sample_contribution 1\n";
	#die $cmd ;
	$cmd .= "\ntouch $canopyStone\n\n";
	my $canMem = int($totMem/5);
	$totMem = 250 if ($totMem<250);
	my $canDep = "";
	if (1){
		$QSBoptHR->{useLongQueue} = 0;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = "0"; #set option how much tmp space is required, and reset afterwards

		my ($jdep,$txtBSUB) = qsubSystem($qsubDir."Canopy2CL.sh",$cmd,$NC,int($canMem/$NC)."G","CAN",$jdeps,"",1,[],$QSBoptHR);
		print "Canopy MGS call send off\n"; 
		$QSBoptHR->{useLongQueue} = 0;
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		$canDep = $jdep;
	} else {
		#die $cmd."\n";
		systemW $cmd;
		print "Finished Canopy clustering.\n";
	}
	return $canDep;
	#exit(0);
}

sub nt2aa(){#takes gene cluster fna and collects respective AA from file
#probably better done with my C program
}

sub cleanUpGC(){#not used any longer
	my ($bdir,$fdir,$cdhID) = @_;
}



#simply rewrites original fasta names to counts used in my genecats
sub rewriteFastaHdIdx{ # #replaces ">MM2__C122;_23" with number from gene catalog
	my ($inf,$gene2num) = @_;
	my $sep="";
	$sep = $_[2] if (@_ >2);
	#my %gene2num = %{$hr};
	#print $gene2num{"A6M74__C100239_L=3213=_3"}."\n";
	my $numHd =0;my $cnts=0;
	my $check = 1;
	print "Rewriting $inf into $inf.tmp to add numeric headers\n";
	open I,"<$inf" or die "Can't open fasta file $inf\n"; 
	open O,">$inf.tmp" or die "can't open tmp fasta out $inf.tmp\n";
	while (my $l = <I>){
		if ($check && $numHd > 100){ 
			if ($cnts == $numHd){
				print "Seems like $inf heads were already reformated to number sheme!\n";
				 close I; close O; system "rm $inf.tmp"; return;
			} else {
				print "Rewriting fna ($inf) headers to numeric scheme\n";
				$check =0 ;
			}
		}
		
		if ($l =~ m/^>/){
			$cnts++;
			#/\x00/
		#print "$numHd  $cnts \n";
			chomp $l;
			my $name = substr($l,1); $name =~ s/\s+$//; #remove trailing spaces..
			my @spl = split /$sep/,$name;
			if ($name =~ m/^\d+$/ && (@spl < 2 || !exists($gene2num->{$spl[0]}{$spl[1]}))){ #$name})){
				$numHd++;
				print O ">".$name."\n";
			} elsif (@spl < 2) {
				die "Separator \'$sep\' not found in gene ID \'$name\', rewritign nt fna names $inf\nAborting..\n";
			} else {
				unless(exists($gene2num->{$spl[0]}{$spl[1]})){die "can not identify \'$name\' gene in index file while rewritign nt fna names $inf\n";}
				#print O ">".$gene2num->{$name}."\n";
				foreach my $pr (@{$gene2num->{$spl[0]}{$spl[1]}}){
					print O ">".$pr."\n";#;.$seq."\n";
				}
			}
		} else {
			print O $l;
		}
	}
	close I; close O;
	systemW "rm -f $inf; mv $inf.tmp $inf";
}
sub protExtract{
	my ($inD,$protXtrF) =@_;
	#gets the AA seqs for each "master" protein
	#1 read 
	#die keys %map;
	my $start = time;
	my $protF = $inD."compl.incompl.$cdhID.prot.faa";
	#my $incl = $inD."$countMatrixP.genes2rows.txt";
	system("rm -f $protF");
	print "Reading gene index\n";
	my ($geneIdxH,$numGenes) = readGeneIdxSpl($inD."$countMatrixP.genes2rows.txt","__");

	rewriteFastaHdIdx($inD.$primaryClusterFNA,$geneIdxH,$smplSep); #"compl.incompl.$cdhID.fna"
	#my %geneIdx = %{$geneIdxH};
	
	print "Writing to new proteins file: \n$protF\n";
	
	#temp DEBUG
	#my @ordG = ('Va48.6M6__C64835_L=218572=_153');
	#if ($numProts != $numGenes){ print "NUmber of genes read ($numProts) not equal to actual number of genes ($numGenes). Not enough mem?\n";}
	my $curSmpl=""; my $ctchStr = ""; my $cnt=0; 
	my @ctchAr=();
	my $ctchStrXtr = "";
	my $fSmpl;
	print "Starting Protein Extraction from source assembly folders\n";
	#collects all genes from a given sample, that is serving as seed gene for clustering
	foreach my $curSmpl (keys %{$geneIdxH}){
		next if ($curSmpl eq "xtraSmpls");
		my $curSmpl2 = $curSmpl;
		#fix MX tag at end of sample
		unless (exists ($map{$curSmpl})){$curSmpl2 = $map{altNms}{$curSmpl} if (exists ($map{altNms}{$curSmpl}));}
		#print "SMPL  $curSmpl\n";
		my $metaGD = getAssemblPath($map{$curSmpl2}{wrdir});
		my $protIn = $metaGD."/".$path2aa;
		die "prot file $protIn doesnt exits\n" unless (-e $protIn || -e $protIn.".gz");
		attachProteins3($curSmpl,$protF,$protIn,$geneIdxH->{$curSmpl},"__");
		$cnt += scalar(keys(%{$geneIdxH->{$curSmpl}}));
		print STDERR $curSmpl2." N=$cnt T=". (time - $start) ."s\n";
		$geneIdxH->{$curSmpl} = undef;
	}
	
	#---------- old routine , don't use any longer
	if (0){
		my ($geneIdxH,$numGenes) = readGeneIdx($inD."$countMatrixP.genes2rows.txt");
		my %linV = %{$geneIdxH}; #represent gene name to matrix ID
		my @ordG = sort keys %linV;
		my $numProts=scalar @ordG;
		for my $k (@ordG){
			last;
			my @tmp = split /__/,$k;
			if (@tmp>1){$fSmpl= $tmp[0]; } else {$ctchStrXtr .= "'$k' ";	next; }#this is an extra protein
			#print "$fSmpl ";
			if ($curSmpl ne $fSmpl){ #this part writes all protein IDs collected for current sample to tmp file
				#my @spl = split(/__/, $k);
				#use faidx to extract all collected gene IDs for current sample
				if ($curSmpl eq ""){
					$curSmpl = $fSmpl; 
				} else {
					unless (exists ($map{$curSmpl})){$curSmpl = $map{altNms}{$curSmpl} if (exists ($map{altNms}{$curSmpl}));}
					unless (exists ($map{$curSmpl})){#also extra protein, but with __ marker in them
						print "sk_ $curSmpl\n"; $ctchStrXtr .= $ctchStr; $ctchStr="";$curSmpl="";next;
					}
					my $metaGD = getAssemblPath($map{$curSmpl}{wrdir});
					my $protIn = $metaGD."/".$path2aa;
					die "prot file $protIn doesnt exits\n" unless (-e $protIn || -e $protIn.".gz");
					attachProteins2(\@ctchAr,$protF,$protIn,$geneIdxH);
					print STDERR $curSmpl." N=$cnt T=". (time - $start) ."s\n";
					$curSmpl = $fSmpl;$ctchStr="";
					@ctchAr=();
				}
			}
			#$ctchStr.="'$k' ";
			push(@ctchAr,$k);
			$cnt++;
			#die "$ctchStr\n" if ($cnt == 10);
		}
		#last round
		unless (exists ($map{$curSmpl})){$curSmpl = $map{altNms}{$curSmpl} if (exists ($map{altNms}{$curSmpl}));}
		unless (
			exists ($map{$curSmpl})){$ctchStrXtr .= $ctchStr; #prob extra protein
		}else{	
			print $curSmpl." $cnt\n";
			my $metaGD = getAssemblPath($map{$curSmpl}{wrdir});
			my $protIn = $metaGD."/".$path2aa;
			attachProteins2(\@ctchAr,$protF,$protIn,$geneIdxH);
		}
	}
	
	print "rewritten $cnt proteins, expected $numGenes\n";
	
	#extra added proteins (not from MATAFILER assembly)
	unlink "$inD/tmp.txt";
	if ($protXtrF ne ""){
		open O,">$inD/tmp.txt";print O $ctchStrXtr;close O;
		#die ("extra\n$inD/tmp.txt\n");
		attachProteins("$inD/tmp.txt",$protF,$protXtrF,$geneIdxH);
	}

	#new cluster numbers and one new file, my format, with Idx
	combineClstr("$inD/compl.incompl.$cdhID.fna","$inD/$countMatrixP.genes2rows.txt") ;

}
sub combineClstr(){
	my ($clstr1,$idx) = @_;
	my $clusMode=2;my $clstr = $clstr1.".ctsv";
	if (-e "$clstr1.clstr"){
		$clusMode = 1;$clstr = $clstr1.".clstr";
	}
	
	#currently can only be run first time!
	print "Combining cluster strings.. $clstr\n$idx\n";
	#tmp out files, copied over to correct locations later!
	open I,"<$clstr"; open O,">$clstr.2"; open Oi,">$clstr.idx2"; open C,"<$idx"; 
	my $chLine = <C>; my $newOil=0;
	#counts if already formated?
	my $evidence=0; my $eviNo=0;
	my $oil = ""; #collects for current cluster genes
	print Oi "#Gene	members\n";
	while (my $line = <I>){
		chomp $line;
		if ($line =~ m/^>/){
			$line =~ m/^>(.*)/;
			#chop $oil;
			print Oi $oil."\n" if ($oil ne "");
			$chLine = <C>;
			my @spl = split(/\t/,$chLine);
			if ($1 eq $spl[0]){
				$evidence++;
				if ($evidence>10 && $eviNo == 0){
					print "It seems like index file was already created, aborting coversion..\n";
					systemW "rm $clstr.idx2 $clstr.2";
					return;
				}
			} else {
				$eviNo ++;
			}
			if ($spl[1] ne $line){
				die "Can;t match \n$line \n$spl[1]\n";
			}
			
			$line = ">$spl[0]";
			$oil = $spl[0]."\t";
			$newOil = 1;
			#die "FND :: $line\n";
		} else {
			$line =~ m/\s+(>.*)\.\.\.\s.*/;
			if ($newOil){
				$oil.=$1;
				$newOil=0;
			} else {
				$oil.=",".$1;
			}
		}
		print O $line."\n";
	}
	#insert last entry
	print Oi $oil."\n";
	systemW "rm $clstr; mv $clstr.2 $clstr";
	systemW "rm $clstr.idx" if (-e "$clstr.idx");
	systemW "mv $clstr.idx2 $clstr.idx";
	print "Done rewriting cluster numbers & creating cluster index\n";
	close I; close O;close C; close Oi;
}



sub clusterFNA($ $ $ $ $ $ $ $ $ $){
	my ($inFNA, $oFNA, $aS,$aL, $ID, $numCor, $gfac, $tmpD, $useMMseqs, $totMemCl) = @_;
	my $cmd = "";
	my $mmS2clstr = getProgPaths("mmS2clstr");
	if ($ID >1){$ID=$ID/100;}
	
	if ($useMMseqs){#mmseq2 clustering  #$clustMMseq
		my $tmpD2 = "$tmpD/mmS/";
		$cmd .= "rm -rf $tmpD2\nmkdir -p $tmpD2\n";
		my $covMin = $aL; $covMin = $aS if ($aS < $aL);
		if ($gfac){
			$cmd .= "$mmseqs2Bin easy-cluster $inFNA $oFNA $tmpD2 -c $covMin --cov-mode 0  --min-seq-id $ID --alignment-mode 3 --threads $numCor --dbtype 2 --spaced-kmer-mode 1 --mask 0 --split-memory-limit ". int(${totMemCl}*0.85)."G --min-aln-len 100 --sort-results 1 --cluster-reassign 0 -v 3 \n";
		} else {
			$cmd .= "$mmseqs2Bin easy-linclust $inFNA $oFNA $tmpD2 -c $covMin --cov-mode 0  --min-seq-id $ID --alignment-mode 3 --threads $numCor --dbtype 2 --spaced-kmer-mode 1 --mask 0 --split-memory-limit ". int(${totMemCl}*0.85)."G --min-aln-len 100 --sort-results 1 -v 3 \n"; #--cluster-reassign 1 
		}
		$cmd .= "$mmS2clstr ${oFNA}_cluster.tsv $oFNA.clstr\n";
		#die $cmd;
		$cmd .= " rm -f $oFNA ${oFNA}_all_seqs.fasta ${oFNA}_cluster.tsv;\n mv ${oFNA}_rep_seq.fasta $oFNA;\n\n";
	} elsif(1) {
		#	$defaultsCDH = "-d 0 -c 0.$cdhID -g 0 -T $numCor -M ".int(($totMem+30)*1024) if (@ARGV>3);
		$cmd .= $cdhitBin."-est -i $inFNA -o $oFNA -n 9 -mask NX -G 1 -r 0 -aS $aS -aL $aL -d 0 -c $ID -g $gfac -T $numCor -M ".int(($totMem+30)*1024)."\n";
	} else {
			die "no longer supported vsearch clustering\n";
		$cmd .= "gunzip $bdir/compl.srt.fna.gz\n" if (-e "$bdir/compl.srt.fna.gz" && !-e "$bdir/compl.srt.fna");
		$cmd .= $vsearchBin." --cluster_fast $bdir/compl.srt.fna --consout $tmpDir/compl.$cdhID.fna --id 0.$cdhID --strand plus --threads $numCor --uc $tmpDir/compl.$cdhID.uc";
	}
	return $cmd;
}


sub writeBucket(){
	my ($OCOMPLar,$O3Par,$O5Par,$OINCar,$bdir,$bnum) = @_;
	if ($toLclustering){return;}
	systemW("mkdir -p $bdir");
	#my @OCOMPL = @{$OCOMPLar}; 
	if (@{$OCOMPLar} ==0){die "no genes found!";}
	if ($toLclustering){#no length sorting
		#open O,">$bdir"."compl.fna"; foreach( @OCOMPL ){print O $_;} close O;
		#open O,">$bdir"."5Pcompl.fna"; foreach( @O3P ){print O $_;} close O;
		#open O,">$bdir"."3Pcompl.fna"; foreach( @O5P ){print O $_;} close O;
		#open O,">$bdir"."incompl.fna";foreach( @OINC ){print O $_;} close O;
	} else {
		open O,">$bdir"."compl.fna" or die "Can't open B0 compl.fna\n"; 
		foreach(sort {length $b <=> length $a} @{$OCOMPLar} ){print O $_;} close O; @{$OCOMPLar} = ();
		my @O5P = @{$O5Par};
		open O,">$bdir"."5Pcompl.fna" or die "Can't open B0 5P.fna\n"; foreach(sort {length $b <=> length $a} @O5P ){print O $_;} close O; @O5P=();
		my @O3P = @{$O3Par};  
		open O,">$bdir"."3Pcompl.fna" or die "Can't open B0 3P.fna\n"; foreach(sort {length $b <=> length $a} @O3P ){print O $_;} close O; @O3P=();
		 my @OINC = @{$OINCar}; 
		open O,">$bdir"."incompl.fna" or die "Can't open B0 incompl.fna\n";foreach(sort {length $b <=> length $a} @OINC ){print O $_;} close O;  @OINC=();
	}
}
sub readFasta($){
  my ($fil) = @_;
  my %Hseq;
  if (-z $fil){ return \%Hseq;}
  open(FAS,"<","$fil") || die("Couldn't open FASTA file $fil.");
    
     my $temp; 
     my $line; my $hea=<FAS>; chomp ($hea);
      my $trHe = ($hea);
      #my @tmp = split(" ",$trHe);
      #$trHe = substr($tmp[0],1);
      # get sequence
    while($line = <FAS>)
    {
      #next if($line =~ m/^;/);
      if ($line =~ m/^>/){
        chomp($line);
        $Hseq{$trHe} = $temp;
        $trHe = ($line);
       # @tmp = split(" ",$trHe);
		#$trHe = substr($tmp[0],1);
		$trHe =~ s/\|//g;
        $temp = "";
        next;
      }
    chomp($line);
    $line =~ s/\s//g;
    $temp .= ($line);
    }
    $Hseq{$trHe} = $temp;
  close (FAS);
    return \%Hseq;
}
sub rewriteClusNumbers($ $ ){
	my ($infna,$newCnt) = @_;
	print "Rewriting Cluster numbers: $newCnt in $infna.clstr\n";
	my $incls = $infna.".clstr";
	my $ocls = $infna.".clstr.new";
	my $mem = 0;
	open I,"<$incls" or die "RewriteClusNum: Can't open i $incls\n";
	open O,">$ocls" or die "RewriteClusNum: Can't open o $ocls\n";
	while(my $line = <I>){
		if ($line =~ m/^>Cluster/){
			#>Cluster 0
			#$line =~ s/>Cluster \d/>Cluster 0/;
			print O ">Cluster $newCnt\n";
			$newCnt++;
		} else {
			$mem++;
			print O $line;
		}
	}
	close I; close O;
	systemW("rm -f $incls\nmv $ocls $incls \n");
	$newCnt--; #to keep accurate track of clus numbers..
	return($newCnt,$mem);
}
sub writeIDs($ $){
	my ($hr, $of) = @_;
	my %ids = %{$hr};
	open O,">$of" or die "can't open $of cluster id file\n";
	foreach my $k (keys %ids){
		print O $k."\t".$ids{$k}."\n";
	}
	close O;
}

sub addCOGgenes{
	my ($inD,$GCd,$inCclN,$outFfna,$outFcls) = @_;
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#look for COG genes
	print "$inD/$COGdir/\n";
	opendir(DIR, "$inD/$COGdir/") or die $!;
	my @cogfiles = sort (grep {/.*\.fna\.clstr$/ } readdir(DIR)); close DIR;#&& -f "$inD/COG/$_" 
	my $logstr = ""; my $inCclNmember = 0;
	#print "@cogfiles"."\n";
	if (@cogfiles < 40){die"\n\nless than 40 FMG genes!\n".@cogfiles."\n@cogfiles\n";}
	my $COGcnt =0; my $inCclN2 = 0; 
	open Oc,">$GCdir/${COGdir}.subset.cats1";  #DEBUG.. switch to .cats later..
	foreach my $cogFNAcl (@cogfiles){
		$inCclN+=100; $COGcnt++;
		my $cogFNA = $cogFNAcl; $cogFNA =~ s/\.clstr//;
		($inCclN2,$inCclNmember) = rewriteClusNumbers("$inD/$COGdir/".$cogFNA,$inCclN);
		#die "$cogFNA\ncogFNA\n";
		systemW("cat $inD/$COGdir/$cogFNAcl >> $outFcls");
		systemW("cat $inD/$COGdir/$cogFNA  >> $outFfna");
		$cogFNA =~ /([^\.]+)\./; my $cog = $1;
		$logstr .= "$cog	$inCclN	$inCclN2\n";
		my @tmpA = ($inCclN .. $inCclN2);
		print Oc $cog . "\t". scalar(@tmpA) . "\t" . join(",",@tmpA)."\n"; 
		$inCclN = $inCclN2;
	}
	close Oc;
	print "\nFound $COGcnt $COGdir genes\n";
	#print "$logstr\n";
	print "added IDs to $outFcls\n\n";
	open O,">$GCdir/LOGandSUB/${COGdir}.clusN.log" or die $!; print O $logstr; close O;
}

sub mergeClsSam(){
	my ($inD,$idP,$GCd) = @_;
	
	my $outFfna = $inD.$primaryClusterFNA;
	my $outFcls = $inD.$primaryClusterCLS;
	
	#single step clustering? only needs to have the FMG genes added..
	if ($allinClust){
		unless (-d "$inD/$COGdir/" || -d "$GCd/B0/$COGdir/"){print"No COG specific genes found\n$inD/$COGdir/\n$GCd/B0/$COGdir/\n";return;}
		my $inCclN = 0;
		print "$outFcls\n\n";
		die "Can't find $outFcls\n" unless (-e $outFcls);
		my $tlW = `tail -n 800 $outFcls| grep '^>' | tail -n 1`;chomp $tlW;
		print $tlW."\n";
		if ($tlW !~ m/>Cluster (\d+)/){
			$tlW = `tail -n 25000 $outFcls| grep '^>' | tail -n 1`;chomp $tlW;
		}
		$tlW =~ m/>Cluster (\d+)/;
		$inCclN = $1;
		
		#main step...
		addCOGgenes($inD,$GCd,$inCclN,$outFfna,$outFcls);
		return;
	}
	
	my @samFs = ($inD."incompl.$idP.align.sam",$inD."35compl.$idP.align.sam");
	my $completeFNA = $inD."compl.$idP.fna";
	my $clFile = $completeFNA.".clstr";
	die "Can't find required input $completeFNA" unless (-e $completeFNA);
	die "Can't find required input $clFile" unless (-e $clFile);

	#system "cp $clFile $clFile.before"; #TODO, remove
	my $logf = "$inD/log/Cluster.log";
	systemW("mkdir -p $inD/log/");
	open LOG,">$logf";
	#read current clusters.. this needs to be extended by sam hits
	my ($hr1,$hr2,$totN,$totM,$idsHr) = readCDHITCls($clFile);
	print LOG "ComplGeneClus	$totN\nComplGeneClusMember	$totM\n";
	my %lnk  = %{$hr2};
	my %cls = %{$hr1};
	my %idsDistr = %{$idsHr};
	writeIDs($idsHr,"$inD/cluster.ids.primary");
	my %Hitin; my $remSeqs;
	my $totCnt = 0;
	my @preTags = ("incompl","P35compl");
	#read (both) sam with hits to complete clusterss.
	my @unclusters;
	my $finalUnclusterd = "$inD/incompl.NAl.fna";
	for (my $k=0;$k<2;$k++){
		my $preTag = $preTags[$k];
		die "Can't find samfile $samFs[$k]\n" unless (-e $samFs[$k]);
		my $unclusFna = "$inD/$preTag.NAl.$idP.fna";
		push(@unclusters,$unclusFna);
		systemW "rm -f $unclusFna" if (-e "$unclusFna");
		#was another format previously..
		systemW "cp $inD/$preTag.NAl.pre.$idP.fna $unclusFna" if (-e "$inD/$preTag.NAl.pre.$idP.fna");
		#die();
		open my $OO,">>","$unclusFna" or die "Can't open $unclusFna\n"; 
		print $OO "\n";
		$hr1 = readSam($samFs[$k],$OO); #,$remSeqs)
	#	print OO $remSeqs;
		close $OO;
		
		#sam read and identified genes without good hit..
		%Hitin = %{$hr1};
		#first, add hits to complete gene clusters
		foreach my $hit (keys %Hitin){
			unless (exists($lnk{$hit})){die "Can't find link $hit\n";}
			my @splClusters = split(/\n/,$cls{$lnk{$hit}}) ;
			my $clCnt = scalar(@splClusters);
			foreach my $spl (split(/\n/,$Hitin{$hit})){
				$cls{$lnk{$hit}} .= "\n$clCnt\t".$spl;
				$totCnt++;$clCnt++;
			}
			#die "" if ($totCnt > 10);
		}
	}

	#write the new clstr file out
	open O,">$outFcls"; my $totCls=0;
	foreach my $cl (keys %cls){
		$totCls++;
		my $ostr = $cl."\n".$cls{$cl};
		print O $ostr."\n";
	}
	close O;
	#print LOG "final Clusters (incomplete + complete) : $totCls\n";
	print "Results $outFcls\n";
	
	
	#now unclustered genes need to be clustered
	systemW("cat ".join(" ",@unclusters)." > $finalUnclusterd");
	my $remClus = "$inD/incompl.rem.$cdhID.fna";
	my $cmd = clusterFNA( "$finalUnclusterd", "$remClus",0.7,0.3,"$cdhID",$numCor,0,$inD,$clustMMseq,$totMem);
	systemW($cmd) unless (-e $remClus);
	
	#old routine, not needed (too convoluted programming flow)
	#my $Cls2ndFNA = secondaryCls($inD,$idP);
	
	
	my ($inCclN,$inCclNmember) = rewriteClusNumbers($remClus,$totCls);
	print "Concatenating Cluster files\n";
	systemW("cat $remClus.clstr >> $outFcls");
	print "Concatenating Cluster Seed fna files\n";
	systemW("cat $completeFNA $remClus > $outFfna");
	print LOG "IncomplComplGeneClus	$inCclN\n";
	print LOG "IncomplComplGeneClusMember	".($inCclNmember+$totM)."\n";
	close LOG;
	
	unless (-d "$inD/$COGdir/" || -d "$GCd/B0/$COGdir/"){print"No COG specific genes found\n";return;}
	
	addCOGgenes($inD,$GCd,$inCclN,$outFfna,$outFcls);
}

sub writeMG_COGs{
	my ($GCd) = @_;
	die "cant' work, since ordering not given. used annotateMGwMotus.pl instead";
	my $FMGd = "$GCd/FMG/";
	system "mkdir -p $FMGd";
	open I,"<$GCd/LOGandSUB/${COGdir}.clusN.log";
	while (my $l = <I>){
		chomp $l;
		my @spl = split /\t/,$l;
		my @range = $spl[1] .. $spl[2];
		my $cmd = "$samBin faidx $GCd/compl.incompl.95.fna ". join (" ", @range) . " > $FMGd/$spl[0].gc.fna\n";
		system $cmd;
		$cmd = "$samBin faidx $GCd/compl.incompl.95.prot.faa ". join (" ", @range) . " > $FMGd/$spl[0].gc.faa";
		system $cmd;
	}
	close I;
}

sub secondaryCls(){
	my ($inD,$cdhID) = @_;
	my $cmd="";
	die "no longer used\nsecondaryCls\n";
	my $outfna = "$inD/incompl.rem.$cdhID.fna";
	if (-e $outfna){return($outfna);}
	#$cmd .= $cdhitBin."-est -i $bdir/P35compl.NAl.$cdhID.fna -o $bdir/P35compl.$cdhID.fna -n 9 -G 0 -M 5000 -aL 0.5 -aS 0.95 $defaultsCDH\n";
	#die("cdf\n");
	#system("cat $inD/P35compl.NAl.$cdhID.fna >> $inD/incompl.NAl.$cdhID.fna\n");
	#size sort
	#my $suc = systemW("cat $inD/P35compl.NAl.$cdhID.fna $inD/incompl.NAl.$cdhID.fna | perl -e 'while (<>) {\$h=\$_; \$s=<>; \$seqs{\$h}=\$s;} foreach \$header (reverse sort {length(\$seqs{\$a}) <=> length(\$seqs{\$b})} keys \%seqs) {print \$header.\$seqs{\$header}}' > $inD/incompl.NAl.srt.$cdhID.fna");
	my $suc = systemW("cat $inD/P35compl.NAl.$cdhID.fna $inD/incompl.NAl.$cdhID.fna  > $inD/incompl.NAl.fna");
	#my $hr = readFasta("$inD/P35compl.NAl.$cdhID.fna");
	#my $hr2 = readFasta("$inD/incompl.NAl.$cdhID.fna");
	#my %remGenes = ( %{$hr}, %{$hr2} ); $hr=0;$hr2=0;
	#my @keys = sort { length($remGenes{$a}) <=> length($remGenes{$b}) } keys(%h);

#	systemW($cdhitBin."-est -i $inD/incompl.NAl.srt.$cdhID.fna -o $outfna -n 9 -G 0 -aL 0.3 -aS 0.8 $defaultsCDH\n") unless (-e $outfna);
	$cmd ="";
	#$cmd .= sortFNA($inD,"incompl.NAl",1,$tmpDir,$numCor);
	#$cmd .= clusterFNA( "$inD/incompl.NAl.fna", "$outfna",0.7,0.3,"$cdhID",$numCor,1,$tmpD);
	systemW($cmd);
	#$cmd .= "rm -f $inD/incompl.NAl.$cdhID.fna\n";
	#$cmd .= "rm -f $inD/P35compl.NAl.$cdhID.fna $inD/35compl.$cdhID.align.sam $inD/incompl.$cdhID.align.sam\n";
	return($outfna);
}

sub readCDHITCls(){
	my ($iF) = @_;
	my %retCls; my %retRepSeq; my %clsIDs;
	open I,"<$iF";
	my $clName = "";
	my $clNum=0; my $totalStore=0;
	
	while (my $line = <I>){
		chomp $line;
		if ($line =~ m/^>/){#open new cluster
			$clName = $line; $clNum++; $totalStore++; next;
		}
		$totalStore++;
		if (exists($retCls{$clName})){
			$retCls{$clName} .= "\n".$line;
		} else {
			$retCls{$clName} = $line;
		}
		if ($line =~ m/\*$/){#cluster seed
			#my @tmp = split(/\s*/,$line);
			$line =~ m/>(.*)\.\.\./;
			#print $1."\n";;
			$retRepSeq{$1} = $clName;
		} else {
			$line =~ m/>(.*)\.\.\. at .\/([0-9\.]+)%/;
			if (!exists $clsIDs{$clName} ){
				$clsIDs{$clName} = $2;
			} else {
				$clsIDs{$clName} .= ",".$2;
			}
		}
	}
	return(\%retCls,\%retRepSeq,$clNum, $totalStore,\%clsIDs);
}



sub geneCatFunc_emapper{
	#tmpD is node-local tmp
	my ($GCd,$tmpD, $ncore,$doClean,$fastaSplits,$stone) = @_;
	
	my $query = "$GCd/compl.incompl.95.prot.faa";
	my $mem = 55; #default mem/job
	my $outD = $GCd."/Anno/Func/emapper/";
	my $curDB = getProgPaths("eggNOGm_path_DB");
	my $emapper = getProgPaths("emapper");
	my $qsubDir2 = "$qsubDir/Funct/";
	my $shrtDB = "emap";
	system "mkdir -p $qsubDir2" unless (-d $qsubDir2);
	system "mkdir -p $outD" unless (-d $outD);
	#use a different dir for qsub jobs to keep main qsub dir clean
	$QSBoptHR->{qsubDir} = $qsubDir2;
	my $doQsub = 1;my $calcDia = 1;my @jdeps;
	my $jdep = "";
	my $tarAnno3 = "$outD/MF.emapper.annotations";
	my $splDir = "$GLBtmp/eggNOGmapper/";
	if (!-e $tarAnno3 || !-s $tarAnno3){
		#setup cluster for diamond focused job
		my @preCons = @{$QSBoptHR->{constraint}};
		push(@{$QSBoptHR->{constraint}}, $avx2Constr);#--constraint=sse4
		my $preHDDspace=$QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = ($mem*1.4) . "G";
		print "Splitting FASTAs into $fastaSplits chunks in $splDir\n";
		my $ar = splitFastas($query,$fastaSplits,$splDir);
		print "Done splitting, submitting $fastaSplits jobs\n";
		my @subFls = @{$ar};
		my $i=0;
		foreach my $f (@subFls){
			#system "mkdir -p $tmpD/$i/" unless (-d "$tmpD/$i/");
			#my $cmd = "emapper.py -m diamond --override --temp_dir $tmpD/$i/ --data_dir $curDB --no_annot --no_file_comments --cpu $ncore -i $f -o $f;"; 			my $outF = "$f.emapper.seed_orthologs";

			my $cmd = "";
			$cmd .= "mkdir -p $tmpD/$i/\n";
			$cmd .= "emapper.py -m diamond --dbmem --override --itype proteins --temp_dir $tmpD/$i/ --data_dir $curDB --no_file_comments --cpu $ncore -i $f -o $f;\n"; 
			my $outF = "$f.emapper";

			if ($calcDia && (!-e $outF || !-s $outF) ){
				#print "$cmd\n";
				if ($doQsub){
					my ($jobName,$mptCmd) = qsubSystem($qsubDir2."D$shrtDB.$i.sh",$cmd,$ncore,($mem/$ncore)."G","eMAP$i","","",1,[],$QSBoptHR); #$jdep.";".
					push(@jdeps,$jobName);
					#die "$jobName\n";
				} else {
					systemW $cmd;
				}
				#print $qsubDir."Diamond.sh\n";
			}
			$i++;
		}
		my $ncore2 = 1;#$ncore;
		@{$QSBoptHR->{constraint}} = @preCons;
		$QSBoptHR->{tmpSpace} = $preHDDspace;
		#my $tarAnno = "$outD/D$shrtDB.emapper.seed_orthologs";
		#my $tarAnno2 = "$outD/D$shrtDB";
		my $cmd = "";
		#$cmd .= "#concatenating separate eggNOG-diamond output\n";
		#$cmd .= "cat $GLBtmp/eggNOGmapper/*.emapper.seed_orthologs > $tarAnno\n";
		#$cmd .=  "#emapper.py --data_dir $curDB --annotate_hits_table $tarAnno --no_file_comments -o $tarAnno2 --cpu $ncore2 --dbmem\n";
		#$cmd .= "\nexit(1)\n";
		#$cmd .= "cat $GLBtmp/eggNOGmapper/*.emapper.annotations > $tarAnno2\n ";
		$cmd .= "head -q -n 1 " . $subFls[0] . ".emapper.annotations > $tarAnno3\n";
		$cmd .= "tail -q -n +2 $GLBtmp/eggNOGmapper/*.emapper.annotations >> $tarAnno3\n";
		my $eSpl = getProgPaths("eggNOGspl_scr");
		$cmd .= "#splitting eggNOG annotations in multiple categories that can be summed up to matrices\n$eSpl $tarAnno3\n";
		#run 
		#$cmd = "" if (-e "$tarAnno3");
		my ($jobName,$mptCmd) = qsubSystem($qsubDir2."CombineEMAP.sh",$cmd,$ncore2,(80/$ncore2)."G","${shrtDB}_comb",join(";",@jdeps),"",1,[],$QSBoptHR); 
		#$jdep.";".
		
		$jdep = $jobName;
	}
# create abundance tables now..
	my $matThr = 4;
	my @emapCats  = ("eggNOGmapper_CAZy","eggNOGmapper_EC","eggNOGmapper_GO","eggNOGmapper_NOG","eggNOGmapper_BIGG","eggNOGmapper_PFAM",
	#"eggNOGmapper_KGP","eggNOGmapper_KGM"
	"eggNOGmapper_KGM", "eggNOGmapper_KGP");
	@jdeps = ();
	foreach my $EMC (@emapCats){
		#"$EMC.geneAss.gz";
		$EMC =~ m/_(\S+)/; my $shrt = "EM.$1";
		my $cmd = "$rareBin sumMat -i $GCd/$countMatrixF.gz -o $outD/$shrt -t $matThr -refD $outD/$EMC.geneAss $rtkFunDelims -extHiera -hieraSrtDown\n";
		$cmd .= "echo \"DONE matrix creation\"\n";
		#my $jobName = "sum$shrt";
		my ($jobName,$mptCmd) = qsubSystem($qsubDir2."$shrt.sh",$cmd,$matThr,(50/$matThr)."G","${shrt}_mat",$jdep,"",1,[],$QSBoptHR); #$jdep.";".
		push(@jdeps,$jobName);
	}
	#clean up and marking stone that all worked out fine..
	my $clnCores = 4;
	my $cmd = "";
	$cmd .= "\ntouch $stone\n\n";
	$cmd .= "$pigzBin -p $clnCores $tarAnno3 $outD/*.geneAss;\n";
	$cmd .= "rm -f -r $GLBtmp/eggNOGmapper  $splDir;\n" ; #unless (-e $tarAnno && -s $tarAnno)
	my ($jobName,$mptCmd) = qsubSystem($qsubDir2."CleanEMAP.sh",$cmd,$clnCores,(70/$clnCores)."G","${shrtDB}_CLN",join(";",@jdeps),"",1,[],$QSBoptHR); #$jdep.";".




}

sub geneCatFunc{
	my ($GCd,$tmpD, $DB, $ncore,$doClean) = @_;
	my $query = "$GCd/compl.incompl.95.prot.faa";
	my $outD = $GCd."/Anno/Func/";
	die "-functAligner has to be \"diamond\" or \"foldseek\"!\n" if ($funcAligner ne "diamond" && $funcAligner ne "foldseek");
	
	#my $DB = "NOG";
	#my $ncore = 40; 
	
	
	my $curDB = $DB; #"NOG";#CZy,ABRc,KGM,NOG
	my $qsubDir2 = "$qsubDir/Funct/";
	system "mkdir -p $qsubDir2" unless (-d $qsubDir2);
	$QSBoptHR->{qsubDir} = $qsubDir2;
	
	my %optsDia = (eval=>1e-8,percID=>25,minPercSbjCov=>0.5,fastaSplits => $fastaSplits,ncore=>$ncore,align=>$funcAligner,
			splitPath=>$GLBtmp,keepSplits=>!$doClean,redo=>$doClean, minAlignLen=>30, minBitScore=>45);
			
			
	my ($allAss,$jdep) = assignFuncPerGene($query,$outD,$tmpD,$curDB,\%optsDia,$QSBoptHR,(!-e "$outD/${curDB}L0.txt")) ;
	my $tarAnno = "${allAss}geneAss.gz";
	my $tmpP2 = "$tmpD/CNT_1e-8_25//";
	#create actual COG table
	my $cmd = "";
	$tarAnno =~ s/\.gz$//;
	
	#$cmd	.= "gunzip $tarAnno.gz\n";
	if ($curDB eq "ABRc"){
		$cmd .= "zcat $tarAnno.gz | sed 's/,/\\|/g' |sed 's/\\t/;/g' | sed 's/;/\\t/' > $tarAnno\n";
	} else {
		$cmd .= "zcat $tarAnno.gz | sed 's/\\t/;/g' | sed 's/;/\\t/' > $tarAnno\n";
	}
	
	#$cmd .= "$rareBin sumMat -i $GCd/$countMatrixF -o $outD/${shrtDB}L1.mat -refD $GCd/NOGparse.NOG.GENE2NOG; gzip $GCd/NOGparse.NOG.GENE2NOG.gz\n";
	my $matThr = 4;
	$cmd .= "$rareBin sumMat -i $GCd/$countMatrixF.gz -o $outD/$curDB -t $matThr -refD $tarAnno $rtkFunDelims \n";
	$cmd .= "rm $tarAnno\n" if (length($tarAnno) > 2);
	#gzip $tarAnno\n";
	#copy interesting files to final dir
	#if ($curDB eq "ABRc"){
#		$cmd.= "zcat $tmpP2/ABRcparse.ALL.cnt.CATcnts.gz > $outD/ABR_res.txt\n" ;
#	}
#	if ($curDB eq "CZy"){
#		$cmd.= "zcat $tmpP2/CZyparse.ALL.cnt.CATcnts.gz > $outD/CZySubstrates.txt\n";
#		$cmd.= "zcat $tmpP2/CZyparse.CZy.ALL.cnt.cat.cnts.gz > $outD/CZyEnzymes.txt\n";
#	}
#	if ($curDB eq "TCDB"){
#		$cmd.= "zcat $tmpP2/TCDBparse.ALL.cnt.CATcnts.gz > $outD/TCDB.cats.txt\n";
#	}
#systemW $cmd."\n";
#die "$outD/${curDB}L0.txt\n";
	my $modCmd = "";
	if ($curDB eq "KGM"){
		$modCmd = calc_modules("$outD/${curDB}L0.txt","$outD/modules/",0.5,0.5,0);#$ModCompl,$EnzCompl);
	}
	#die "$modCmd\n";

	$cmd = "" if (-e "$outD/${curDB}L0.txt");
	#die "$cmd.$modCmd\n";
	if (0 && -e $tarAnno){ #already exists
		systemW $cmd ;
	} else {
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		($cmd,$jdep) = qsubSystem($qsubDir2."${curDB}_matrix.sh",$cmd.$modCmd,$matThr,"12G","${curDB}_mat",$jdep,"",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}


}

sub readSam($$){
	my ($iF,$OO) = @_;#$fref,
	#my %fas = %{$fref};
	open I,"<$iF" or die "Can't open $iF\n";
	my $cnt = 0; my $totLines = 0;
	my $add2file=0; my $add2cls=0;
	my $bwt2sam = 0; #assumme bowtie2, switch to minimap2
	my %ret; #my $fasStr = "";
	print $OO "\n";
	while (my $line = <I>){
		chomp $line; $totLines++;
		my @sam = split(/\t/,$line);
		my $qu = $sam[0]; my $ref = $sam[2];
		my $refL = length($sam[9]);
		next if ($sam[1] & 0x2048); #completely ignore..
		if ($sam[1] & 0x4 ){#not a hit, minimap2
			#print "$sam[1]\n";
			my $nqu = ">".$qu;
			if ($refL > 100){ #too short hits don't need to be attached..
				print $OO $nqu."\n".$sam[9]."\n";#$fas{$nqu}."\n";
				$add2file++;
			}
			next;
		}
		#die "$sam[0] $sam[1] $sam[2] $sam[3]\n$line\n";
		#if ($qu =~ m/MM28__C41733_L=413;_1/){print "TRHERE\n";}
		my $xtrField = join("\t",@sam[11..$#sam]);
		#die $xtrField;
		my $pid = 100; my $acc=1;
		if ($bwt2sam==1) {#bowtie2
			if ($xtrField =~ m/XM:i:(\d+)\s.*XO:i:(\d+)\s.*XG:i:(\d+)/){
				$pid = $1/($refL-$2-$3);
				$acc = 0 if (($2+$3)/$refL > 0.1 || $pid > 0.05);
			} else {$bwt2sam=0;}
		}
		if ($bwt2sam==0) {#mini2
			if ($xtrField =~ m/NM:i:(\d+).*de:f:([0-9\.]+)/){
				$pid = $2; #$mismatches = $1;
				$acc = 0 if ($pid > 0.05);
			} else {$bwt2sam=1;}
		}
		#print "$1 $2 $3 $refL ".$1/$refL." ".($2+$3)/$refL."\n";
		#95% id || 90% seq length
		if ( !$acc){ #not good enough hit criteria, attach to fasta
			my $nqu = ">".$qu;
			if ($refL > 100){ #too short hits don't need to be attached..
				print $OO $nqu."\n".$sam[9]."\n";#$fas{$nqu}."\n";
				$add2file++;
			}
			#die $nqu."\n".$sam[9]."\n";
			next;
		}
		#my $mism = $1; my $gaps = $2+$3;
		#$cnt++;
		my $clsStr = $refL."nt, >".$qu."... at +\/". int(( (1.0-$pid)*100) ) .".00%";
		if (exists($ret{$ref})){
			$ret{$ref} .= "\n".$clsStr
		} else {
			$ret{$ref} = $clsStr;
		}
		$add2cls++;
		#if ($qu =~ m/MM28__C41733_L=413;_1/){print "print\n";}

		#print $ret{$ref}."\n";
		#die if ($cnt == 10);
	}
	close I;
	print LOG $add2cls." hits to clusters, $add2file added FNAs to be reclustered (of ".$totLines." lines) in $iF\n";
	print $add2cls." hits to clusters, $add2file added FNAs to be reclustered (of ".$totLines." lines) in $iF\n";
	#die;
	return (\%ret);
}

sub FOAMassign{
	my ($GCd,$tmpD, $DB) = @_;
	my $query = "$GCd/compl.incompl.95.prot.faa";
	my $fastaSplits=10;
	my $ar = splitFastas($query,$fastaSplits,$GLBtmp."DB/");
	my @subFls = @{$ar};
	my @jdeps; my @allFiles;
	my $N = 20;my $jdep=""; my $colSel = 4;
	my $tmpSHDD = $QSBoptHR->{tmpSpace};
	$QSBoptHR->{tmpSpace} = "150G"; #set option how much tmp space is required, and reset afterwards

	for (my $i =0 ; $i< @subFls;$i++){
		my $tmpOut = "$tmpD/$DB.hmm.dom.$i";
		my $outF = "$GCd/assig.$DB.$i";
		my $cmd = "mkdir -p $tmpD\n";
		if ($DB eq "FOAM"){
			$cmd .= "$hmmBin3 --cpu $N -E 1e-05 --noali --domtblout $tmpOut $FOAMhmm $subFls[$i] > /dev/null\n";
			
		} elsif ($DB eq "ABR") {
			#--tblout=<output tab-separated file>
			$cmd .= "$hmmBin3 --cpu $N --domtblout $tmpOut --cut_ga $ABresHMM $subFls[$i] > /dev/null\n";
			$colSel = 5;#select gene name
		}
		$cmd .= "sort $tmpOut > $tmpOut.sort\n";
		#DEBUG
		#$cmd .= "cp $tmpOut.sort $GCd\n";
		
		$cmd .= "python $hmmBestHitScr $tmpOut.sort > $tmpOut.sort.BH\n";
		$cmd .= "awk '{printf (\"%s\\t%s\\n\", \$1,\$$colSel)}' $tmpOut.sort.BH |sort > $outF\n";
		$cmd .= "rm -f -r $tmpOut* $subFls[$i]\n";
		my $jobName = "$DB"."_$i";
		#die $cmd."\n";
		my ($cmdRaw,$jdep) = qsubSystem($qsubDir."$DB$i.sh",$cmd,$N,"1G",$jobName,"","",1,[],$QSBoptHR);
		push(@jdeps,$jdep);
		push(@allFiles,$outF);
		#if ($i==5){die;}
	}
	$QSBoptHR->{tmpSpace} = $tmpSHDD; 
	#last job that converges all
	my $assigns = "$GCd/$DB.assign.txt";
	my $cmd= "cat ".join(" ",@allFiles). " > $assigns\n";
	$cmd .= "rm -f ".join(" ",@allFiles) . "\n";
	#tr [:blank:] \\t
	$cmd .= "$rareBin sumMat -i $GCd/$countMatrixF.gz -o $GCd/$DB.mat -t 1 -refD $assigns $rtkFunDelims \n";
	print "@jdeps\n";
	$tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
	($cmd,$jdep) = qsubSystem($qsubDir."collect$DB.sh",$cmd,1,"40G","$DB"."Col",join(";",@jdeps),"",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	#return $jdep,$outF;
}



# minlengthpercent=0    (mlp) Smaller contig must be at least this percent of larger contig's len    gth to be absorbed.
# minoverlappercent=0   (mop) Overlap must be at least this percent of smaller contig's length to     cluster and merge.
# minoverlap=200        (mo) Overlap must be at least this long to cluster and merge.


# /g/bork3/home/hildebra/bin/bbmap/./dedupe.sh in=/g/scb/bork/hildebra/SNP/GCs/SimuB/B0/compl.fna out=/g/scb/bork/hildebra/SNP/GCs/SimuB/X0/ddtest.compl.fna exact=f threads=20 outd=/g/scb/bork/hildebra/SNP/GCs/SimuB/X0/drop.fna minidentity=95 storename=t renameclusters=t usejni=t cluster=t k=18 -Xmx50g


















