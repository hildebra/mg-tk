#!/usr/bin/env perl
#The Metagenomic Assembly, Genomic Recovery and Assembly Independent Mapping Tool (formerly MATAFILER, now mg-tk)
#main mg-tk routine
# (c) Falk Hildebrand, 2016-2024
#examples
#./mg-tk.pl map2tar test/refCtg.fasta,test/refCtg.fasta test1,test2
#./mg-tk.pl map2tar test/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa TEC2
#./mg-tk.pl -map dir/map

use warnings;
use strict;
use File::Basename;
use Cwd 'abs_path';
use POSIX;
use Getopt::Long qw( GetOptions );
use vars qw($CONFIG_FILE);


#load MF specific modules
use Mods::GenoMetaAss qw(readMap getDirsPerAssmblGrp lcp readFastHD prefixFAhd prefix_find gzipopen fileGZe
			readFasta writeFasta systemW getAssemblPath  filsizeMB resetAsGrps
			iniCleanSeqSetHR checkSeqTech is3rdGenSeqTech hasSuppRds 
			getRawSeqsAssmGrp getCleanSeqsAssmGrp addFileLocs2AssmGrp);
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile jgi_depth_cmd inputFmtSpades inputFmtMegahit createGapFillopt  
			buildMapperIdx mapperDBbuilt decideMapper  checkMapsDoneSH greaterComputeSpace);
use Mods::SNP qw(SNPconsensus_vcf );
use Mods::TamocFunc qw (cram2bsam getSpecificDBpaths getFileStr displayPOTUS bam2cram checkMF checkMFFInstall);
use Mods::phyloTools qw(fixHDs4Phylo);
#use Mods::Binning qw (runMetaBat runCheckM runSemiBin runMetaDecoder );
use Mods::Subm qw (qsubSystemWaitMaxJobs qsubSystem emptyQsubOpt findQsubSys qsubSystemJobAlive MFnext add2SampleDeps numUserJobs);


#some useful HPC commands..
#bjobs | awk '$3=="CDDB" {print $1}' |xargs bkill
#bjobs | grep 'dWXmOT' | cut -f12 -d' ' | xargs -t -i bkill {}
#squ | grep 'r' | cut -f12 -d' ' | xargs -t -i scontrol update TimeLimit=84:00:00 jobid={}
#squ | grep 'dencyNev' | cut -f12 -d' ' | xargs  -t -i scancel {}
#bhosts | cut -f1 -d' ' | grep -v HOST_NAME | xargs -t -i ssh {} 'killall -u hildebra'
#hosts=`bhosts | grep ok | cut -d" " -f 1 | grep compute | tr "\\n" ","`; pdsh -w $hosts "rm -rf /tmp/hildebra"



#local subs
sub announce_MGTK;
sub smplStats; sub checkDrives; 
sub isLastSampleInAssembly;
sub uploadRawFilePrep; sub unploadRawFilePostprocess;
sub seedUnzip2tmp; sub cleanInput; #unzipping reads; removing these at later stages ; remove tmp dirs
sub manageFiles;sub clean_tmp; 

sub sdmClean; sub sdmOptSet;  #qual filter reads
sub mergeReads; #merge reads via flash
sub removeHostSeqs; sub krakenTaxEst;sub prepKraken;
sub loop2C_check;

sub metagAssemblyRun;
sub createPsAssLongReads; #pseudo assembler
sub prepPreAssmbl; #hybrid pacbio/ill assemblies
sub genePredictions; sub run_prodigal; #gene prediction
sub mapReadsToRef;  sub scndMap2Genos;
sub runContigStats;#sub bam2cram;
sub mocat_reorder; sub postSubmQsub; 
sub detectRibo;  sub RiboMeta;sub riboSummary;

sub runOrthoPlacement;
sub runDiamond; sub DiaPostProcess; sub IsDiaRunFinished;
sub nopareil; sub calcCoverage;sub d2metaDist; 
sub metphlanMapping; sub mergeMP2Table;
sub mOTU2Mapping; sub mergeMotu2Table; sub prepMOTU2;
sub genoSize; sub check_map_done; sub check_depth_done;
sub postprocess;
sub setDefaultMFconfig;
sub getCmdLineOptions;
sub setupHPC;


#------- version history MATAFILER --------
#.22: fixing of assembly.txt paths, if output folders were later copied around
#.23: upgrade to metaphlan3 & motus3 & general update to workflow that these no longer need sdm filtered reads
#.24: new sdm version, better HDD space usage control, reworked EBIupload routines, sdm cores back to 4, bug fixes..
#.25: added SemiBin and MetaDecoder
#.26: added pprodigal, MUSCLE5, addapted cc.bin resource usage, semibin integration in MGS.pl, checkM2 by default used
#.27: added "loopTillComplete" functionality
#.28: extended loopTillComplete, better sq control
#.29: added more control of tmp node space, better reporting on scratch usage & paired samples
#.30: changed installer massively, now split into 3 packs: matafiler.yml, gtdbtk.yml, Rbase.yml. motus.yml is only for motus, as massive dep problems with this one..
#.31 further cleanup to installer, sdm does now trimmomatic job (saves 1 step), sdm v2.07, tmp_space refinements, SNP calling performance improved, flow improved
#.32: more automatation to ssd & RAM usage
#.33: fixes metagStats, phylotools, cc.bin and strainWithin.pl
#.34: cut5PR1 and cut5PR2 mapping headers added. sdm 2.10
#.35: bugfixes to flow control (should be less iteration MF now), replaced bedtools with mosdetph
#.36: rm mosdepth, replaced with samtools depth; integrated kma, but bowtie2 remains default aligner
#.37: bug fix rmDup correctly activated
#.38: first implementation of automatic sample lock; furthermore bam->fastq functionality (single end) implemented
#.39: smplLock improvements, genecat improvements, prevent rm bwt idx
#.40: refactor how internally file paths are transferred and kept up-to-date between different subfunctions
#.41: 31.1.24: bugfixes related to refactor. full implementation of eggNOGmapper. bugfixes in MGS and genecat workflows
#.42: 2.3.24: reworked how nodeTmpSpace & Ram are denoted internally; starting binner now in same round as other steps (one less round of starting MF)
#.43: 17.3.24: first version of hybrid ill-PB assemblies working! also reworked getProgPaths to be faster
#.44: 20.3.24: hybrid ill-PB implemented for assembly groups. PB assembler implemented for assembly group. added proper subfunction to manage assembly mode
#.45: 13.4.24: reactivated ribofinder module
#.46: 23.4.24: hybrid assembly reworked. refactor code. less false subm dependencies.
#.47: 24.4.24: gene coverage by default gz'd, more safety checks that related files are gz'd. semiBin2 by default used now. SemiBin2 can use support read coverage profiles in addition to primary reads (increases MAG recovery).
#.48: stabilized hybrid assm. Contigstats uses .gz as input.
#.49: 17.5.24: tentative strobealign added
#.50: 18.5.24: bugfixes workflow, autodetect to rm locks, code cleanups
#.51: code refactor
#.52: 6.6.24: start implementation (.3di). more checks on correct pipe excecution.
#,53: 27.7.24: routine to check if some essential programs have been correctly installed and are avaialble. added flag -checkInstall. Moved geneCat.pl and MGS.pl to secScripts/ dir to avoid confusion on where to start..
#.54: 15.8.24: renamed MATAFILER to mg-tk
#.55: 23.10.24: integrated -getAssemblConsSNPsuppRds -SNPconsMinDepth flags and functionalities
#.56: 25.10.24: small bugfixes to submission logic, to submit less jobs that would fail in any case
#.57:26.10.24: enabled multi-input files for minimap2/strobealign, added "-mapperLargeRef" flag
#.58: 30.12.24: firstXrdsRd & firstXrdsWr function added, required sdm 3.08
#.59: 21.1.25: hostile integration
#.60: 28.2.25: updates to loop2complete mechanic (debugging), clusterMAGs version 0.25
my $MATFILER_ver = 0.60;


#operation mode?
my $ARGV0 = "";$ARGV0 = $ARGV[0] if (@ARGV > 0);


#----------------- defaults ----------------- 


my %jmp=();
my $logDir = ""; #this is the local logdir
my $sharedTmpDirP = ""; #e.g. /scratch/MG-TK/
my $nodeTmpDirBase = "";#/tmp/MG-TK/
my $baseDir = ""; my $baseOut = "";


#control broad flow
my $FROM1=0; my $TO1=999999999999;
#counter on where MF is in map
my $JNUM=0;
my $doSubmit=1; 
my $rewrite=0;
my $baseID = "";


my $loop2completion = "0" ; #
my $loop2c_winsize=0;
my $loop2completion_ini=0;

#config is more overall configuration for MATAFILER
my %MFconfig; 

#MFcontstants: object to store essential paths/file endings
my %MFcontstants;

#MFopt: global object with options for MG-TK. Added in MF v0.5, slowly rebuild MF around this system
my %MFopt; 

#keep track of DBs that the metagenome will be filtered against..
my @filterHostDB = ();
#track secondary mapping and ref DBs
my %map2ndTogRefDB;my %make2ndMapDecoy;
my @bwt2outD =(); my @DBbtRefX = (); my @DBbtRefGFF=(); my @bwt2ndMapNmds;
my @scaffTarExternalOLib1; my @scaffTarExternalOLib2;

my @EBIjobs = (); #keeps track of $MFconfig{uploadRawRds} jobs, submits postprocessing (md5)
#----------- map all reads to a specific reference - options ---------

#progStats: object to track progress of programs/submissions
my %progStats;#count up progress of submitted jobs in current run

#HDDspace: object to handle HDD usage: Always format as "XXG" XX = space requirements in Gb. Excecption: "-1"
my %HDDspace;

my %inputFileSizeMB; #stores file size/sample
my %locStats; #keeps statistics of samples in hash (from already finished samples


#say hello to user 
announce_MGTK();
setDefaultMFconfig();
getCmdLineOptions;
checkMF(1);



#set up further dependencies for MF

if ($loop2completion =~ m/(\d+):(\d+)/){
	$loop2c_winsize = int($2);$loop2completion=$1;$loop2completion_ini=$1;
	print "Loop2completion=$loop2completion; Window size=$loop2c_winsize\n";
} elsif ($loop2completion ne "0") {
	$loop2completion = 6 ;$loop2completion_ini=6; #set to std number of iterations..
}


#dirs from config file--------------------------
#can be overwritten by $map{opt}{GlbTmpD} $map{opt}{NodeTmpD}
$sharedTmpDirP = getProgPaths("globalTmpDir",0) unless ($sharedTmpDirP ne "");
$nodeTmpDirBase = getProgPaths("nodeTmpDir",0) unless ($nodeTmpDirBase ne "");


#programs of global (pun) importance --------------------------
my $smtBin = getProgPaths("samtools");#
my $pigzBin  = getProgPaths("pigz");
my $avx2Constr =  getProgPaths("avx2_constraint",0);


#set up link to submission system on cluster
my $QSBoptHR = setupHPC();




# the map and some base parameters (base ID, in path, out path) can be (re)set
my %map; my %AsGrps; my %DOs;#DOs only required for metabat, to use all mappings within an assembly group
my ($hr,$hr2) = readMap($MFconfig{mapFile},0,\%map,\%AsGrps,$MFconfig{oldStylFolders});
%AsGrps = %{$hr2}; %map = %{$hr};
if ($MFopt{DoMetaBat2}){
	my ($hrD,$hrM) = getDirsPerAssmblGrp(\%map,\%AsGrps);
	%DOs = %{$hrD};
}

#$baseDir = $map{inDir} if (exists($map{inDir} ));
$baseOut = $map{opt}{outDir} if (exists($map{opt}{outDir} ) && $map{opt}{outDir}  ne "");
$baseID = $map{opt}{baseID} if (exists($map{opt}{baseID} ) && $map{opt}{baseID}  ne "");
die "provide an outdir in the mapping file\n" if ($baseOut eq "");
die "provide a baseID in the mapping file\n" if ($baseID eq "");
#overwrite tmp dirs??
if ($map{opt}{GlbTmpD} ne ""){print "Taking Global temp dir from map: $map{opt}{GlbTmpD} \n";$sharedTmpDirP = $map{opt}{GlbTmpD} ;}
if ($map{opt}{NodeTmpD} ne ""){print "Taking Node temp dir from map: $map{opt}{NodeTmpD} \n";$nodeTmpDirBase = $map{opt}{NodeTmpD} ;}


my $runTmpDirGlobal = "$sharedTmpDirP/$baseID/";
my $runTmpDBDirGlobal = "$runTmpDirGlobal/DB/";
unless (-d $runTmpDBDirGlobal){
	system "mkdir -p $runTmpDBDirGlobal" ;
	#and check that this dir exists...
	sleep(1);
	die "Can't create $runTmpDBDirGlobal\n" unless (-d $runTmpDBDirGlobal);
}
my $globaldDiaDBdir = $runTmpDBDirGlobal."DiamDB/";
system "mkdir -p $globaldDiaDBdir" unless (-d $globaldDiaDBdir);


#----------- map all reads to a specific reference, preparation ---------
my $map2ndMpde=0;#0=map2tar;2=map2DB;3=map2GC
#-----------   scaffolding external contigs parameters
my $scaffTarExternal = "";my $scaffTarExternalName = ""; 
my $scaffTarExtLibTar = ""; my $bwt2ndMapDep = ""; 

map2ndPrep();



#some base stats kept in vars
my $sequencer = "hiSeq";#plattform the algos have to deal with
#my $DBpath=""; 
my $assDir=""; 
my $continue_JNUM = 0; #debug, set to 0 for full run
my $prevAssembly = ""; my $shortAssembly = "";#files with full length and short length assemblies
#my $mmpuOutTab = "";
#fixed stones
my $STOpreAssmblDone = "preassmblDone.sto"; #marks preAssebmly (hybrid assemblies etc) is done
my $STOassmbleDone = "ass.done.sto"; #marks assembly is done


#fixed dirs for specific set of samples
my $dir_MP2 = $baseOut."pseudoGC/Phylo/MP2/"; #metaphlan 2 dir
my $dir_mOTU2 = $baseOut."pseudoGC/Phylo/mOTU2/"; #mOUT 2 dir
my $dir_TaxTar = $baseOut."pseudoGC/Phylo/TaxaTarget/"; #taxaTar dir
my $dir_ContigStats = "/assemblies/metag/ContigStats/";


my $dir_RibFind = $baseOut."pseudoGC/Phylo/RiboFind/"; #ribofinder dir
my $dir_KrakFind = $baseOut."pseudoGC/Phylo/KrakenTax/$MFopt{globalKraTaxkDB}/"; #kraken dir
system("mkdir -p $baseOut") unless (-d $baseOut);
my $globalLogDir = $baseOut."LOGandSUB/"; #this is the gloabl logdir (across all samples in current run)
system("mkdir -p $globalLogDir/sdm") unless (-d "$globalLogDir/sdm");
open $QSBoptHR->{LOG},">",$globalLogDir."qsub.log";# unless ($doSubmit == 0);
my $collectFinished = $baseOut."runFinished.log\n";
my $globalNPD = $baseOut."NonPareil/";
system "mkdir -p $globalNPD" unless (-d "$globalNPD");
system "cp $MFconfig{mapFile} $globalLogDir/inmap.txt";

print $globalLogDir."qsub.log\n";

my $presentAssemblies = 0; my $totalChecked=0;
my @samples = @{$map{opt}{smpl_order}}; my @allSmplNames;
my @allFilter1; my @allFilter2; my @inputRawFQs; 
my @EmptySample;

my $statStr = ""; my $statStr5 = "";
my %sampleSDMs; 
my $curSmpl = "";
#my %assmblGrpLog; #security that no assmblGrp exists twice.. #actually not needed..


my $baseSDMopt = getProgPaths("baseSDMopt_rel"); #"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_inifilter_relaxed.txt";
if ($MFopt{useSDM} ==2 ){$baseSDMopt = getProgPaths("baseSDMopt");}#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_inifilter.txt";}
my $baseSDMoptMiSeq = getProgPaths("baseSDMoptMiSeq_rel");#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_miSeq.txt";	
if ($MFopt{useSDM} ==2 ){$baseSDMoptMiSeq = getProgPaths("baseSDMoptMiSeq");}#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/data/sdm_opt_miSeq_relaxed.txt";	}

my @unzipjobs; 
my $sdmjNamesAll = "";
my $waitTime = 0;
my @grandDeps; #used for loop2completion , collects dependencies 

my $emptCmd = "sleep 333";
#qsubSystem($globalLogDir."emptyrun.sh","sleep 333",1,"1G",0,"empty","","",1);


#set up kraken human filter
my $krakDeps = ""; my $krakenDBDirGlobal = $runTmpDirGlobal;
if ($MFopt{DoKraken} && $MFopt{globalKraTaxkDB} eq ""){die "Kraken tax specified, but no DB specified\n";}
$krakDeps = prepKraken() if ( ($MFopt{humanFilter}>0 && $MFopt{humanFilter}<3) || ($MFopt{DoKraken}) || $MFopt{DoEukGenePred});
#profiling prep
my $mOTU2Deps = prepMOTU2(); prepMetaphlan();
#redo d2s intersample distance?
if ($MFopt{DoCalcD2s}) {$MFopt{DoCalcD2s} = !-e "$baseOut/d2StarComp/d2meta.stone";}



if ($TO1 > @samples){
	print "Reset range of samples to ". @samples."\n";
	$TO1 = @samples;
}
my $from = $FROM1; my $to = $TO1;
#die "\"@samples\"\n";
if ($loop2c_winsize > 0){
	$to = $from + $loop2c_winsize; $to = $TO1 if ($to > $TO1);
}


#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#die $to."\n";
#for loop that goes over every single sample in the .map
for ($JNUM=$from; $JNUM<$to;$JNUM++){
	
	qsubSystemWaitMaxJobs($MFconfig{checkMaxNumJobs}, $MFconfig{killDepNever});
	
	
	#set up initial local paths for a given sample
	my $dir2rd=""; my $curDir = "";my $curOutDir = ""; 
	$curSmpl = $samples[$JNUM];
	#key IDs for sample
	my $cAssGrp = $curSmpl;my $cMapGrp = $map{$curSmpl}{MapGroup};
	if ($map{$curSmpl}{AssGroup} ne "-1"){ $cAssGrp = $map{$curSmpl}{AssGroup};}
	my @sampleDeps = (); #catalogues all job dependencies created in this loop

	#local flow control
	if (exists($jmp{$JNUM})){loop2C_check($cAssGrp,\@sampleDeps);next;}

	#print "SMPL::$curSmpl\n";
	$dir2rd = $map{$curSmpl}{dir};
	$dir2rd = $map{$curSmpl}{prefix} if ($dir2rd eq "");
	#print "$curSmpl  $map{$curSmpl}{dir}  $map{$curSmpl}{rddir}  $map{$curSmpl}{wrdir}\n";	next;
	my $SmplName = $map{$curSmpl}{SmplID};
	if ($dir2rd eq "" ){#very specific read dir..
		if ($map{$curSmpl}{SupportReads} ne ""){
			$curDir = "";	$curOutDir = "$baseOut$SmplName/";	
		} else {
			die "Can;t find valid path for $SmplName\n";
		}
		$dir2rd = $SmplName;
	} else {
		$curDir = $map{$curSmpl}{rddir};	$curOutDir = $map{$curSmpl}{wrdir};	
	}
	$logDir = "$curOutDir/LOGandSUB/";
	system("mkdir -p $logDir") unless (-d $logDir);

	#ignore samples .. for various reasons ------------------------------------------------------------------------------------
	if ($MFconfig{ignoreSmpl} ne ""){
		if ($MFconfig{ignoreSmpl} =~ m/$SmplName/){print "\n ======= Ignoring sample $SmplName =======\n";
		loop2C_check($cAssGrp,\@sampleDeps);next;}
	}
	my $smplLockF = "$logDir/$MFcontstants{DefaultSampleLock}";
	if (-e $smplLockF){
		if ($MFconfig{rmSmplLocks}){ system "rm -f $smplLockF";
		} else {print "\n    >>>>>>>>>> Sample $SmplName is locked! <<<<<<<<<<  \n";loop2C_check($cAssGrp,\@sampleDeps);next;}
	} 
	$QSBoptHR->{LOCKfile} = $smplLockF; #set lockfile to be created if any job is submitted
	

	push (@allSmplNames,$SmplName);
	#die "$curDir\n$curOutDir\n$baseOut\n";
	
	my $samplReadLength = $MFconfig{defaultReadLength}; #some default value
	if (exists $map{$curSmpl}{readLength} && $map{$curSmpl}{readLength} != 0){
		$samplReadLength = $map{$curSmpl}{readLength};
	}
	my $samplReadLengthX = $MFconfig{defaultReadLengthX}; #for any supplementary reads (eg PacBio)
	if (exists $map{$curSmpl}{readLengthX} && $map{$curSmpl}{readLengthX} != 0){
		$samplReadLengthX = $map{$curSmpl}{readLengthX};
	}
	
	
	my $contigStatsUsed = 0;

	%locStats = ();
	$locStats{hasPaired} = 0;	$locStats{hasSingle} = 0; 

	$totalChecked++;
	
	#die $SmplName."  $samplReadLength\n";
	print "\n======= $SmplName - $JNUM - $dir2rd =======\n" unless($MFconfig{silent});
	
	
	
	#set up dirs ------------------------------------------------------------------------------------
	my $smplTmpDir = "$runTmpDirGlobal$SmplName/"; #curTmpDir
	my $nodeSpTmpD = "$nodeTmpDirBase/$SmplName";
	my @checkLocs = ($smplTmpDir);
	push(@checkLocs,$curDir)  if ($curDir ne "");
	$QSBoptHR->{LocationCheckStrg} = checkDrives(\@checkLocs);

	my $mapOut = "$smplTmpDir/mapping/";
	#$DBpath="$curOutDir/readDB/";
	my $finalCommAssDir = "$curOutDir/assemblies/metag/";
	my $finalCommAssDirSingle = $finalCommAssDir; #this is only used for checking..
	my $finalMapDir = "$curOutDir/mapping/";
	my $KrakenOD = $curOutDir."Tax/kraken/$MFopt{globalKraTaxkDB}/";
	
	
	$AsGrps{$cAssGrp}{AssemblSmplDirs} .= $curOutDir."\n";
	my $AssemblyGo=0; my $MappingGo=0; #controls if assemblies / mappings are done in respective groups
	
	#complicated flow control for multi sample assemblies
	die "cAssGrp eq  \"$cAssGrp\" ".$baseID."\n" if ($cAssGrp eq "");
	my $assmGrpTag = "AssmblGrp_$cAssGrp";
	#die "AssemblGrp exists already!!: \"$assmGrpTag\"\n" if (exists($assmblGrpLog{$assmGrpTag}));
	#$assmblGrpLog{$assmGrpTag} = 1;
	$assDir="$runTmpDirGlobal/$assmGrpTag/";
	$finalCommAssDir = "$baseOut/$assmGrpTag/metag/" if ( !exists($AsGrps{$cAssGrp}{CntAimAss}) || $AsGrps{$cAssGrp}{CntAimAss}>1);
	#my $metaGpreAssmblDir = "$runTmpDirGlobal/pre$assmGrpTag/";
	my $metaGpreAssmblDir = "$smplTmpDir/pre$assmGrpTag/"; #needs to be sample specific (to allow for different coverage)

	
	
	#assign job name (dependency) only ONCE
	if ( !exists($AsGrps{$cAssGrp}{CntAss}) || $AsGrps{$cAssGrp}{CntAss} == 0){
		$AsGrps{$cAssGrp}{AssemblJobName} = "";#"_XXASpl$cAssGrp"."XX_" ;
		$AsGrps{$cAssGrp}{CSfinJobName} = "";#"_XXCSpl$cAssGrp"."XX_" ;
	}
	
	$AsGrps{$cAssGrp}{CntAss} ++;
	print "AssmblyGrp: " . $AsGrps{$cAssGrp}{CntAss} .":".$AsGrps{$cAssGrp}{CntAimAss}. "; " if ($AsGrps{$cAssGrp}{CntAimAss} > 1 && !$MFconfig{silent});
	if ($AsGrps{$cAssGrp}{CntAss}  >= $AsGrps{$cAssGrp}{CntAimAss} ){
		#print "running assembluy\n";
		$AssemblyGo = 1;
		if ($AsGrps{$cAssGrp}{CntAimAss}<=1){
			$assDir="$smplTmpDir/assemblies/";
		}
	}
	
	
	#mapping groups?
	$AsGrps{$cMapGrp}{CntMap} ++;
	print "MapGroup: ".$AsGrps{$cMapGrp}{CntMap} .":".$AsGrps{$cMapGrp}{CntAimMap}.";\n"  unless ($MFconfig{silent});
	if (!exists($AsGrps{$cMapGrp}{CntMap})){ die "Can;t find CntMap for $cMapGrp";}
	if ($AsGrps{$cMapGrp}{CntMap}  >= $AsGrps{$cMapGrp}{CntAimMap}  ) {
		#print "running mapping";
		$MappingGo = 1;
	}
	#die "MG: $MappingGo\n";
	if ($MFopt{DoAssembly} ==0 ){$AssemblyGo=0;}
		#mapping related
	my $cramthebam = 1;
	my $bamcramMap = "bam"; if ($cramthebam){$bamcramMap = "cram";}

	
	#----------------------------  2  ----------------------------------
	#set up dirs for this sample
	my $metagAssDir = $assDir."metag/";
	my $geneDir = $metagAssDir."genePred/";
	#system("mkdir -p $metagAssDir $geneDir");
	my $metaGassembly=$metagAssDir."scaffolds.fasta.filt"; 
	my $finalCommScaffDir = "$finalCommAssDir/scaffolds/";
	my $metaGscaffDir = "$metagAssDir/scaffolds/";
	my $STOfinScaff = "$finalCommScaffDir/scaffDone.sto";
	my $pseudoAssFile = "$metagAssDir/longReads.fasta.filt";
	my $pseudoAssFileFinal = "$finalCommAssDir/longReads.fasta.filt";
	my $finAssLoc = "$finalCommAssDir/scaffolds.fasta.filt";
	my $ContigStatsDir  = "$curOutDir/$dir_ContigStats/";
	my $coveragePerCtg = "$ContigStatsDir/Coverage.percontig.gz";
	my $suppCoveragePerCtg = "$ContigStatsDir/Cov.sup.percontig.gz";
	my $nonParDir = $curOutDir."nonpareil/";
	#SNP calling on assembly related files
	my $SNPdir = "$curOutDir/SNP/";
	my $contigsSNP = "$SNPdir/contig.SNPc.$MFopt{SNPcallerFlag}.fna"; #keep without gz, although will be gz'd
	my $genePredSNP = "$SNPdir/genes.shrtHD.SNPc.$MFopt{SNPcallerFlag}.fna.gz";
	my $genePredAASNP = "$SNPdir/proteins.shrtHD.SNPc.$MFopt{SNPcallerFlag}.faa.gz";
	my $vcfSNP = "$SNPdir/allSNP.$MFopt{SNPcallerFlag}.vcf";	$vcfSNP = "" if (!$MFopt{saveVCF});
	my $vcfSNPsupp = "$SNPdir/allSNP.$MFopt{SNPcallerFlag}-sup.vcf";	$vcfSNPsupp = "" if (!$MFopt{saveVCF});
	my $CRAMmap = "$finalMapDir/$SmplName-smd.cram";
	my $SupCRAMmap = "$finalMapDir/$SmplName.sup-smd.cram";
	my $inputRawFile = "$curOutDir/input_raw.txt";
	my $binningDir = "$finalCommAssDir/Binning/";
	my @smplIDs = ("");#smpls in current assembly group
	if ($MFopt{DoMetaBat2}){ 
		@smplIDs = @{$DOs{$cAssGrp}{SmplID}};
	}
	my $BinningOut = "$binningDir/MB2/$smplIDs[-1]";
	$BinningOut = "$binningDir/SB/$smplIDs[-1]" if ($MFopt{DoMetaBat2} == 2);
	$BinningOut = "$binningDir/MD/$smplIDs[-1]" if ($MFopt{DoMetaBat2} == 3);
	
	
	#stones
	my $STOcram = "$CRAMmap.sto";
	my $STOsupCram = "$SupCRAMmap.sto";
	my $STOmapFinal="$finalMapDir/done.sto";
	my $STOsnpCons = "$contigsSNP.SNP.cons.stone"; #	my $SNPstone = $ofasConsDir."SNP.cons.stone";
	my $STOsnpSuppCons = "$contigsSNP.SNP.supp.cons.stone";
	
	# collect stats on seq qual, assembly etc
	if ($MFconfig{alwaysDoStats}){
		my ($statsHD,$curStats,$statsHD5,$curStats5) = smplStats($curOutDir,$assDir,$SmplName);
		#add to global string that is later written
		if ($statStr eq ""){
			$statStr.="SMPLID\tDIR\t".$statsHD."\n".$curSmpl."\t$dir2rd\t".$curStats."\n";
			$statStr5.="SMPLID\tDIR\t".$statsHD5."\n".$curSmpl."\t$dir2rd\t".$curStats5."\n";
		} else {$statStr.=$curSmpl."\t$dir2rd\t".$curStats."\n"; $statStr5.=$curSmpl."\t$dir2rd\t".$curStats5."\n";}
	}
	#die;
	
	#detect what already exists..
	my $efinAssLoc = 0; $efinAssLoc = 1  if (-s $finAssLoc && -e "$finalCommAssDir/$STOassmbleDone");
	#activate if two assemblies for single sample required, e.g. hybrid assemblies
	my $ePreAssmbly = 0; $ePreAssmbly = 1 if (-s $finAssLoc && -e "$finalCommAssDir/$STOpreAssmblDone");
	my $ePreAssmblPck = 0; $ePreAssmblPck = 1 if (-s $finAssLoc && -e "$metaGpreAssmblDir/moved.sto");
	my $doPreAssmFlag = 0; my $postPreAssmblGo =0 ;
	($doPreAssmFlag,$postPreAssmblGo,$ePreAssmblPck) = prepPreAssmbl($finalCommAssDir,$metaGpreAssmblDir,$finalMapDir, "$smplTmpDir/preAssmblData/",
				$ContigStatsDir, $curSmpl, $cAssGrp, $efinAssLoc, $ePreAssmblPck);#moves files to new locations
				
	my $eCovAsssembly = 0; $eCovAsssembly = 1 if (-e $coveragePerCtg);
	my $eSuppCovAsssembly = 0; $eSuppCovAsssembly = 1 if (-e $suppCoveragePerCtg);
	#will be created in contigstats step (not related to bowtie & sortbam)
	my $emetaGassembly = 0; $emetaGassembly = 1 if (-e $metaGassembly && -e "$metagAssDir/$STOassmbleDone"); 
	my $eFinMapCovGZ = 0; $eFinMapCovGZ = 1 if (-e $STOcram && -e "$finalMapDir/$SmplName-smd.bam.coverage.gz");#"$finalMapDir/$SmplName-smd.bam.coverage.gz"; 
	#$MFopt{mapSupport2Assembly}
	my $locMapSup2Assembly =0; $locMapSup2Assembly =1 if ($MFopt{mapSupport2Assembly} && $map{$curSmpl}{"SupportReads"} ne "");
	my $eFinSupMapCovGZ = 0; $eFinSupMapCovGZ = 1 if ($locMapSup2Assembly && -e $STOsupCram && -e "$finalMapDir/$SmplName.sup-smd.bam.coverage.gz");
	#die "$eFinSupMapCovGZ  $finalMapDir\n";
	my $dfinalCommAssDir = 0 ; $dfinalCommAssDir = 1 if (-d $finalCommAssDir);
	my $eFinalMapDir = 0; $eFinalMapDir = 1 if (-s $STOmapFinal);
	#upload2EBI 
	my $DoUploadRawReads = 0; $DoUploadRawReads = 1 if ($MFconfig{uploadRawRds} ne ""); 
	
	
	#print "$emetaGassembly   XX";
	#die "$eFinMapCovGZ $STOcram && -e $finalMapDir/$SmplName-smd.bam.coverage.gz\n";

	#first take care of moving finished assemblies..
	if (!$efinAssLoc && $emetaGassembly){
		print "moving finished assembly to final location\n";
		#print "$metaGassembly\n";
		systemW "mkdir -p $finalCommAssDir" unless ($dfinalCommAssDir);
		systemW "rsync -r  --remove-source-files $metagAssDir/* $finalCommAssDir/";
		$efinAssLoc = 1; $emetaGassembly = 0;
	}
	

	
	
	my $locRedoSNPcalling =0; 
	my $locRedoAssMapping = $MFopt{redoAssMapping};

	#check if current assembly group is the same as before!
	my $locRewrite = 0; my $locRedoAssembl = 0;
	if ($efinAssLoc && -e "$finalCommAssDir/smpls_used.txt"){
		my $currAssmlCnt = 0;#`cat $finalCommAssDir/smpls_used.txt | grep -v '^\\\$' | wc -l`;
		open I,"<$finalCommAssDir/smpls_used.txt " or die "Smple used: $!\n"; while (<I>){$currAssmlCnt++ unless m/^$/;} close I;
		chomp $currAssmlCnt; $currAssmlCnt = int($currAssmlCnt);
		if ($currAssmlCnt < $AsGrps{$cAssGrp}{CntAimAss}){
			print "$cAssGrp assembl count has changed! (from $currAssmlCnt to $AsGrps{$cAssGrp}{CntAimAss})\n$finalCommAssDir\nRemoving assembly and all processed reads\n";
			unless ($MFconfig{OKtoRWassGrps}) {print "Stopping MATAFILER, human intervention needed.. use the flag \"-OKtoRWassGrps 1\" to allow MATAFILER to delete files\n"; die;}
			$locRewrite=1;
		}
	}
	
	if ($eFinMapCovGZ && (!$ePreAssmblPck && !$ePreAssmbly && !$efinAssLoc && !$emetaGassembly) ){#impossible, so reason must be severe!
		print "Mapping exists, but no assembly, removing mapping..\n$finalMapDir\n$finAssLoc\n$metaGassembly\n";
		die unless ($MFconfig{OKtoRWassGrps});
		$locRewrite = 0; $locRedoAssembl = 0;
	}
	
	if ($efinAssLoc && $finAssLoc ne "$finalCommAssDirSingle/scaffolds.fasta.filt" && -s "$finalCommAssDirSingle/scaffolds.fasta.filt"){
		print "Something wrong.. assembly group assembly and single assembly present:\n$finAssLoc\n$finalCommAssDirSingle/scaffolds.fasta.filt\n";
		#die;
		$locRedoAssMapping=1;$locRedoSNPcalling=1;
		system "rm -fr $finalCommAssDirSingle; mkdir -p $finalCommAssDirSingle;\n";
		$eCovAsssembly=0;$eFinMapCovGZ=0;$eFinalMapDir=0;$eFinSupMapCovGZ=0;$eSuppCovAsssembly=0;$eSuppCovAsssembly=0;
	}
	
	if ($locStats{totRds}==0 && !defined($locStats{uniqAlign}) && -e $inputRawFile && $eFinMapCovGZ){ #do a deeper look
		my $line = getFileStr("$inputRawFile",0); #open I,"<$inputRawFile"; my $line = <I>; close I; chomp($line);
		#my @spl = split /,/,$line; my $inFileSize = -s $spl[0];
		print "weird empty: $locStats{totRds} $line \nredoing..\n";
		die;
		$locRewrite = 1;$locRedoAssembl = 1;
	}
	if ( ($MFconfig{skipWrongPairedSmpls} || $MFconfig{OKtoRWassGrps}) && -e "$logDir/sdmReadCleaner.sh.etxt" && `tail -n 70 $logDir/sdmReadCleaner.sh.etxt | grep 'invalid paired read' ` ne ""){
		print "$logDir/sdmReadCleaner.sh.etxt problems! Delete outdir\n";
		if ($MFconfig{OKtoRWassGrps}){
			$locRewrite=1 ;
		}elsif ($MFconfig{skipWrongPairedSmpls}){
			loop2C_check($cAssGrp,\@sampleDeps);next;
		}
	}
	
	
	#--------------------------  DELETION SECTION  -----------------------------------
	#DELETION SECTION
	#redo run - or parts thereof	
		
	if ($MFconfig{OKtoRWassGrps} && ($rewrite || $locRewrite)){
		print "Deleting previous results.. rerun MATAFILER for sample\n";
		system ("rm -r -f $assDir $finalCommAssDir");
		system("rm -f -r $curOutDir $smplTmpDir $collectFinished ");
		#next; #too deep, needs a complete new round over dir..
		#$efinAssLoc = 0;	$eFinMapCovGZ = 0;	$emetaGassembly =  0;
		$efinAssLoc =0 ;$emetaGassembly = 0;  $dfinalCommAssDir =0;
		$eCovAsssembly = 0; $eSuppCovAsssembly=0; $eFinSupMapCovGZ=0; $eFinMapCovGZ = 0;$eFinalMapDir = 0;$locRewrite = 0; $locRedoAssembl = 0;$eSuppCovAsssembly=0;
	} 
	
	#automatically delete mapping, if assembly no longer exists..
	#print "locRedoAssMapping : $locRedoAssMapping\n";
	if ($MFopt{map2Assembly} ){
		if ($eFinMapCovGZ && !$eFinalMapDir){$locRedoAssMapping = 1 ; print "R0 ";}
		if (!$efinAssLoc && !$ePreAssmbly && !$emetaGassembly ){$locRedoAssMapping = 1 ;}#print "R1";}
		if (!$MappingGo && $eFinalMapDir){$locRedoAssMapping = 1 ;print "R2 ";}
		if (-e $STOcram && !-e "$finalMapDir/$SmplName-smd.bam.coverage.gz" && !-e "$finalMapDir/$SmplName-smd.bam.coverage"){$locRedoAssMapping = 1 ;print "R3 ";}
		#if ($eFinMapCovGZ && (exists($locStats{uniqAlign}) && $locStats{uniqAlign} > 20) && -s $CRAMmap <300){$locRedoAssMapping = 1 ;print "R4";}
		#print "$CRAMmap :: $locRedoAssMapping\n";
		print "redo assem mapping!" . " -s $CRAMmap \n" if ($locRedoAssMapping && -e $CRAMmap); #. (-s $CRAMmap)
		
		#die "$STOcram && !-e $finalMapDir/$SmplName-smd.bam.coverage.gz";
		
	}

	#die "locRedoAssMapping : $locRedoAssMapping $finalMapDir   : !$efinAssLoc && !$emetaGassembly\n" if ($locRedoAssMapping);
	#my $sizemap = -s $CRAMmap;#print "size map: " . $sizemap . "\n";
	#delete assembly
	if ($MFopt{redoAssembly} || $locRedoAssembl){
		print "Removing assembly ... \n" if ($emetaGassembly);
		system "rm -fr $finalCommAssDir";
		$efinAssLoc = 0;	
		$locRedoAssMapping=1;
	}
	#delete mapping to assembly
	if ($locRedoAssMapping){
		#die "locDel\n";
		print "Deleting previous assembly mapping, size map: ". (-s $CRAMmap) . "\n$CRAMmap\n" if ($MappingGo && $eFinalMapDir);
		#die "$finAssLoc && !-e $metaGassembly\n";
		system "rm -fr $finalMapDir $mapOut $ContigStatsDir/Coverage.* $ContigStatsDir/Cov.sup.*";
		$eFinMapCovGZ = 0;	$emetaGassembly =  0;
		$eCovAsssembly = 0; $eSuppCovAsssembly=0; $eFinSupMapCovGZ=0; $eFinalMapDir = 0;
		#are there SNPs called? remove as well..
		$locRedoSNPcalling=1; 
	}
	#Case: primary assembly mapping was done, support reads were not yet mapped.. need to redo binning 
	#print "$locMapSup2Assembly && !$eFinSupMapCovGZ) && ($MFopt{map2Assembly} && $eFinMapCovGZ \n";
	if ( ($locMapSup2Assembly && !$eFinSupMapCovGZ) && ($MFopt{map2Assembly} && $eFinMapCovGZ ) ){
		print "redoing binning due to support mapping not included..\n";
		system("rm -rf  $binningDir/");
	}
	#debug case: binning was empty
	if ($MFopt{DoMetaBat2} && ( $MFopt{BinnerRedoAll} || ($MFopt{BinnerRedoEmpty} && -e $BinningOut && !-s $BinningOut) ) ){
		#die;
		print "redoing binning due to empty bins (flag -redoEmptyBins 1) ..\n";
		system "rm -rf $binningDir";
	}

	if ((!$eFinMapCovGZ && $eCovAsssembly)  #redo only contigstats related to coverage..
				|| ($eSuppCovAsssembly && !$eFinSupMapCovGZ) || $MFconfig{redoCS}){
		print "redoing contig stats global..\n";
		system("rm -rf $finalCommAssDir/ContigStats/ $ContigStatsDir $binningDir/");
		$eCovAsssembly = 0; $eSuppCovAsssembly=0; #contigstats needs redoing..
	}
	system "rm -r $KrakenOD" if ($MFopt{RedoKraken} && -d $KrakenOD);
	if ($MFopt{RedoRiboFind}){system "rm -rf $curOutDir/ribos";}
	if ($MFopt{RedoRiboAssign}){system "rm -rf $curOutDir/ribos//ltsLCA";}
	if ($MFopt{DoRibofind} && -e "$curOutDir/LOGandSUB/RiboLCA.sh.etxt"){
		#my $LCAetxt = `cat $curOutDir/LOGandSUB/RiboLCA.sh.etxt`;
		#if ($LCAetxt =~ m/ParseError thrown: Unexpected character .\@. found/){system "rm -rf $curOutDir/ribos";}
	}
	if ($locRedoSNPcalling){system "rm -fr $SNPdir";}
	if ($MFopt{redoSNPcons}){		system "rm -rf $genePredSNP* $contigsSNP* $genePredAASNP* $smplTmpDir/SNP $logDir/SNP";
	} elsif ($MFopt{redoSNPgene}){		system "rm -rf $genePredSNP* $genePredAASNP* ";
	}
	my $boolGenePredOK=0;
	if ($MFopt{DoEukGenePred}){
		$boolGenePredOK = 1 if (-s "$finalCommAssDir/genePred/proteins.bac.shrtHD.faa" || ($MFopt{pseudoAssembly} && -e "$finalCommAssDir/genePred/proteins.bac.shrtHD.faa"));
	} else {
		$boolGenePredOK = 1 if (-s "$finalCommAssDir/genePred/proteins.shrtHD.faa" || ($MFopt{pseudoAssembly} && -e "$finalCommAssDir/genePred/proteins.shrtHD.faa") );
	}
	#central flag
	my $boolAssemblyOK=0;
	$boolAssemblyOK=1 if ($boolGenePredOK && $efinAssLoc );#&& (!$MFopt{map2Assembly} || $eFinMapCovGZ ) );
	#die "$boolGenePredOK && $efinAssLoc && (!$MFopt{map2Assembly} || $eFinMapCovGZ ) $locRedoAssMapping\n";
			#&& (-s "$finalMapDir/$SmplName-smd.bam" || -s "$finalMapDir/$SmplName-smd.cram")
	#die "$boolAssemblyOK\n$finalCommAssDir/genePred/proteins.shrtHD.faa\n$finalMapDir/$SmplName-smd.bam.coverage.gz\n";
	

	if ( ( !$boolAssemblyOK && $MFconfig{unfiniRew}==1 ) ){
		die "Deleting previous results..\n";
		system("rm -f -r $curOutDir $smplTmpDir $collectFinished");
		$efinAssLoc = 0;	$eFinMapCovGZ = 0;	$emetaGassembly =  0;
	} elsif ($boolAssemblyOK && $eCovAsssembly) {
		#check that assembly path fits..
		getAssemblPath($curOutDir,$finalCommAssDir);
	}

	
#--------------------- secondary map deletions & flags --------------------------
	my $boolScndMappingOK = 1; my $iix =0;
	my $boolScndCoverageOK = 1;
	if ($MFopt{MapRewrite2nd}){ 
		print "rewriting secondary map\n";
		foreach my $bwt2outDTT (@bwt2outD){
			my $expectedMapCovGZ = "$bwt2outDTT/$bwt2ndMapNmds[$iix]"."_".$SmplName."-0-smd.bam.coverage.gz"; #$bamcramMap : 2nd map only has .bam output
			system "rm -f $expectedMapCovGZ*";
			$eFinMapCovGZ = 0;	
		}
	}
	#die "eFinMapCovGZ $eFinMapCovGZ\n";
	foreach my $bwt2outDTT (@bwt2outD){
		my $expectedMapCovGZ = "$bwt2outDTT/$bwt2ndMapNmds[$iix]"."_".$SmplName."-0-smd.bam.coverage.gz";
		my $expectedMapBam = "$bwt2outDTT/$bwt2ndMapNmds[$iix]"."_".$SmplName."-0-smd.bam";
		$iix++;
		
		#print $expectedMapCovGZ."\n";
		if ( -e "$expectedMapCovGZ" && -e $expectedMapBam && $MappingGo  ){
			$boolScndMappingOK=1;
		}else{
			$boolScndMappingOK=0; $boolScndCoverageOK=0;
			#be clean
			system "rm -f $expectedMapBam*";
			last;
		}
		if ($MFopt{mapModeCovDo} && (!-e $expectedMapCovGZ.".median.percontig" || !-e $expectedMapCovGZ.".percontig"|| !-e $expectedMapCovGZ.".pergene")){
			$boolScndCoverageOK=0;
			system "rm -f $expectedMapCovGZ.*";
		}
	}
	if (@bwt2outD == 0 ){$boolScndMappingOK = 1 ; $boolScndCoverageOK=1;}#|| !$MappingGo);	
	#die "$boolScndMappingOK\n";
	#print $boolScndMappingOK."\n$boolAssemblyOK\n";

#--------------------- other flags --------------------------
	#contamination flag: redo/do contaminant removal to eg make sure human contamination is correctly logged
	my $calcContamination = 0;
	$calcContamination = 1 if ($locStats{contamination} eq "?\t" && $efinAssLoc && $MFopt{completeContaStats});
	#print "Conta: $calcContamination   \"$locStats{contamination}\"\n";

	#Kraken flag
	my $calcKraken =0;
	$calcKraken = 1 if ($MFopt{DoKraken} && (!-d $KrakenOD || !-e "$KrakenOD/krakDone.sto"));
	if (!$calcKraken && $MFopt{DoKraken}){
		opendir D, $KrakenOD; my @krkF = grep {/krak\./} readdir(D); closedir D;
		foreach my $kf (@krkF){
			$kf =~ m/krak\.(.*)\.cnt\.tax/; my $thr = $1;# die $thr."  $kf\n";
			system "mkdir -p $dir_KrakFind/$thr" unless (-d "$dir_KrakFind/$thr"); #system "mkdir -p $dir_RibFind/SSU/" unless (-d "$dir_RibFind/SSU/"); system "mkdir -p $dir_RibFind/LSU/" unless (-d "$dir_RibFind/LSU/");
			system "cp $KrakenOD/$kf $dir_KrakFind/$thr/$SmplName.$thr.krak.txt";
		}
	} else {$progStats{KrakTaxFailCnts}++;}
	
	#system "rm -f $curOutDir/ribos//ltsLCA/LSU_ass.sto $curOutDir/ribos//ltsLCA/Assigned.sto"; fix for new LSU assignments
	my $calcGenoSize=0; $calcGenoSize=1 if ($MFopt{DoGenoSizeEst} && 	!-e "$curOutDir/MicroCens/MC.0.result");
	my $calcRibofind = 0; my $calcRiboAssign = 0;
	$calcRibofind = 1 if ($MFopt{DoRibofind} && (!-e "$curOutDir/ribos//SSU_pull.sto"|| !-e "$curOutDir/ribos//LSU_pull.sto" || ($MFopt{doRiboAssembl} && !-e "$curOutDir/ribos/Ass/allAss.sto" ))); #!-e "$curOutDir/ribos//ITS_pull.sto"|| 
	$calcRiboAssign = 1 if ($MFopt{DoRibofind} && ( !-e "$curOutDir/ribos//ltsLCA/Assigned.sto"  || !-e "$curOutDir/ribos//ltsLCA/SSU_ass.sto") );		#!-e "$curOutDir/ribos//ltsLCA/ITS_ass.sto"||  #ITS no longer required.. unreliable imo
	#die "$calcRibofind\n\n";
	#die "$calcRibofind $calcRiboAssign\n";
	RiboMeta($calcRibofind,$calcRiboAssign,$curOutDir,$SmplName);

	my ($calcDiamond,$calcDiaParse) = IsDiaRunFinished($curOutDir);
	#die "XX $calcDiamond $calcDiaParse\n";
	my $calcMetaPhlan=0;
	#die $dir_MP2."$SmplName.MP2.sto";
	if ($MFopt{DoMetaPhlan}){
		if (!-e $dir_MP2."$SmplName.MP2.sto"){
			$calcMetaPhlan=1 ;
			$progStats{metaPhl2FailCnts}++;
		} else { $progStats{metaPhl2ComplCnts}++;}
	}
	my $calcMOTU2=0;
	if ($MFopt{DoMOTU2}){
		if (!-e $dir_mOTU2."$SmplName.Motu2.sto"){
			$calcMOTU2=1;
			$progStats{mOTU2FailCnts}++;
		} else { $progStats{mOTU2ComplCnts}++;}
	}
	my $calcTaxaTar = 0;
	if ($MFopt{DoTaxaTarget}){
		if (!-e $dir_TaxTar."$SmplName.TaxTar.sto"){
			$calcTaxaTar=1 ;
			$progStats{taxTarFailCnts}++;
		} else { $progStats{taxTarComplCnts}++;}
	}


	#TODO
	my $allMapDone =0;#used for SNP calling and Binning - but binning requires info if all maps are finished from all samples
	$allMapDone = 1 if (-e "$finalMapDir/$SmplName-smd.$bamcramMap" && $eCovAsssembly && !$ePreAssmbly && ($eSuppCovAsssembly || !$locMapSup2Assembly) && $AsGrps{$cAssGrp}{MapDeps} !~ m/[^;]/ );
	
	
	my $calcBinning = 0;
	if ($MFopt{DoMetaBat2} && $boolAssemblyOK && $AssemblyGo && $AsGrps{$cAssGrp}{MapDeps} !~ m/[^;]/ &&  (!-e "$BinningOut.cm" && !-s "$BinningOut.cm2") ) {
		$calcBinning=$MFopt{DoMetaBat2};
		#die "$MFopt{DoMetaBat2} && $boolAssemblyOK && $AssemblyGo && $AsGrps{$cAssGrp}{MapDeps} !~ m/[^;]/ &&  (!-e $BinningOut.cm || !-s $BinningOut.cm2\n";
	}
	
	#not complete yet? Then delete..
	if ($MFconfig{redoFails} && ($calcRibofind||$calcDiamond || $calcDiaParse ||$calcMOTU2 || $calcMetaPhlan || $calcTaxaTar)){
		die "now recalc $curSmpl\n";
		system ("rm -r -f $assDir $finalCommAssDir");
		system("rm -f -r $curOutDir $smplTmpDir $collectFinished ");
	}



	#die "$metaGassembly\n$finAssLoc\n$nodeSpTmpD\n";

#	#-----------------------  FLAGS  ------------------------  
#	#and some more flags for subprocesses
	my $nonPareilFlag = !-s "$nonParDir/$SmplName.npo" && $MFopt{DoNonPareil} ;
	my $scaffoldFlag = 0; if ( !-e $STOfinScaff && $map{$curSmpl}{"SupportReads"} =~ m/mate/i ){$scaffoldFlag = 1 ;}# print "SUPP:: $map{$curSmpl}{SupportReads}\n";}
	my $assemblyFlag = 0; $assemblyFlag = 1 if ( $MFopt{DoAssembly} && !$boolAssemblyOK && !$efinAssLoc && !$emetaGassembly);

	#die "$assemblyFlag = 1 if ( $MFopt{DoAssembly} && !$boolAssemblyOK && !$efinAssLoc && !-e $metaGassembly\n $doPreAssmFlag\n";
	my $calcReadMerge = 0;
	$calcReadMerge = 1 if ($MFopt{doReadMerge} && ($MFopt{calcOrthoPlacement} || $calcDiamond || $calcGenoSize));
	my $mapAssFlag = 0; $mapAssFlag = 1 if ($MFopt{map2Assembly} && !$eFinMapCovGZ  );
	my $calcCoverage = 0; $calcCoverage =1 if (!$eCovAsssembly && $MFopt{map2Assembly});
	
	#only for support reads (from hybrid assemblies)
	my $mapSuppAssFlag =0;$mapSuppAssFlag = 1 if ($locMapSup2Assembly && !$eFinSupMapCovGZ && $efinAssLoc  );#hasSuppRds(\%AsGrps,$cAssGrp,$curSmpl ) );
	my $calcSuppCoverage = 0; $calcSuppCoverage =1 if ($MFopt{mapSupport2Assembly} && !$eSuppCovAsssembly && $map{$curSmpl}{"SupportReads"} ne "");
	#die "$mapSuppAssFlag = 1 if ($locMapSup2Assembly && !$eFinSupMapCovGZ && $efinAssLoc \n";

	my $assemblyBuildIndexFlag=0; 
	$assemblyBuildIndexFlag=1 if ($MFopt{DoAssembly}  && !$assemblyFlag  && $MFopt{map2Assembly} && ($mapAssFlag || $mapSuppAssFlag ) 
							# && $AssemblyGo -> do already in first round, store dependency
							&& ! mapperDBbuilt($finAssLoc,$MFopt{MapperProg}) ); 
	#die "$assemblyBuildIndexFlag  $MFopt{DoAssembly}  && !$assemblyFlag  && $MFopt{map2Assembly} && ($mapAssFlag || $mapSuppAssFlag ) \n". mapperDBbuilt($finAssLoc,$MFopt{MapperProg})  ."\n";
	#print "build $assemblyBuildIndexFlag   $MFopt{DoAssembly} && !$assemblyFlag && $MFopt{map2Assembly} && $mapAssFlag && $MFopt{MapperProg}\n";
	#requires only bam/cram && assembly
	my $calcConsSNP=0; $calcConsSNP =1 if ($MFopt{DoConsSNP} && (!-e  $genePredSNP || -s $genePredSNP < 100 || ($MFopt{saveVCF} && ! fileGZe($vcfSNP) )));
	
	my $calcSuppConsSNP=0; $calcSuppConsSNP =1 if ($MFopt{DoSuppConsSNP} && (!-e  $STOsnpSuppCons  ));
	
	
	#die "genePredSNP $genePredSNP\n";
	my $calc2ndMapSNP = 0; $calc2ndMapSNP = 1 if ($MFopt{Do2ndMapSNP});
	my $pseudAssFlag = 0; $pseudAssFlag = 1 if ($MFopt{pseudoAssembly} && $map{$curSmpl}{ExcludeAssem} eq "0" && (!-e $pseudoAssFileFinal.".sto" || !$boolGenePredOK));
	my $dowstreamAnalysisFlag = 0; 
	#$unpackZip simulates dowstreamAnalysisFlag, just to get sdm running..
	$dowstreamAnalysisFlag=1 if ( $MFconfig{unpackZip} || $calcContamination || $MFopt{calcOrthoPlacement} || $scaffTarExternal ne "" || $assemblyFlag  
		|| $pseudAssFlag || $scaffoldFlag  || $nonPareilFlag || $calcGenoSize || $calcDiamond || $calcDiaParse
		|| $MFopt{DoCalcD2s} || $calcKraken || $calcRibofind || $calcRiboAssign || $calcMOTU2 ||  $calcMetaPhlan || $calcTaxaTar );
	my $requireRawReadsFlag = 0;#only list modules that really need raw reads
	$requireRawReadsFlag = 1 if ( !$boolScndMappingOK || $calcContamination);
	
	my $seqCleanFlag = 0; $seqCleanFlag =1 if ($dowstreamAnalysisFlag && !-e "$smplTmpDir/seqClean/filterDone.stone" );
	my $porechopFlag = 0;
	$porechopFlag = 1 if ($MFopt{usePorechop} && $dowstreamAnalysisFlag && !-e "$smplTmpDir/rawRds/poreChopped.stone");
	#die "$assemblyFlag\t$seqCleanFlag\t$boolScndMappingOK\n";
	my $calcUnzip=0;
	$calcUnzip=1 if ($calcDiamond || $porechopFlag || $seqCleanFlag  || $mapAssFlag || $mapSuppAssFlag || (!$MFopt{useUnmapped} && !$boolScndMappingOK)); 
	#print "chk1 $mapSuppAssFlag $calcSuppCoverage $eSuppCovAsssembly\n" ;

	if ($scaffTarExternal ne "" &&  $map{$curSmpl}{"SupportReads"} !~ /mate/i && $scaffTarExtLibTar ne $curSmpl ){print"scNxt\n";loop2C_check($cAssGrp,\@sampleDeps);next;}
	
	if (!$assemblyFlag && $AssemblyGo && !$efinAssLoc && $emetaGassembly){
		die "assem copy\n"; #should actually never be here..
		push(@{$AsGrps{$cAssGrp}{AssCopies}}, $assDir."/metag/*",$finalCommAssDir);
	}


#	#-----------------------  END FLAGS  ------------------------  

	#some more flow control..
	if ( 
		!$DoUploadRawReads && $boolScndMappingOK && !$MFopt{DoCalcD2s} &&
		!$calcConsSNP && !$calcSuppConsSNP && !$calcBinning && !$calc2ndMapSNP && $boolAssemblyOK && $boolScndCoverageOK 
		 && !$calcCoverage && !$calcSuppCoverage && !$dowstreamAnalysisFlag
		#&& !$calcRibofind && !$calcRiboAssign && !$MFopt{calcOrthoPlacement} && !$calcGenoSize && !$calcDiamond && !$calcDiaParse && 
		#!$calcMetaPhlan && !$calcTaxaTar && !$calcMOTU2 && !$calcKraken && $scaffTarExternal eq ""
	){
		#free some scratch
		system "rm -rf $smplTmpDir &" if ($MFconfig{rmScratchTmp} );
		#system "rm -f $smplLockF" if (-e $smplLockF);
		
		if ( ($boolAssemblyOK || ($doPreAssmFlag && $ePreAssmbly && !$ePreAssmblPck)) && !$locRedoAssMapping ){#causes a lot of overhead but mainly to avoid unpacking reads again..
			print "present: $curOutDir \n"; $presentAssemblies ++;#= $AsGrps{$cAssGrp}{CntAimAss};
			#base is present, but is the additions? 

			if ($MFopt{map2Assembly} && $MFopt{doBam2Cram} && #.cram.sto to check that everything went well
				$AssemblyGo && $MappingGo && $eFinMapCovGZ && $allMapDone ){
							#(-e "$finalCommAssDir/scaffolds.fasta.filt.$MFcontstants{mini2IdxFileSuffix}" || -e "$finalCommAssDir/scaffolds.fasta.filt${bwt2IdxFileSuffix}.1.bt2") &&
				#delete bwt2 index to save some space..
				#need to ensure that all mapping in mapping group are fine..
				system "rm -f $finalCommAssDir/scaffolds.fasta.filt$MFcontstants{bwt2IdxFileSuffix}* $finalCommAssDir/scaffolds.fasta.filt.fai $finalCommAssDir/scaffolds.fasta.filt.$MFcontstants{mini2IdxFileSuffix} $finalCommAssDir/scaffolds.fasta.filt.$MFcontstants{kmaIdxFileSuffix}";
				system "rm -f $finalMapDir/${SmplName}*smd.bam.bai $finalMapDir/${SmplName}*smd.bam.coverage.c* $finalMapDir/${SmplName}*smd.bam.coverage.gen* $finalMapDir/${SmplName}*smd.bam.coverage.m* $finalMapDir/${SmplName}*smd.bam.coverage.p*";
				
			}
			if ($MFopt{map2Assembly} && !$MFopt{mapSaveCram} && -s $BinningOut && (-s "$BinningOut.cm2" || -e "$BinningOut.cm") &&
				!$calcBinning && !$calcConsSNP && !$calcSuppConsSNP && $eFinMapCovGZ && $efinAssLoc && $eFinalMapDir && -s $CRAMmap){
#				die;
				print "deleting $CRAMmap to save space..\n";
				system "rm -rf $CRAMmap";
			}
		}
		
		
		print "next due to sample finished";
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps);next;

	}
	

	#report for debugging:
	#print "!$calcConsSNP && !$calcBinning && !$calc2ndMapSNP && $boolAssemblyOK && $boolScndCoverageOK \n	&& $boolScndMappingOK && !$MFopt{DoCalcD2s} && !$DoUploadRawReads \n	&& !$calcRibofind && !$calcRiboAssign && !$MFopt{calcOrthoPlacement} && !$calcGenoSize && !$calcDiamond && !$calcDiaParse && \n	!$calcMetaPhlan && !$calcTaxaTar && !$calcMOTU2 && !$calcKraken && $scaffTarExternal eq \n $allMapDone\n $eFinMapCovGZ && $eCovAsssembly && $calcCoverage\n";
	








#-----------------------------------------------------------------------------------------
#                   actual job submissions starts from here
#-----------------------------------------------------------------------------------------

#unzipping & sorting of files into scratch dir..
#also does either trimomatic or porechop..
#die "calcUnzip $calcUnzip\n" ;#if ($calcUnzip);
#	my ($jdep,$cfp1ar,$cfp2ar,$cfpsar,$WT,$rawFiles, $mmpuNum, $libInfoRef, $inputRawSize) = 
	#my %seqSet = (pa1 => \@pa1, pa2 => \@pa2, pas => \@pas, paX1 => \@paX1, paX2 => \@paX2, paXs => \@paXs, libInfo => \@libInfo, libInfoX => \@libInfoX,
	#		totalInputSizeMB => $totalInputSizeMB,rawReads => $rawReads,mmpu => $mmpu, WT => $WT);
	#unzip and change ifastap & cfp1/cfp2
	my $curUnzipDep = ""; 
	if (0 && $MFconfig{maxUnzpJobs} >0 && @unzipjobs > $MFconfig{maxUnzpJobs}){#only run X jobs in parallel, lest the cluster IO breaks down..
		$curUnzipDep = $unzipjobs[-($MFconfig{maxUnzpJobs})];#join(";",@last_n);
		$waitTime=0;
	}
	my ($jdep,$hrefSeqSet) = 
			seedUnzip2tmp($curDir,$curSmpl,$curUnzipDep,$nodeSpTmpD,
			$smplTmpDir,$waitTime,$AsGrps{$cMapGrp}{CntMap},$calcUnzip,$finalMapDir,
			$porechopFlag,$inputRawFile);
	my %seqSet = %{$hrefSeqSet};
	push (@unzipjobs,$jdep) unless ($jdep eq "");
	$seqSet{samplReadLength} = $samplReadLength; $seqSet{samplReadLengthX} = $samplReadLengthX;
	#die "@{$seqSet{pa1}}\n@{$seqSet{pas}}\n";

	#$seqSet{"curSDMopt"} = $curSDMopt; $seqSet{"curSDMoptSingl"} = $curSDMoptSingl;
	
	if ($jdep eq "EMPTY_DO_NEXT"){push(@EmptySample,$curSmpl);loop2C_check($cAssGrp,\@sampleDeps);next;}
	push(@inputRawFQs,$seqSet{"rawReads"});
	if($scaffTarExtLibTar eq $curSmpl){
		@scaffTarExternalOLib1 = @{$seqSet{"pa1"}}; @scaffTarExternalOLib2 = @{$seqSet{"pa2"}};
		next unless ($map{$curSmpl}{"SupportReads"} =~ /mate/i );
		#print "@scaffTarExternalOLib1\n";
		#die;
	}
	
	
	#$mmpuOutTab .= $dir2rd."\t".$seqSet{"mmpu"}."\n";
	$waitTime = $seqSet{"WT"};
	$AsGrps{$cMapGrp}{SeqUnZDeps} .= $jdep.";";$AsGrps{$cAssGrp}{UnzpDeps} .= $jdep.";";$AsGrps{$cAssGrp}{readDeps} = $AsGrps{$cAssGrp}{UnzpDeps};
	my $UZdep = $jdep;
	#my @cfp1 = @{$seqSet{"pa1"}}; my @cfp2 = @{$seqSet{"pa2"}}; #stores the read files used plus which library they come from
	#if (@cfp1!=@cfp2){print "Fastap path not of equal length:\n@@cfp1\n@@cfp2\n"; die();}
	#push(@{$AsGrps{$cMapGrp}{RawSeq1}},@{$seqSet{"pa1"}}); push(@{$AsGrps{$cMapGrp}{RawSeq2}},@{$seqSet{"pa2"}});
	#push(@{$AsGrps{$cMapGrp}{RawSeqS}},@{$seqSet{"pas"}});push(@{$AsGrps{$cMapGrp}{Libs}},@{$seqSet{"libInfo"}});
	
	
	#empty links for assembler and nonpareil
	#my($arp1,$arp2,$singAr,$matAr,$sdmjN) = ([],[],[],[],"");
	my $sdmjN = ""; #job on main  reads
	my $cleanSeqSetHR = iniCleanSeqSetHR(\%seqSet);	

	#my($arp1,$arp2,$singAr,$matAr,$sdmjN) = ($seqSet{"pa1"},$seqSet{"pa2"},$seqSet{"pas"},[],$jdep);
	#empty links and objects for merging of reads
	my($mergJbN) = ("");
	# punsh the whole thing through sdm.. 
	if ( (!$boolAssemblyOK||$calcContamination) && $MFopt{useSDM}!=0 ){#&& !$is3rdGen) {
		#$ifastp
		#($cleanSeqSet{arp1},$cleanSeqSet{arp2},$cleanSeqSet{singAr},$cleanSeqSet{matAr},$sdmjN) = sdmClean($curOutDir,\%seqSet, 
		($cleanSeqSetHR,$sdmjN)  = sdmClean($curOutDir,\%seqSet, $cleanSeqSetHR, 
				$smplTmpDir."mateCln/",$smplTmpDir."seqClean/",$jdep,$dowstreamAnalysisFlag,0) ;
		# check for support reads as well..
		my $sdmjN2 = ""; #job on support  reads
		($cleanSeqSetHR,$sdmjN2) = sdmClean($curOutDir,\%seqSet, $cleanSeqSetHR,
				$smplTmpDir."mateCln/",$smplTmpDir."seqClean/",$jdep,$dowstreamAnalysisFlag,1) ;
		$sdmjN .= ";$sdmjN2" if ($sdmjN2 ne "");
	}  
	#adds raw and cleaned read file location to the whole assembly group
	my $assGrpHR = addFileLocs2AssmGrp(\%AsGrps, $cAssGrp,$SmplName, $cleanSeqSetHR, \%seqSet);   %AsGrps = %{$assGrpHR};
	
	
	
	#make sure that no mapping is started when postAssmbl not yet done..
	#needs to happen after reads registered (for assembly of multiple samples..
	#$MappingGo=0 
	if ( !$efinAssLoc ){ #there should be no support mapping, unless the hybrid assembly has finished!
		$mapSuppAssFlag=0;$eSuppCovAsssembly=0;
	}
	if ($ePreAssmblPck && !$efinAssLoc && !$postPreAssmblGo){$mapAssFlag =0; $mapSuppAssFlag=0;}
	#print "mapCHK $mapSuppAssFlag $mapAssFlag $ePreAssmblPck && ! $efinAssLoc && ! $postPreAssmblGo\n";
	if ( !$mapAssFlag &&  ( ( $ePreAssmblPck && !$postPreAssmblGo && !$efinAssLoc ) #nothing to do untile $doPreAssmFlag releases
			|| (!$ePreAssmblPck && $ePreAssmbly && !$postPreAssmblGo && !$efinAssLoc )  )
	){ #last sample (assembly) should not map while other maps are still running..
		print "next due to waiting for preassemblies";
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); loop2C_check($cAssGrp,\@sampleDeps); next;
	}


	#filter human or other hosts..
	$sdmjN = removeHostSeqs($cleanSeqSetHR,$nodeSpTmpD,$sdmjN,1) if ($dowstreamAnalysisFlag && (!$boolAssemblyOK || $calcContamination));
	#merge reads?
	($cleanSeqSetHR,$mergJbN) = mergeReads($cleanSeqSetHR,$sdmjN,$smplTmpDir."merge_clean/",$calcReadMerge,$dowstreamAnalysisFlag);
	
	#raw files only required for mapping reads to assemblies, so delete o/w
	#$cfp1ar,$cfp2ar,
	if (!$MFopt{DoAssembly} && $MFconfig{importMocat}==0 && $MFconfig{removeInputAgain} && !$requireRawReadsFlag){ $sdmjN = cleanInput(\%seqSet, $sdmjN,$smplTmpDir);}
	$AsGrps{$cAssGrp}{readDeps} .= ";$mergJbN";

	#keeps track of all sdm jobs
	$sdmjNamesAll .= ";".$sdmjN;
	
	#print "@{${$cleanSeqSetHR}{arp1}} ZZZ \n";
	#die "${$cleanSeqSetHR}{readTec} FF\n";
	
	push(@allFilter1,@{${$cleanSeqSetHR}{arp1}});
	push(@allFilter2,@{${$cleanSeqSetHR}{arp2}});
	#die;
	if ($MFconfig{unpackZip}){
		print "next due to onlyFilterZip 1\n";MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps);next;
	}
	
	
	
	
	
	#-----------------------------------------------------------------
	#---------functions only dependent on UN/FILTERED reads-----------
	#-----------------------------------------------------------------
	#take long reads and filter for very long reads (454 etc might be long enough)
	#then do gene predictions *instead* of gene predictions on assemblies
	#die $pseudAssFlag."\n";
	
	#create new dependency flag for both unzip and sdm dep:
	my $primaryDep = "$UZdep;$sdmjN";
	
	if ($pseudAssFlag ){
		my ($psAssDep, $psFile, $metagDir) = createPsAssLongReads($cleanSeqSetHR, $mergJbN.";".$primaryDep, $pseudoAssFile, $finalCommAssDir, $SmplName);#pseudoAssFileFinal
		if ($psAssDep ne ""){
			$AsGrps{$cAssGrp}{pseudoAssmblDep} = $psAssDep;
		}
		#predict genes on assembly
		my $prodRun = genePredictions($psFile,$metagDir."/genePred/",$psAssDep,$finalCommAssDir,"",$smplTmpDir,1);
		if ($prodRun ne ""){
			$AsGrps{$cAssGrp}{pseudoAssmblDep}  = $prodRun;
		}
		push(@{$AsGrps{$cAssGrp}{PsAssCopies}}, $assDir."/metag/*",$finalCommAssDir);
		$AsGrps{$cAssGrp}{readDeps} .= ";$prodRun";
	}
	
	if ($MFopt{calcOrthoPlacement}){
		runOrthoPlacement($cleanSeqSetHR,$curOutDir."orthos/",$nodeSpTmpD."/orthos/",
					$mergJbN.";".$primaryDep);
	}
	if ($calcDiamond || $calcDiaParse){
		my ($djname,$djCln) = runDiamond($cleanSeqSetHR,$curOutDir."diamond/",$globaldDiaDBdir,$nodeSpTmpD."/diaRefDB/",
					$mergJbN.";".$primaryDep,$MFopt{reqDiaDB}); #GlbTmpPath
		$AsGrps{$cAssGrp}{DiamDeps} = $djname.";";
		$AsGrps{global}{DiamDeps} .= ";$djname";
		$AsGrps{global}{DiamCln} = $djCln unless($djCln eq "");
		$AsGrps{$cAssGrp}{readDeps} .= ";$djname";
	}
	
	#DEBUG
	#my ($arp1,$arp2,$singAr,$mergRdsHr) = ("","","","");
	#non pareil (estimate community size etc)
	if ($nonPareilFlag){
		nopareil(${$cleanSeqSetHR}{arp1},$nonParDir, $globalNPD, $SmplName,$primaryDep);
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps);next;
	}
	if ($calcGenoSize){#use microbeCensus to get avg genome size
		my $gsJdep = genoSize($cleanSeqSetHR,$curOutDir."MicroCens/",$mergJbN.";".$primaryDep);
		$AsGrps{$cAssGrp}{readDeps} .= ";$gsJdep";
	}
	
	#kraken (estimate taxa abundance
	if ($calcKraken){
		my $krJdep = krakenTaxEst($cleanSeqSetHR,$KrakenOD, $nodeSpTmpD."krak/",$SmplName,$primaryDep);
		$AsGrps{$cAssGrp}{readDeps} .= ";$krJdep";
		
	}
	if ($calcRibofind || $calcRiboAssign){
#		die "STOP ribo\n";
		my $ITSrun = detectRibo($cleanSeqSetHR,$nodeSpTmpD."ITS/",$curOutDir."ribos/",$primaryDep,$SmplName,$runTmpDirGlobal); #GlbTmpPath  \@cfp1,\@cfp2
		#$AsGrps{$cAssGrp}{ITSDeps} .= $ITSrun.";";
		$AsGrps{$cAssGrp}{readDeps} .= ";$ITSrun";
	}
	
	#metaphlan2 - taxa abudnance estimates
	if ($calcMetaPhlan){
		my $MP2jname = metphlanMapping($cleanSeqSetHR,$nodeSpTmpD."MP2/",$dir_MP2,$SmplName,$MFopt{MapperCores},$primaryDep); #\@cfp1,\@cfp2
		$AsGrps{$cAssGrp}{readDeps} .= ";$MP2jname";
	}
	
	if ($calcTaxaTar){
		my $taxTarjname = TaxaTarget($cleanSeqSetHR,$nodeSpTmpD."TaxTar/",$dir_TaxTar,$SmplName,$MFopt{MapperCores},$primaryDep); 
		$AsGrps{$cAssGrp}{readDeps} .= ";$taxTarjname";
	}
	
	#mOTU2  - taxa abundance estimates
	if ($calcMOTU2){
		my $MP2jname = mOTU2Mapping($cleanSeqSetHR,$nodeSpTmpD."Motu2/",$dir_mOTU2,$SmplName,$MFopt{MapperCores},$primaryDep.";".$mOTU2Deps); #\@cfp1,\@cfp2
		$AsGrps{$cAssGrp}{readDeps} .= ";$MP2jname";
	}
	
	my $SmplNameX = $SmplName;
	if ($AsGrps{$cAssGrp}{CntAss} > 1){$SmplNameX .= "M".$AsGrps{$cAssGrp}{CntAss};}
	#my @tmp = @{$arp1};die ("ASSflag: ".$assemblyFlag."\n@tmp\n");
	#die "$assemblyFlag\n $AssemblyGo\n";
	
	
	#add dependencies to wait on for loopUntil
	add2SampleDeps(\@sampleDeps, [$AsGrps{$cAssGrp}{readDeps},$AsGrps{$cAssGrp}{DiamDeps},$AsGrps{$cMapGrp}{SeqUnZDeps},$mergJbN,$sdmjN]);
	
	#needs assmebly group centric list of ALL filtered/raw input reads..
	#push(@{$AsGrps{$cAssGrp}{FilterSeq1}},@{${$cleanSeqSetHR}{arp1}});push(@{$AsGrps{$cAssGrp}{FilterSeq2}},@{${$cleanSeqSetHR}{arp2}});push(@{$AsGrps{$cAssGrp}{FilterSeqS}},@{${$cleanSeqSetHR}{singAr}});push(@{$AsGrps{$cAssGrp}{ReadTec}},${$cleanSeqSetHR}{readTec});
	
#	${$AsGrps{$cAssGrp}{CleanSeqs}}{$SmplName} = $cleanSeqSetHR;
#	${$AsGrps{$cAssGrp}{RawSeqs}}{$SmplName} = \%seqSet;
	#DEBUG
	#my ($ar1,$ar2,$ars,$liar,$rear) = getRawSeqsAssmGrp(\%AsGrps, $cAssGrp, 1);
	#my ($ar1,$ar2,$ars,$rear) = getCleanSeqsAssmGrp(\%AsGrps, $cAssGrp, 1);
	#print "@{$ar1}  ..  @{$ars}\n";die;
	#die;






	#-----------------------------------------------------------------
	#------------------------ ASSEMBLY -------------------------------
	#-----------------------------------------------------------------
	$AsGrps{$cAssGrp}{SeqClnDeps} .= $sdmjN.";" if ($assemblyFlag);
	if ( ($assemblyFlag || $scaffoldFlag || $scaffTarExternal ne "") && $AssemblyGo){ #assembly does not exist
		die "Can't do assembly and pseudoassembly on the same sample!\n" if ($pseudAssFlag || $MFopt{pseudoAssembly});
		#print "preAsmChk: $ePreAssmbly, $ePreAssmblPck, $doPreAssmFlag, $postPreAssmblGo\n";
		#die;
		metagAssemblyRun( $cAssGrp,$cleanSeqSetHR,"$nodeSpTmpD/ass",$metagAssDir ,$shortAssembly, $SmplNameX,$scaffoldFlag,$metaGscaffDir,
					$assemblyFlag,$AssemblyGo,$ePreAssmbly, $doPreAssmFlag, $postPreAssmblGo);
		push(@{$AsGrps{$cAssGrp}{AssCopies}}, $assDir."/metag/*",$finalCommAssDir);

		#call genes, depends on assembly
		$AsGrps{$cAssGrp}{prodRun} = genePredictions($metaGassembly,$geneDir,$AsGrps{$cAssGrp}{AssemblJobName},$finalCommAssDir,"",$smplTmpDir,1);
	}
	if ($assemblyBuildIndexFlag && $AsGrps{$cAssGrp}{AssemblJobName} eq ""){ #in this case asembly was done, but index was never built
		buildAssemblyMapIdx($finAssLoc, $cAssGrp, $mapAssFlag,$mapSuppAssFlag,$SmplName);
	}
	
	if (!$assemblyFlag || ($ePreAssmbly && $doPreAssmFlag) ){   # gene predictions on assembly
		$metaGassembly = $finAssLoc; print "No Assembly routines required\n" if ($MFopt{DoAssembly}!=0);
		if (!$boolGenePredOK && $AssemblyGo ){
			$geneDir = $finalCommAssDir."/genePred/";
			$AsGrps{$cAssGrp}{prodRun} = genePredictions($metaGassembly,$geneDir,$AsGrps{$cAssGrp}{AssemblJobName},$finalCommAssDir,"",$smplTmpDir,1);
		}
	}
	if (0&&$AssemblyGo && $AsGrps{$cAssGrp}{PostAssemblCmd} ne ""){#no assembly required, but maybe still other dependent jobs (i.e. mapping)
		print "assembly exists, but postassembly jobs unfinished\n" if (!$assemblyFlag);
		#die "POSTCMD: ".$AsGrps{$cAssGrp}{PostAssemblCmd}."\n";
		postSubmQsub("$logDir/MultiMapper.sh",$AsGrps{$cAssGrp}{PostAssemblCmd},$AsGrps{$cAssGrp}{AssemblJobName},$AsGrps{$cAssGrp}{AssemblJobName});
		$AsGrps{$cAssGrp}{PostAssemblCmd} = "";#always add in dep on read extraction
		#$metaGassembly = $finAssLoc;
		#die "!$boolGenePredOK && $AssemblyGo && !$assemblyFlag\n";
	} 
	
	
	
	#-----------------   secondary mapping -----------------
	#2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd
	# maps on given reference genome(s)
	if (@bwt2outD>0 && $MappingGo && (!$boolScndMappingOK || !$boolScndCoverageOK || $calc2ndMapSNP) ){#map reads to specific tar
		scndMap2Genos($SmplName,$cleanSeqSetHR,$cMapGrp,$cAssGrp,$curOutDir,$nodeSpTmpD,$smplTmpDir,
			\@sampleDeps,$samplReadLength,$calc2ndMapSNP,$boolScndCoverageOK);
	} elsif ($MFopt{mapModeActive}) {
		#still needs delays in cleaning command
		$AsGrps{$cMapGrp}{ClSeqsRm} .= ";".$smplTmpDir;
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps); next;
	}
	
	
	#%%%%%%%%%%%%%%%%   functions dependent on assembly -> submit post-assembly   #%%%%%%%%%%%%%%%%
	#  ---------     primary assembly mapping -----------------------   map reads to Assembly      ------------------------------
	if ($MappingGo && !$eFinMapCovGZ && $MFopt{map2Assembly} && ($MFopt{DoAssembly} || $mapAssFlag)){ #mapping to the assembly (can be multi-sample assembly as well)
		my $moveMappings = 0; $moveMappings =1 if (!$eFinMapCovGZ && -e "$mapOut/$SmplName-smd.bam.coverage.gz");
		my $mapNow =0;
		$mapNow = 1 if (!$moveMappings && ($AssemblyGo ||  $efinAssLoc || ($ePreAssmbly && $doPreAssmFlag) ) );#controls if several samples are assembled together and this needs to be waited for
		#print "main map\n";
		my $unAlDir = "$mapOut/unaligned/";$unAlDir = "" if (!$MFopt{SaveUnalignedReads}); #in most cases we don't need unaligned reads..
		
		my %dirset = 	(nodeTmp=>$nodeSpTmpD,outDir => "$finalMapDir/", unalDir => $unAlDir,
						sbj => $metaGassembly, assGrp => $cAssGrp,  smplName => $SmplName,
						glbTmp => $smplTmpDir."/toMGctgs/",glbMapDir => $mapOut,mapSupport => 0,
						readTec => ${$cleanSeqSetHR}{readTec}, submit => 1,submNow => $mapNow,
						sortCores => $MFopt{bamSortCores}, mapCores => $MFopt{MapperCores}, cramAlig => $cramthebam);
		# primary mapping (onto de novo assembly)
		my ($map2Ctgs,$delaySubmCmd,$mapOptHr) = mapReadsToRef(\%dirset, \%AsGrps,
			$AsGrps{$cAssGrp}{AssemblJobName}.";$jdep");#\@libsCFP);
		my ($map2Ctgs_2,$delaySubmCmd_2,$mapStat)  = bamDepth(\%dirset,$map2Ctgs,$mapOptHr);
	
		
#		die;
		$delaySubmCmd .= "\n".$delaySubmCmd_2;
		$AsGrps{$cAssGrp}{MapDeps} .= $map2Ctgs_2.";";
		$AsGrps{$cAssGrp}{BinDeps} .= $map2Ctgs_2.";";
		my $cpyStrm = "MapCopies";$cpyStrm = "MapCopiesNoDel" if ($mapSuppAssFlag || $eFinSupMapCovGZ);
		$cpyStrm = "nothing" if ($map2Ctgs_2 eq "" );# deactivate copying if no job was submitted..
		if (!${$mapOptHr}{immediateSubm} ){#$map2Ctgs_2 ne ""){
			#store command for later..
			$AsGrps{$cAssGrp}{PostAssemblCmd} .= $delaySubmCmd;
			push(@{$AsGrps{$cAssGrp}{$cpyStrm}},$mapOut."/*",$finalMapDir);
			#print "\n\nmapcopydel\n";
		}elsif ($moveMappings){ 
			#this part is checking only if files were not copied..
			#just do it.. 
			print "Moving mappings from globaltmp to finaldir\n";
			system "mkdir -p $finalMapDir;" unless (-d $finalMapDir);
			systemW "rsync -r   --remove-source-files $mapOut/* $finalMapDir/";
		} elsif ($map2Ctgs_2 ne "") { #make sure files get copied
			push(@{$AsGrps{$cAssGrp}{$cpyStrm}},$mapOut."/*",$finalMapDir);
		}
		#die "$curOutDir/mapping\n" . "$mapOut/$SmplName-smd.bam.coverage.gz"
	}
	if ($mapSuppAssFlag){ #supplementary mappings (eg PacBio in hybrid assemblies)
		print "mapping support reads!!\n";
		my $mapNow = 1;
		my $unAlDir = "$mapOut/unaligned_supp/";$unAlDir = "" if (!$MFopt{SaveUnalignedReads});
		my %dirset = 	(nodeTmp=>$nodeSpTmpD,outDir => "$finalMapDir/", unalDir => $unAlDir,
						sbj => $metaGassembly, assGrp => $cAssGrp,  smplName => $SmplName,
						glbTmp => $smplTmpDir."/toMGctgsSupp/",glbMapDir => $mapOut, mapSupport => 1,
						readTec => "", submit => 1,submNow => $mapNow,
						sortCores => $MFopt{bamSortCores}, mapCores => $MFopt{MapperCores}, cramAlig => $cramthebam);
		# primary mapping of support reads (onto de novo assembly)
		my ($map2Ctgs,$delaySubmCmd,$mapOptHr) = mapReadsToRef(\%dirset, \%AsGrps,
			$AsGrps{$cAssGrp}{AssemblJobName}.";$jdep");#\@libsCFP);
		my ($mapSup2Ctgs_2,$delaySubmCmd_2,$mapStat)  = bamDepth(\%dirset,$map2Ctgs,$mapOptHr);
			$delaySubmCmd .= "\n".$delaySubmCmd_2;
		$AsGrps{$cAssGrp}{MapDeps} .= $mapSup2Ctgs_2.";";
		$AsGrps{$cAssGrp}{BinDeps} .= $mapSup2Ctgs_2.";";
		if ($mapSup2Ctgs_2 eq "" && !$eFinSupMapCovGZ && -e "$mapOut/$SmplName.sup-smd.bam.coverage.gz"){  # just copy over..
			print "Moving supplementary mappings from globaltmp to finaldir\n";
			system "mkdir -p $finalMapDir;" unless (-d $finalMapDir);
			systemW "rsync -r --remove-source-files $mapOut/* $finalMapDir/";
		} else {
			push(@{$AsGrps{$cAssGrp}{MapSupCopies}},$mapOut."/*",$finalMapDir);
		}
#		die;
	}
	
	#moves finished assemblies & mappings, deletes temp dirs, logic for when to do that:
	my $rmRdsFlag=0; $rmRdsFlag=1 if ($MappingGo && ( $AssemblyGo || $efinAssLoc) );
	#die "$rmRdsFlag ($MappingGo && ( $AssemblyGo || $efinAssLoc) )\n";
	my $cln1 = manageFiles($cAssGrp, $cMapGrp, $rmRdsFlag,  $doPreAssmFlag, $curOutDir, $jdep, $smplTmpDir, $AssemblyGo);
	add2SampleDeps(\@sampleDeps, [$jdep , $AsGrps{$cAssGrp}{MapDeps} , $AsGrps{$cAssGrp}{scndMapping},$AsGrps{$cAssGrp}{prodRun} ]);
	
	#die;
	
	
	# calc statsitics concercing readqual, mappings, genes & contigs
	if ($pseudAssFlag || ($AssemblyGo && $MFopt{DoAssembly}) || 
				isLastSampleInAssembly($finalCommAssDir,$curOutDir) ) {
		my $subprts = $MFconfig{defaultContigSubs}."gFG"; $subprts .= "m" if ($MFopt{DoBinning});
		$subprts .= "k" if ($MFopt{kmerAssembly} );$subprts .= "4" if ($MFopt{kmerPerGene});
		my ($contRun,$tmp33,$tmpCDd) = runContigStats($curOutDir ,$cln1.";".$AsGrps{$cAssGrp}{prodRun},$finalCommAssDir,$subprts,1,$samplReadLength,$samplReadLengthX, $nodeSpTmpD,1,6, $curSmpl) unless ($contigStatsUsed);
		#run contig stats
		postSubmQsub("$logDir/MultiContigStats.sh",$AsGrps{$cAssGrp}{PostClnCmd},$AsGrps{$cAssGrp}{CSfinJobName},$contRun);
		$AsGrps{$cAssGrp}{PostClnCmd} = "";$AsGrps{$cAssGrp}{CSfinJobName} = $contRun;
		$jdep = $contRun;
	} elsif ( (exists($AsGrps{$cAssGrp}{MapDeps}) && $AsGrps{$cAssGrp}{MapDeps} =~ m/[^;\s]/ && $MappingGo) || ($eFinMapCovGZ) ) {
		#die "test23  $AsGrps{$cAssGrp}{MapDeps}\n";
		#calculate solely abundance / gene, has to be run after clean & assembly contigstat step and after mapping has started (at all!)
		my ($jn,$delaySubmCmd2,$tmpCDd) = runContigStats($curOutDir ,$cln1 . ";".$AsGrps{$cAssGrp}{CSfinJobName},$finalCommAssDir,$MFconfig{defaultContigSubs},1,$samplReadLength,$samplReadLengthX,$nodeSpTmpD,$AssemblyGo,1, $curSmpl);
		$AsGrps{$cAssGrp}{PostClnCmd} .= $delaySubmCmd2;
		$jdep = $jn;
		$AsGrps{$cAssGrp}{BinDeps} .= ";$jdep" if ($jdep ne "");
	}
	add2SampleDeps(\@sampleDeps, [$cln1,$jdep]);

	#Binning, SNP calling: only after copying files from tmp and running contig stats
	if ( $calcBinning &&$AssemblyGo ){  #$allMapDone rm: this is checked now via $AsGrps{$cAssGrp}{MapDeps}
		my $binnerTmp = $nodeSpTmpD;
		$binnerTmp = $smplTmpDir if ($MFopt{useBinnerScratch});
		my $binDep  = submitGenomeBinner($binnerTmp,$metaGassembly,$BinningOut, $cAssGrp,$smplIDs[-1]);
		add2SampleDeps(\@sampleDeps, [$binDep]);
	}

	
	#die "@{$AsGrps{$cAssGrp}{MapCopies}}\n";
	
	if ( ($calcConsSNP || $calcSuppConsSNP) && $allMapDone){
		#die "conssnp:: $calcConsSNP $allMapDone $finalMapDir\n";
		#my $ofas = "$curOutDir/SNP/genePred/genes.shrtHD.SNPc.fna";
		my %SNPinfo = (gff => "$finalCommAssDir/genePred/genes.gff",
						assembly => $metaGassembly,
						mapD => "$finalMapDir",
						SNPcaller => $MFopt{SNPcallerFlag},
						ofas => $contigsSNP, #primary file of contigs
						genefna => $genePredSNP,genefaa => $genePredAASNP,
						vcfFile => $vcfSNP,vcfFileSupp => "$vcfSNPsupp",gffFile => "$finalCommAssDir/genePred/genes.gff",
						nodeTmpD => $nodeSpTmpD,scratch => "$smplTmpDir/SNP/",
						smpl => $SmplName,bamcram => $bamcramMap,minDepth => $MFopt{consSNPminDepth},
						depthF => $coveragePerCtg,firstInSample => 1, #($i == 0 ? 1 : 0)
						bpSplit => 1e6,runLocal => 1,SeqTech => $map{$curSmpl}{SeqTech},SeqTechSuppl => "",
						cmdFileTag => "ConsAssem",maxCores => $MFopt{maxSNPcores},memReq => $MFopt{memSNPcall},
						dependency => $AsGrps{$cAssGrp}{BinDeps},split_jobs => $MFopt{SNPconsJobsPsmpl},
						overwrite => $MFopt{redoSNPcons},
						STOconSNP => $STOsnpCons, STOconSNPsupp => "",
						minCallQual => $MFopt{SNPminCallQual},
					);
		$SNPinfo{SeqTechSuppl} = "PB" if (  $map{$curSmpl}{"SupportReads"} =~ m/PB:/);
		if ($calcSuppConsSNP){
			$SNPinfo{STOconSNPsupp} = $STOsnpSuppCons   ; #trigger for also looking at cons SNP for support reads
		}
		my $consSNPdep = createConsSNP(\%SNPinfo);
		add2SampleDeps(\@sampleDeps, [$consSNPdep]);
		#push(@sampleDeps, $consSNPdep) if (defined $consSNPdep && $consSNPdep ne "");
	}

	MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
	
	### loop2complete functionality
	loop2C_check($cAssGrp,\@sampleDeps);

	#print "END2\n@sampleDeps\n".@sampleDeps."\n";
}






print "\n\n###################################\n".$baseOut."\nFINISHED MG-TK submission loop\n";

postprocess();


print "###################################\n\n";






close $QSBoptHR->{LOG};
exit(0);







































#--------------------------------------------------------------
#get ref genomes
#print("/g/bork5/hildebra/bin/cdbfasta/cdbfasta $PaulRefGenomes");
#if (!-f $thisRefSeq){system("/g/bork5/hildebra/bin/cdbfasta/cdbyank -a $GID $PaulRefGenomes.cidx > $thisRefSeq");}

#--------------------------------------------------------------
#and the assembly
#$cmd = "/g/bork5/hildebra/dev/Perl/assemblies/./multAss.pl $DBpath2 $assDir2 $assDir2/tmp/ $Ref 4 $CutSeq";

#####################################################
#
#####################################################

sub loop2C_check(){
	my ($cAssGrp,$sampleDeps_AR) = @_;
	if ($loop2completion ){
		push (@grandDeps, @{$sampleDeps_AR});
		if ($JNUM == ($to-1)){
			#print "L2C:: $loop2completion  @{$sampleDeps_AR}\n";
			qsubSystemJobAlive( \@grandDeps,$QSBoptHR ,1 );
			#reset some key params...
			@grandDeps = ();
			resetAsGrps(\%AsGrps);
			
			#$totalChecked -= ($JNUM - $from);
			$JNUM=($from-1);$loop2completion--;
			print "\n\n-------------------------------------------\n-------------------------------------------\n";
			print "Repeating samples loop: going into iteration  - $loop2completion: \n";
			print "Reanalyzing samples $from till $to\n";
			if ($loop2c_winsize > 0 && !$loop2completion){
				my $tmpStr = "Changing sample window from $from -> $to to ";
				$from = $to;$to = $from + $loop2c_winsize; 
				if ($to > $TO1){
					$to = $TO1 ;
					if ($from-1 >= $TO1){
						print "Last loop, breaking..\n";$from=$TO1;
					}
				} else {$loop2completion = $loop2completion_ini;}
				print "$tmpStr$from -> $to \n";#" . ($JNUM + $loop2c_winsize) . "\n";
			}
			$statStr = ""; 
			print "-------------------------------------------\n-------------------------------------------\n";
		}
	}
}
sub postprocess{

	#print "\n\n###################################\nMain Loop done\n######################################\n";
	#print "$sharedTmpDirP\n";
	#global clean up cmds (like DB removals from scratch)
	print "Postprocessing:\n";

	#print input files, sorted by samples
	if (@inputRawFQs == @allSmplNames){
		open O,">$baseOut/Input_raw.txt";
		for (my $i=0;$i<@allSmplNames; $i++){
			print O "$allSmplNames[$i]\t$inputRawFQs[$i]\n";
		}
		close O;
	}

	if ($statStr ne ""){
		open O,">$baseOut/metagStats.txt";
		print O $statStr;
		close O;
		print "Stats in $baseOut/metagStats.txt\n";
	}
	if (0 && $statStr5 ne ""){
		open O,">$baseOut/metagStats_500.txt";
		print O $statStr5;
		close O;
		#print "Stats in $baseOut/metagStats_500.txt\n";
	}



	#my $dir_MP2 = $baseOut."pseudoGC/Phylo/MP2/"; #metaphlan 2 dir
	#my $dir_RibFind = $baseOut."pseudoGC/Phylo/RiboFind/"; #metaphlan 2 dir
	#my $dir_KrakFind = $baseOut."pseudoGC/Phylo/KrakenTax/"; #metaphlan 2 dir

	#merging of per sample assignments..
	riboSummary();
	mergeMP2Table($dir_MP2);
	mergeMotu2Table($dir_mOTU2);
	unploadRawFilePostprocess();
	DiaPostProcess("",$baseOut);
	d2metaDist(\@allSmplNames,\@allFilter1,\@allFilter2,$sdmjNamesAll,$baseOut."/d2StarComp/");

	if ($MFopt{DoKraken} ){
		if( $progStats{KrakTaxFailCnts}){
			print "$progStats{KrakTaxFailCnts} samples with incomplete KrakenTax\n";} 
		else {print "All samples have assigned Kraken Taxonomy\n";
			my $mergeTblScript = getProgPaths("metPhl2Merge");#"/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
			opendir D, $dir_KrakFind; my @krkF = grep {/0.*/ && -d $dir_KrakFind.$_} readdir(D); closedir D;
			my $mrgCmd = "";
			foreach my $kf (@krkF){
				#$kf =~ m/krak\.(.*)\.cnt\.tax/; my $thr = $1;# die $thr."  $kf\n";
				$mrgCmd .= "$mergeTblScript $dir_KrakFind/$kf/*.krak.txt > $dir_KrakFind/Krak.$kf.mat\n";
			}
			#die $mrgCmd."\n$dir_KrakFind\n";
			system "$mrgCmd";
		}
	}

	#open O,">$baseOut/MMPU.txt";print O $mmpuOutTab;close O;

	#final message with instructions on how to build a gene cat..
	if ($MFopt{DoAssembly}){
		my $warnMsg = "";
		$warnMsg = "(may be inaccurate due to loop2complete)" if ($loop2completion_ini);
		print "Found ".$presentAssemblies." of ".$totalChecked." samples already assembled $warnMsg\n";
		
	}
	my $totalScratchUse=0;
	foreach my $smpl (keys %inputFileSizeMB){
		#input size for 1) raw 2) filtered 3) cram + assembl + consensus
		$totalScratchUse += $inputFileSizeMB{$smpl} * 3 + ($inputFileSizeMB{$smpl}/15*2);
	}
	
	
	print "Estimated scratch use: " . int($totalScratchUse/1024)."G\n";
		if (@EmptySample>0){
		print "Found Empty samples:\n".join(",",@EmptySample) ."\n\n";
	}


	if ($presentAssemblies >0 &&$presentAssemblies == $totalChecked){
		my $gcScr = getProgPaths("geneCat_scr");
		my $GCsub = $baseOut."/GeneCat_pre.sh";
		my $gcmd = "";
		$gcmd .= "#creates gene catalog in the specified outdir with specified cores, attempting to reuse existing dirs (in case catalog creation failed):\n";
		my $sugGCmem = 100; 
		$sugGCmem = 200 if ($presentAssemblies > 100);
		$sugGCmem = 700 if ($presentAssemblies > 500);
		$sugGCmem = 1200 if ($presentAssemblies > 1000);
		$sugGCmem = 2500 if ($presentAssemblies > 5000);
		my $sugGCcores = 12;
		$sugGCcores = 24 if ($presentAssemblies > 100);
		$sugGCcores = 32 if ($presentAssemblies > 500);
		$sugGCcores = 48 if ($presentAssemblies > 1000);
		$sugGCcores = 72 if ($presentAssemblies > 5000);

		$gcmd .= "$gcScr -map $MFconfig{mapFile} -GCd [insert outdir] -mem $sugGCmem -cores $sugGCcores -clusterID 95 -doStrains $MFopt{DoConsSNP} -continue 1 -binSpeciesMG $MFopt{DoMetaBat2} -useCheckM2 $MFopt{useCheckM2} -useCheckM1 $MFopt{useCheckM1} -MGset GTDB \n";
		print "\n\nNext step, create a genecatalog with call to (but modify .sh first!): \nsbatch $GCsub\n";
		$QSBoptHR->{doSubmit} = 0;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		my ($jobN, $tmpCmd) = qsubSystem($GCsub,$gcmd,1,"80G","GeCat","","",1,[],$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$QSBoptHR->{doSubmit} = 1;
	}

}





sub spaceInAssGrp{
#determine how much space is used in total assembly group..
	my ($curSmplX) = @_;
	my $inputSizeloc = 0;
	#print "map: $map{$curSmplX}";
	my @curMems = (); 
	@curMems = @{$map{$curSmplX}{AG_members}} if (defined $map{$curSmplX}{AG_members});
	my $missedSmpls=0;my $totalSmpls = scalar(@curMems);
	foreach my $memsSmpls (@curMems){
		#print "$memsSmpls $inputFileSizeMB{$memsSmpls}\n";
		if (exists($inputFileSizeMB{$memsSmpls})){
			$inputSizeloc += $inputFileSizeMB{$memsSmpls} ;
		} else {
			#print "Warning: file size estimation not registered for $memsSmpls\n";
			$missedSmpls++;
		}
	}
	$inputSizeloc += $inputFileSizeMB{$curSmplX} if (!@curMems);
	#estimate to account for missed samples..
	if ($missedSmpls){$inputSizeloc *= (($totalSmpls-$missedSmpls)/$totalSmpls);} 
	return $inputSizeloc;
}

sub submitGenomeBinner{
	my ($nodeSpTmpD,$metaGassembly,$MetaBat2out,$cAssGrp,$smplIDs1) = @_; #$finalCommAssDir,
	#$MFopt{DoMetaBat2} = 1: metabat2, 2: SemiBin 3: MetaDecoder
	
	#die;
	#if (!-s "$MetaBat2out.cm2"){system "rm $MetaBat2out.cm2";}
	my $BinDir = $MetaBat2out;$BinDir =~ s/[^\/]+$//; #\/[^\/]+\/
	#my $finalCommAssDir = $MetaBat2out;$finalCommAssDir =~ s/\/Binning[^\/]+$//;
	my $MB2coresL = $MFopt{BinnerCores};
	my $HDDspL = $HDDspace{metabat2};
	my $inputSizeloc = spaceInAssGrp($curSmpl);

	my $nodeSpTmpD2= $nodeSpTmpD."/Bin$JNUM$cAssGrp$MFopt{DoMetaBat2}/";
	
	my $mb2Qual = getProgPaths("mb2qualCheck_scr");
	#if ($MFconfig{rmBinFailAssmbly}){print"Warning submitGenomeBinner::\n\nremoving $finalCommAssDir\n\n";system"rm -r $finalCommAssDir";return;}
	my @paths = @{$DOs{$cAssGrp}{wrdir}};#split /,/,$allPaths;
	my $smplIncl = scalar(@paths);
	
	my $CM1done = 0; my $CM2done = 0; my $eBinAssStat=0;
	$CM1done = 1 if (-e "$MetaBat2out.cm" );
	$eBinAssStat =1 if (-e "$MetaBat2out.assStat");
	my $eStone = 0;$eStone = 1 if (-e "$BinDir/Binning.stone");
	$CM2done = 1 if (-s "$MetaBat2out.cm2" );
	
	#die "($CM2done && $MFopt{useCheckM2})";
	if ($eBinAssStat && (( $MFopt{useCheckM1} && $CM1done) || ($CM2done && $MFopt{useCheckM2}) ) ){
		return;
	}
	#die "$MetaBat2out\n";
	my $totMem = 120;#80;
	$totMem = 90 if ($CM1done || !$MFopt{useCheckM1});
	#die "$HDDspL $totMem\n";
	if ($smplIncl >20 || $inputSizeloc > 1e5 ){$HDDspL = greaterComputeSpace($HDDspL,480);$totMem = greaterComputeSpace($totMem,400);
	}elsif ($smplIncl >10 || $inputSizeloc > 5e4 ){$HDDspL = greaterComputeSpace($HDDspL,280);$totMem = greaterComputeSpace($totMem,260);
	}elsif  ($smplIncl >5 || $inputSizeloc > 2e4 ){$HDDspL = greaterComputeSpace($HDDspL,180);$totMem = greaterComputeSpace($totMem,200);
	}elsif  ($smplIncl >3 || $inputSizeloc > 2e4 ){$HDDspL = greaterComputeSpace($HDDspL,120);$totMem = greaterComputeSpace($totMem,180);
	}elsif  ($smplIncl >1 || $inputSizeloc > 1e4 ){$HDDspL = greaterComputeSpace($HDDspL,100);$totMem = greaterComputeSpace($totMem,120);
	}
	$totMem = greaterComputeSpace($totMem,220) if ($MFopt{useCheckM1});
	$totMem = greaterComputeSpace($totMem,$MFopt{BinnerMem}) if ($MFopt{BinnerMem} > 0);
	
	
	#problems with NRP cluster: some jobs require more mem and no nodes available for these..
	if ($MFopt{useBinnerScratch} ){$HDDspL=0;}
	
	#die "$totMem\n";
	my $MBcmd = "mkdir -p $nodeSpTmpD2\n";
	my $BinnerName = "none";
	
	$MBcmd .= "\n\n#checking that all required mappings are done\n". checkMapsDoneSH(\@paths) ."#checks done\n\n";
	#die "$MBcmd\n\n binner\n";
	
	my $seqTec = "hiSeq";
	if (exists($map{$curSmpl}{SeqTech})){$seqTec = $map{$curSmpl}{SeqTech};}
	if ($MFopt{DoAssembly} == 5){$seqTec = "hybrid";}
	#die "$seqTec\n";
	

	#execute calcs later in perl script..
	my $BinnerScr = getProgPaths("Binner_scr");
	$MBcmd .= "$BinnerScr -binner $MFopt{DoMetaBat2} -binD $BinDir -smplID \"$smplIDs1\" -tmpD \"$nodeSpTmpD2\" -assmbl $metaGassembly -assmblGrp $cAssGrp -cores $MB2coresL -smplDirs " . join(",",@paths) . " -seqTec \"$seqTec\" \n";
	
	
	if ($MFopt{DoMetaBat2} == 1){
		$BinnerName = "MB2";
	} elsif ($MFopt{DoMetaBat2} == 2){#SemiBin
		$BinnerName = "SB";
		$QSBoptHR->{useGPUQueue} = 1;
	}elsif ($MFopt{DoMetaBat2} == 3){#MetaDecoder
		$BinnerName = "MD";
	}

	system "rm $MetaBat2out*" if (-e $MetaBat2out && -s $MetaBat2out == 0 && !$CM1done && !$CM2done); #hard flag to just redo calculations..
	$MBcmd = "" if (-s $MetaBat2out);
	
	#die $MBcmd;
	
	#my $MetaBat2out = "$finalCommAssDir/Binning/MB2/$smplIDs[-1]";
	my $postCmd = "";
	my $numCoreCHKM = 3; #fixed --pplacer_threads in routine already..
	if ( ( $MFopt{useCheckM1} && !$CM1done) || (!$CM2done && $MFopt{useCheckM2})  || !$eBinAssStat){
		$postCmd = "\n\nrm -rf $nodeSpTmpD2; mkdir -p $nodeSpTmpD2;\n";
		$postCmd .= "$mb2Qual $metaGassembly $MetaBat2out $nodeSpTmpD2 $MB2coresL $MFopt{useCheckM2} $MFopt{useCheckM1} $MFopt{DoMetaBat2}\n" ;
	} 
	$postCmd .= "\nrm -rf $nodeSpTmpD2\n";
	#die "$MBcmd$postCmd";
	my $jobName = "Bin$JNUM";
	my $preHDDspace=$QSBoptHR->{tmpSpace};$QSBoptHR->{tmpSpace} = $HDDspL;
	#die "$AsGrps{$cAssGrp}{BinDeps}  , $cAssGrp\n";
	my $deps = ""; $deps = $AsGrps{$cAssGrp}{BinDeps} if (defined($AsGrps{$cAssGrp}{BinDeps}));
	my ($jobName2, $tmpCmd) = qsubSystem($paths[-1]."LOGandSUB/Binner$BinnerName.sh", $MBcmd.$postCmd,
			$MB2coresL, int($totMem/$MB2coresL)."G" , $jobName, $deps,"",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} = $preHDDspace;$QSBoptHR->{useGPUQueue} = 0;
#	die $jobName2;
	return $jobName2;
}

sub createConsSNP{
	my ($SNPinfohr) = @_;
	my %SNPinfo = %{$SNPinfohr};
	my $preHDDspace=${$QSBoptHR}{tmpSpace};
	$SNPinfo{qsubDir} = "$logDir/SNP/" unless (exists($SNPinfo{qsubDir}));
	$SNPinfo{JNUM} = $JNUM;
	my $ASFS = filsizeMB($SNPinfo{assembly});
#	${$QSBoptHR}{tmpSpace}  = int($inputFileSizeMB{$curSmpl}*15/1024)+15  ."G"; #*20 for SNPconsensus_vcf #= $HDDspace{SNPcall};
	${$QSBoptHR}{tmpSpace}  = int($ASFS*400/1024)+15  ."G"; #*20 for SNPconsensus_vcf #= $HDDspace{SNPcall};

	#die "${$QSBoptHR}{tmpSpace}\n";
	$SNPinfo{QSHR} = $QSBoptHR;

	# $SNPinfo{jdeps} not used
	#system "rm -f $SNPinfo{mapD}/multi.vcf*" if ($SNPinfo{overwrite});
	my $overwrite = 0;
	unless (exists($SNPinfo{MAR})){
		my @mapping = ($SNPinfo{mapD}."/".$SNPinfo{smpl}."-smd.".$SNPinfo{bamcram});
		$SNPinfo{MAR} = \@mapping;
	}
	if ($SNPinfo{STOconSNPsupp} ne "" && !exists($SNPinfo{MARsupp})){
		my @mapping = ($SNPinfo{mapD}."/".$SNPinfo{smpl}.".sup-smd.".$SNPinfo{bamcram});
		$SNPinfo{MARsupp} = \@mapping;
		#$finalMapDir/$SmplName.sup-smd.bam.coverage.gz
		
	} 
#	my ($ovcf,$jdep) = SNPconsensus_vcf2(\%SNPinfo);
	my ($ovcf,$jdep) = SNPconsensus_vcf(\%SNPinfo);

	${$QSBoptHR}{tmpSpace} = $preHDDspace;
	#SNPconsensus_fasta($ovcf,\%SNPinfo,$jdep,$QSBoptHR);
	
	return $jdep;
}


sub DiaPostProcess(){
	my ($cmd,$baseOut) = @_;
	#also do higher level summaries of diamond to DB mappings
	
	return if (0 || !$MFopt{DoDiamond} ); #temporary deactivated
	
	#die;
	my $mrgDiScr = getProgPaths("mrgDia_scr");
	my @DBS = split/,/,$MFopt{reqDiaDB};
	foreach my $DB (@DBS){
		$progStats{$DB}{SearchCompl} =0 unless (exists($progStats{$DB}{SearchCompl}));
		$progStats{$DB}{SearchIncomplete}=0 unless (exists($progStats{$DB}{SearchIncomplete}));
		my $countFile = "$baseOut/pseudoGC/FUNCT/$DB.compl";
		my $refDone = 0; $refDone = int(getFileStr($countFile)) if (-e $countFile);
		next if (!exists($progStats{$DB}{SearchCompl}) || $progStats{$DB}{SearchCompl} <= $refDone );
		print "$DB :: $progStats{$DB}{SearchCompl} ($refDone previously done)\n";
		$cmd .= "$mrgDiScr $baseOut $DB\necho $progStats{$DB}{SearchCompl} > $countFile\n" if ($progStats{$DB}{SearchCompl}>= 1 );
		#`echo $progStats{$DB}{SearchCompl} > $countFile`;
		#print "$DB: complete $progStats{$DB}{SearchCompl} >= incomplete $progStats{$DB}{SearchIncomplete}\n"
	}
	$cmd .= $AsGrps{global}{DiamCln};
	print "\nMerg diamond:$cmd\n\n";
	systemW $cmd;

}
sub mergeMotu2Table($){
	my ($inD) = @_;
	return unless ($MFopt{DoMOTU2});
	if ($progStats{mOTU2FailCnts}){print "$progStats{mOTU2FailCnts} / ". ($progStats{mOTU2ComplCnts}+$progStats{mOTU2FailCnts}) ." samples with incomplete mOTU2 assignments\n"; return;}
	print "\nAll samples ($progStats{mOTU2ComplCnts}) have mOTU2 assignments.\n";
	my $outD = $inD;
	$outD =~ s/[^\/]+\/?$//;

	my $m2mrgSto = "$outD/m2.Smpl.cnts.stone";
	my $oldM2Cnt = 0;
	if (-e $m2mrgSto){$oldM2Cnt = `cat $m2mrgSto`; chomp $oldM2Cnt;}
	#print "$m2mrgSto $oldM2Cnt $progStats{mOTU2ComplCnts}\n";
	if ($oldM2Cnt >= $progStats{mOTU2ComplCnts} && -e "$outD/m2.class.txt"){return;}

	
	#if (-e "$outD/m2.$taxLvlN[0].txt"){return;}
	
	print "Merging motu2 files into tax tables..\n";
	my $mrgMOTU2 = getProgPaths("mrgMOUT2_scr");
	my $mrgCmd = "$mrgMOTU2 $inD $progStats{mOTU2ComplCnts}\n";
	#die "$mrgCmd\n";
	my ($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/MOTU2merg.sh",$mrgCmd,1,"80G","MOTU2mrg","","",1,[],$QSBoptHR) ;
}



#calculates the hmm based freq estimates and divergence from these -> used for dist matrix
sub d2metaDist{
	my ($arSmpls,$arPaths1,$arPaths2,$deps,$outPath) = @_;

	
	my @paths = @{$arPaths1}; my @Smpls = @{$arSmpls};
	if(!$MFopt{DoCalcD2s}){return;}
	if (@paths < 1){print "Not enough samples for d2s!\n";return;}
	my $d2metaBin = getProgPaths("d2meta");#"/g/bork3/home/hildebra/bin/d2Meta/d2Meta/d2Meta.out";
	print "Calculating kmer distances for ".@paths." samples\n";
	system "mkdir -p $outPath/LOGandSUB";
	my $d2MK = 6;
	my $smplFile = "$outPath/mapd2s.txt";
	open O,">$smplFile";
	for (my $i=0;$i<@paths;$i++){
		print O "$paths[$i] $Smpls[$i]\n";
	}
	close O;
	my $cmd = "";
	$cmd .= "cd $outPath\n";
	$cmd .= "$d2metaBin $d2MK $smplFile Q\n"; #or Q for fastQ
	$cmd .= "touch $outPath/d2meta.stone\n";
	my $jobN = ""; my $tmpCmd;
	unless (-e "$outPath/d2meta.stone"){
		$jobN = "_d2met";
		($jobN, $tmpCmd) = qsubSystem($outPath."/LOGandSUB/d2Met.sh",$cmd,1,"80G",$jobN,"$deps","",1, $QSBoptHR->{General_Hosts},$QSBoptHR) ;
	}
	return $jobN;
}


sub postSubmQsub(){#("$logDir/MultiMapper.sh",$AsGrps{$cAssGrp}{PostAssemblCmd},$AsGrps{$cAssGrp}{AssemblJobName},$SmplName);
	my ($outf,$cmdM,$placeh,$jobN) = @_;
	#die;
	#return; #deactivated since doesn't work on slurm
	return if ($cmdM eq "" || !$doSubmit);
	
	
	#print "$cmdM\n";
	if (0){
		print "Replacing Job $placeh with $jobN\n";
		$cmdM =~ s/$placeh/$jobN/g;
		$cmdM =~ s/-w "done\(\)"//g;
		$cmdM =~ s/--dependency=afterok[^:]+: / /g; #slurm
		#die "$cmdM\n";
	}
	my @spl = split /\n;/,$cmdM;
	foreach my $sp (@spl){
		if ($sp =~ /sbatch.*([\S].sh)/){
			if (-e $1){
				my $fil = `cat $1`; my @spl2 = split /\n/,$fil;
				open O,">$1" or die "Can't open $1\n";
				foreach my $li (@spl2){
					if ($li =~ m/#SBATCH --dependency=afterok:/){
						$li =~ s/$placeh/$jobN/g;
					}
					print O $li;
				}
			}
		}
	}

	#die $AsGrps{$cAssGrp}{PostAssemblCmd}."\n";
	open O,">$outf"; print O $cmdM; close O;
	#print "$outf\n";
	#die;
	system "bash $outf";
}


sub RiboMeta($ $ $ $){
	my ($calcRibofind,$calcRiboAssign,$curOutDir,$SmplName) = @_;
	if ($calcRibofind || $calcRiboAssign){
		if ($MFopt{RedoRiboThatFailed} ){
			system "rm -r $curOutDir/ribos/";
			$calcRibofind = 1;
		}
		$progStats{riboFindFailCnts} ++ ;
	} elsif ($MFopt{DoRibofind} && !$calcRiboAssign) { #all done, copy files to central dir for postprocessing..
		my @RFtags = ("SSU","LSU");#"ITS",
		foreach my $RFtag (@RFtags){
			system "mkdir -p $dir_RibFind/$RFtag/" unless (-d "$dir_RibFind/$RFtag/"); #system "mkdir -p $dir_RibFind/SSU/" unless (-d "$dir_RibFind/SSU/"); system "mkdir -p $dir_RibFind/LSU/" unless (-d "$dir_RibFind/LSU/");
			my $fromCp = "$curOutDir/ribos/ltsLCA/${RFtag}riboRun_bl.hiera.txt"; my $toCpy = "$dir_RibFind/$RFtag/$SmplName.$RFtag.hiera.txt";
			if ($MFopt{checkRiboNonEmpty}){
				#pretty hard check
				my $numLines=0;
				if (-e "$fromCp.gz"){$numLines = `zcat $fromCp.gz | wc -l`;
				} else {$numLines = `wc -l $fromCp`;} $numLines =~ /(\d+)/; $numLines=$1;
				#die $numLines."\n";
				if ($numLines<=1){$calcRiboAssign=1;$calcRibofind=1;
					system "rm -r $curOutDir/ribos//ltsLCA $curOutDir/ribos/*.sto ";last;
				}
			}
			if (!-e $toCpy  || ( (-e $toCpy  || -e "$toCpy.gz" ) && -e $fromCp && -s $fromCp != -s $toCpy)){
				unlink "$toCpy" if -e ($toCpy);
				if (-e "$fromCp.gz"  ){
					#system "zcat $fromCp.gz > $toCpy" ;
					system "rm -f $toCpy.gz;ln -s $fromCp.gz $toCpy.gz" if (!-e "$toCpy.gz");
				} elsif (-e $fromCp) {
					system "gzip $fromCp";
					system "rm -f $toCpy.gz;ln -s $fromCp.gz $toCpy.gz" if (!-e "$toCpy.gz");
				} elsif($RFtag eq "SSU") {#just redo.. SSU is only essential thing
					system "rm -rf $curOutDir/ribos\n"; 
					$calcRibofind = 1; $calcRiboAssign=1;
					last;
				}
			}
			#system "gzip $fromCp" unless (-e "$fromCp.gz");
			system "rm -f $curOutDir/ribos/ltsLCA/inter${RFtag}riboRun_bl.fna" if (-e "$curOutDir/ribos/ltsLCA/inter${RFtag}riboRun_bl.fna");
		}
		$progStats{riboFindComplCnts} ++; #completed already
	} 
}
sub riboSummary{
	return if (!$MFopt{DoRibofind} );
	if( $progStats{riboFindFailCnts}){
		print "$progStats{riboFindFailCnts} / $progStats{riboFindComplCnts} samples with incomplete RiboFind\n";
		return;
	} 
	my $prevItems = 0;
	if ( -e "$dir_RibFind/SSU.cnt.stone"){
		my $tmp = `cat $dir_RibFind/SSU.cnt.stone`;
		chomp $tmp;  $prevItems = $tmp+0;
	}
	my $mergeMiTagScript = getProgPaths("mrgMiTag_scr");#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/miTagTaxTable.pl";
	print "All samples ($progStats{riboFindComplCnts}) have assigned RiboFinds\n";
	my $ITSpres=0;
	my @lvls = ("domain","phylum","class","order","family","genus","species","hit2db");
	#unless (!-d "$dir_RibFind/ITS/"){$mrgCmd .= "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/ITS.miTag $dir_RibFind/ITS/ \n"; $ITSpres=1;}
	my $mrgCmd = "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/SSU.miTag $dir_RibFind/SSU/ \n" unless (!-d "$dir_RibFind/SSU/");
	$mrgCmd .= "echo \"$progStats{riboFindComplCnts}\" > $dir_RibFind/SSU.cnt.stone\n";
	my $mrgCmd2 = "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/LSU.miTag $dir_RibFind/LSU/ \n" unless (!-d "$dir_RibFind/LSU/");
	$mrgCmd2 .= "echo \"$progStats{riboFindComplCnts}\" > $dir_RibFind/LSU.cnt.stone\n";
	#die $mrgCmd."\n";
	my $of_exist = 1;
	foreach my $lvl (@lvls){ 
		$of_exist=0 unless ((!$ITSpres || -e "$dir_RibFind/ITS.miTag.$lvl.txt") && -e "$dir_RibFind/LSU.miTag.$lvl.txt" && -e "$dir_RibFind/SSU.miTag.$lvl.txt"); 
		#die "$dir_RibFind/ITS.miTag.$lvl.txt\n$dir_RibFind/LSU.miTag.$lvl.txt\n$dir_RibFind/SSU.miTag.$lvl.txt\n" if (!$of_exist);
	}
	if ($prevItems < $progStats{riboFindComplCnts}){
		print "Redoing ribo tables, as more samples currently available ($progStats{riboFindComplCnts}, prev:$prevItems)\n";
		$of_exist = 0;
	}
	#die;
	#system $mrgCmd."\n" ; die;
	if ($of_exist == 0){
		my ($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/SSUmerge.sh",$mrgCmd,1,"80G","SSUmrg","","",1,[],$QSBoptHR) ;
		($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/LSUmerge.sh",$mrgCmd2,1,"80G","LSUmrg","","",1,[],$QSBoptHR) ;
	}
	#system "$mrgCmd";
}

sub detectRibo(){
	my ($cleanSeqSetHR , $tmpP,$outP,$jobd,$SMPN,$glbTmpDDB) = @_;
	my $ar1 = ${$cleanSeqSetHR}{arp1}; my $ar2 = ${$cleanSeqSetHR}{arp2}; my $sa1 = ${$cleanSeqSetHR}{singAr}; 
	#print "DEP: $jobd\n";
	my $numCore = 12;
	my $numCore2 = 12;
	my $cLSUSSUscript = getProgPaths("cLSUSSU_scr");#"perl /g/bork3/home/hildebra/dev/Perl/16Stools/catchLSUSSU.pl";
	#my $lambdaIdxBin = getProgPaths("lambdaIdx");
	my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";
	
	my $srtMRNA_path = getProgPaths("srtMRNA_path");
	
	my @re1 = @{$ar1}; my @re2 = @{$ar2}; my @singl = @{$sa1};
	#print "ri"; 
	if (@re1 > 1 || @re2 > 1){
		#print"\nWARNING::\nOnly the first read file will be searched for ribosomes\n";
		#first create new tmp file
	}
	
	#copy DB to server
	my $DBrna = "$glbTmpDDB/rnaDB/";
	my $DBrna2 = "$glbTmpDDB/LCADB/";
	if ($MFopt{globalRiboDependence}->{DBcp} eq "" ){

		my $DBcmd = "";
		my $DBcores = 1;
		$MFopt{globalRiboDependence}->{DBcp}="alreadyCopied";
		#my $ITSfilePref = $1;
		my @DBs = split(/,/,getProgPaths("srtMRNA_DBs"));
		my @DBsIdx = @DBs;my @DBsTestIdx = @DBs; my $filesCopied = 1;
		#die @DBs."@DBs\n";
		if ( !-d $DBrna ){
			$DBcmd .= "mkdir -p $DBrna\n";
		}
		for (my $ii=0;$ii<@DBsIdx;$ii++){
			$DBsIdx[$ii] =~ s/\.fasta$/\.idx/;
			$DBsTestIdx[$ii] =~ s/\.fasta$/\.idx\.kmer_0\.dat/;
			if ( -e "$DBrna//$DBs[$ii]"  && -e "$DBrna//$DBsTestIdx[$ii]"  ){
				next;
			}
			die "\nCould not find expected sortmerna file:\n$srtMRNA_path/rRNA_databases/$DBs[$ii]\n" if ( !-e "$srtMRNA_path/rRNA_databases/$DBs[$ii]"  );
			if ( !-e "$srtMRNA_path/rRNA_databases/$DBsTestIdx[$ii]"  ){
				$DBcmd .= "$srtMRNA_path./indexdb_rna --ref $srtMRNA_path/rRNA_databases/$DBs[$ii],$srtMRNA_path/rRNA_databases/$DBsIdx[$ii]\n";
			}
			$DBcmd .= "\ncp $srtMRNA_path/rRNA_databases/$DBs[$ii] $srtMRNA_path/rRNA_databases/$DBsIdx[$ii]* $DBrna\n";
		}
		#die @DBs."@DBs\n";
		my $ITSDBfa = getProgPaths("ITSdbFA",0);
		if ($ITSDBfa ne ""){ #only do if not empty.. otherwise ignore (not required)
			my $ITSDBpref = $ITSDBfa;$ITSDBpref =~ s/\.fa.*$//;
			my $ITSDBidx = $ITSDBfa; $ITSDBidx =~ s/\.fa.*$/\.idx/;
			$ITSDBpref =~ m/\/([^\/]+)$/;
			$ITSDBpref=~ m/(^.*)\/[^\/]+/;
			
			if (!-e "$ITSDBpref.idx.kmer_0.dat"){ #ITS DBs
				#die "$ITSDBpref\n$$ITSDBpref.idx.kmer_0.dat\n";
				#$DBcmd .= "cp /g/bork3/home/hildebra/DB/MarkerG/ITS_fungi/sh_general_release_30.12.2014.* $DBrna\n";
				if (!-e "$ITSDBfa"){
					print "Missing $ITSDBfa  ITS DB file!\n"; exit(32);
				}
				if (!-e "$ITSDBidx.kmer_0.dat"){
					$DBcmd .= "\n$srtMRNA_path./indexdb_rna --ref $ITSDBfa,$ITSDBidx\n";
				}
				$DBcmd .= "\ncp ${ITSDBpref}* $DBrna\n";
				#has to be noted that this doesn't need to happen again
				#print "ribo DB already present\n";
				
				#$DBcmd = "";
			} 
		}
		#and get flash DBs as well over to that dir
		my @DBn = ("LSUdbFA","LSUtax","SSUdbFA","SSUtax");#,"PR2dbFA","PR2tax"); #"ITSdbFA","ITStax",
		my $LCAar = getProgPaths(\@DBn,0);
		my @LCAdbs = @{$LCAar}; 
#		die "@LCAdbs\n".@LCAdbs."\n";
#			("/g/bork3/home/hildebra/DB/MarkerG/SILVA/123/SLV_123_LSU_sorted_097_centroids.fasta*","/g/bork3/home/hildebra/dev/lotus//DB//SLV_123_LSU.tax",
#			"/g/bork3/home/hildebra/DB/MarkerG/SILVA/123//SLV_123_SSU_sorted_0.97_centroids*", "/g/bork3/home/hildebra/dev/lotus//DB//SLV_123_SSU.tax",
#			"/g/bork3/home/hildebra/DB/MarkerG/SILVA/128/SLV_128_LSU.fa*","/g/bork3/home/hildebra/dev/lotus//DB//SLV_128_LSU.tax",
#			"/g/bork3/home/hildebra/DB/MarkerG/SILVA/128//SLV_128_SSU.fa*", "/g/bork3/home/hildebra/dev/lotus//DB//SLV_128_SSU.tax",
#			"/g/bork3/home/hildebra/DB/MarkerG/ITS_combi/ITS_comb.fa*","/g/bork3/home/hildebra/DB/MarkerG/ITS_combi/ITS_comb.tax",
			#"/g/bork3/home/hildebra/dev/lotus//DB//UNITE/sh_refs_qiime_ver7_99_02.03.2015.fasta", "/g/bork3/home/hildebra/dev/lotus//DB//UNITE/sh_taxonomy_qiime_ver7_99_02.03.2015.txt",
#			"/g/bork3/home/hildebra/DB/MarkerG/PR2/gb203_pr2_all_10_28_99p.fasta*", "/g/bork3/home/hildebra/DB/MarkerG/PR2/PR2_taxonomy.txt");
		#check first if the lambda DB was already built
		my $doCopyDBtoScratch = 0;
		$DBcmd .= "\nmkdir -p $DBrna2\n" if (!-d $DBrna2);
		for (my $kk=0;$kk<@LCAdbs; $kk+=2){
			my $DB = $LCAdbs[$kk];
			#print "XX $DB\n";
			#first test if already copied to scratch..
			$LCAdbs[$kk] =~ m/\/([^\/]+)$/;
			my $SLVtestNme = $1;
			if (!-e $DBrna2."/$SLVtestNme.lambda/index.lf.drp"|| !-e $DBrna2."/$SLVtestNme" ){
				$doCopyDBtoScratch = 1;
			} else {next;}
			#print "$DB\n";
			unless (-e $DB){print "$DBn[$kk] not found, won't use it\n";next;}
			die "wrong DB checked for index built:$DB\n" if ($DB =~ m/\.tax$/);
#			if (!-f $DB.".dna5.fm.sa.val"  ) { #old lambda 1.0x style
			if (!-f "$DB.lba.gz" ){#new lambda3 style  # !-f $DB.".lambda/index.lf.drp"  ) { #new 1.9x style
				print "Building LAMBDA index anew (may take up to an hour)..\n";
				$DBcores = 12;
				#$DBcmd .= "\n$lambdaIdxBin -p blastn -t $DBcores -d $DB\n";
				$DBcmd .= "$lambdaBin mkindexn -t $DBcores -d $DB;\n";
				$DBcmd .= "$pigzBin -p $DBcores $DB.lba;\n";

				#die "$DBcmd\n";
			}
			if ($doCopyDBtoScratch){
				#my $jstr = join("* ",@LCAdbs);$jstr =~ s/\s\*//g;	$jstr =~ s/^\*//g;
				$DBcmd.= "\ncp -r $LCAdbs[$kk]* ".$LCAdbs[$kk+1]."* $DBrna2\n";
			}
		}
		#die "$DBcmd\n$DBrna2\n";
		#die "$SLVlsuNme\n";
		#print "$DBrna2/SLV_128_LSU.tax || !-e $DBrna2/gb203_pr2_all_10_28_99p.fasta  $DBrna2/ITS_comb.fa.dna5.fm.sa.val\n";
		
		my $jN = "_RRDB$JNUM"; my $tmpCmd;
		#die $DBcmd."$DBrna/ITS_comb.idx.kmer_0.dat";
		if ($DBcmd ne ""){
			#die $DBcmd."\n";
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			($jN, $tmpCmd) = qsubSystem($logDir."RiboDBprep.sh",$DBcmd,$DBcores,"10G",$jN,"","",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$MFopt{globalRiboDependence}->{DBcp} = $jN;
		}
	}
	#die;
	
	#first detect LSU/SSU/ITS in metag
	my $tmpDY = "$tmpP/SMRNA/";
	my $cmd = "";
	#die "$MFconfig{readsRpairs}\n";
	my $readConfig = 1;
	my $numCoreL = int($numCore / 2); 
	if (@re1 > 0){
		$cmd .= "\n$cLSUSSUscript -R1 '".join(",",@re1)."' -R2 '". join(",",@re2)."' ";
		if (@singl>0){
			$cmd .= " -RS '".join(",",@singl) . "' "; #tmpP < scratch too slow
		} else {
			$cmd .=" -RS '-1' "; #tmpP < scratch too slow
		}
		$cmd .= "-tmpDir $tmpDY -alignDir $outP -cores $numCore -smplID $SMPN -assmblRibos $MFopt{doRiboAssembl} -DBdir $DBrna\n";#"$tmpDY $outP $numCore $SMPN $MFopt{doRiboAssembl} $DBrna\n\n";
	} else {
		$readConfig = 0; 
		$cmd .= "\n$cLSUSSUscript -R1 '-1' -R2 '-1' -RS '".join(",",@singl)."' -tmpDir $tmpDY -alignDir $outP -cores $numCore -smplID $SMPN -assmblRibos $MFopt{doRiboAssembl} -DBdir $DBrna\n\n"; #tmpP < scratch too slow
	}
	my $sto1 = "$outP/RibFnd.sto";my $stoLCAL = "$outP//ltsLCA/LSU_ass.sto";	my $stoLCAS = "$outP//ltsLCA/SSU_ass.sto";
	$cmd .= "touch $sto1\n" unless ($cmd eq "");
	#this part assigns tax
	my $tmpDX = "$tmpP/LCA/"; 
	my $numCoreL2 = int($numCore2) ;  
	#my $readConfig =$MFconfig{readsRpairs};
	$readConfig=2 if (@singl>0 && $readConfig);
	#ver 2
	#my $cmd2 = "$lotusLCA_cLSU $outP $SMPN $numCoreL2 $DBrna2 $tmpDX $readConfig\n\n";
	#ver 3
	my $cfgstr = "";
	$cfgstr = "-config $MFconfig{configFile} " if ($MFconfig{configFile} ne "" );
	my $lotusLCA_cLSU = getProgPaths("lotusLCA_cLSU_scr");#"perl /g/bork3/home/hildebra/dev/Perl/16Stools/lotus_LCA_blast2.pl";
	my $cmd2 = "";
	$cmd2 .= "rm -rf $tmpDX/*;\nmkdir -p $tmpDX\n\n";
	$cmd2 .= "$lotusLCA_cLSU -dir $outP -smplID $SMPN -cores $numCoreL2 -DBdir $DBrna2 -tmpD $tmpDX -keepReads $MFopt{riboStoreRds} -pairedRds $readConfig $cfgstr -maxReadNum $MFopt{riboLCAmaxRds} -simMode 2 \n\n";
	
	
	#die $cmd2."\n";
	my $jobName="";
	my $Scmd = "";
	my $allLCAstones = 0; 
	$allLCAstones = 1 if ( -e $stoLCAL && -e $stoLCAS );#&& -e "$outP//ltsLCA/ITS_ass.sto");

	if (-d $outP  && -e "$outP/SSU_pull.sto" && -e "$outP/LSU_pull.sto" &&  #&& -e "$outP/ITS_pull.sto"
		-s "$outP//ltsLCA/LSUriboRun_bl.hiera.txt" && $allLCAstones &&
		($MFopt{doRiboAssembl} && -e $outP."/Ass/allAss.sto" ) ){
		#really everything done
		$jobName = $jobd;
	} else {
		$jobd.=";".$MFopt{globalRiboDependence}->{DBcp} unless ($MFopt{globalRiboDependence}->{DBcp} eq "alreadyCopied");
		my $tmpCmd; my $mem = "3G";
		#better to double check calcRiboFind
		my $calcRiboFind=0;$calcRiboFind = 1 if( !-e "$outP/SSU_pull.sto"|| !-e "$outP/LSU_pull.sto" || ($MFopt{doRiboAssembl} && !-e $outP."/Ass/allAss.sto"));
		if (  $calcRiboFind ){ #!-e "$outP/ITS_pull.sto"||
			$jobName = "_RF$JNUM"; 
			#die "RIBOFIND\n$outP/SSU_pull.sto\n"; 
			my $tmpSHDD = $QSBoptHR->{tmpSpace};
			my $curSHFF = int($inputFileSizeMB{$SMPN}/1024*11)+10  ;
			my $predefSHDD = $HDDspace{Ribos}; $predefSHDD =~ s/G$//;
			if ($QSBoptHR->{tmpSpace} < $predefSHDD){ $QSBoptHR->{tmpSpace} = $HDDspace{Ribos};}#overwrite with larger val
			$QSBoptHR->{tmpSpace}= $curSHFF . "G";
			
			#die "SPACE:: $QSBoptHR->{tmpSpace} $jobd\n$outP\n$re1[0]\n";

			($jobName, $tmpCmd) = qsubSystem($logDir."RiboFinder.sh",$cmd,$numCore,$mem,$jobName,$jobd,"",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} = $tmpSHDD; 
		} else {
			$jobName = $jobd;
		}
		if (!-e "$outP//ltsLCA/LSUriboRun_bl.hiera.txt" || !-e "$outP//ltsLCA/SSUriboRun_bl.hiera.txt" || !$allLCAstones ){ #|| !-e "$outP//ltsLCA/ITSriboRun_bl.hiera.txt" 
			$jobd=$jobName; $mem="4G";
			$jobd .= ";".$MFopt{globalRiboDependence}->{DBcp} unless ($MFopt{globalRiboDependence}->{DBcp} eq "alreadyCopied");
			$QSBoptHR->{useLongQueue} = 0;
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = $HDDspace{Ribos};
			($jobName, $tmpCmd) = qsubSystem($logDir."RiboLCA.sh",$cmd2,$numCore2,$mem,"_RA$JNUM",$jobName,"",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$QSBoptHR->{useLongQueue} = 0;
		}
	}
	return $jobName;
}




sub IsDiaRunFinished($){
	my ($curOutDir) = @_;
	my @alldbs = split /,/,$MFopt{reqDiaDB};
	if (!$MFopt{DoDiamond}){return (0,0);}
	my $secCogBin = getProgPaths("secCogBin_scr");
	if ($MFopt{rewriteDiamond} && @alldbs == $MFopt{maxReqDiaDB}){ 
		system "rm -r $curOutDir/diamond/" if (-d "$curOutDir/diamond/");
	}
	my $cD = 0; my $pD = 0; #dia_calc, dia_parse
	if ($MFopt{rewriteDiamond}){$MFopt{redoDiamondParse} = 1;}
	if ($MFopt{redoDiamondParse} && @alldbs == $MFopt{maxReqDiaDB}){
		system "rm -r $curOutDir/diamond/CNT*";
	}
	foreach my $term (@alldbs){
		#print $term."   $cD, $pD\n";
		if ($MFopt{redoDiamondParse} ){#&& ( -e "$curOutDir/diamond/dia.$term.blast.gz.stone" || -e "$curOutDir/diamond/dia.$term.blast.srt.gz.stone") ){ 
			system "rm -f $curOutDir/diamond/dia.$term.blast.*.stone" ;
			system "$secCogBin -i $curOutDir/diamond/XX -DB $term -eval $MFopt{diaEVal} -mode 4";
			system "rm -fr $curOutDir/diamond/ABR/" if ($term eq "ABR");
		}
		if ($MFopt{rewriteDiamond} ){system "rm -f $curOutDir/diamond/dia.$term.blast*" ;$pD=1; $cD=1;}
		#die "$MFopt{rewriteDiamond}\n";
#print "$curOutDir/diamond/dia.$term.blast.gz\n";
		#die "$curOutDir/diamond/dia.$term.blast.gz\n$curOutDir/diamond/dia.$term.blast.srt.gz";
		#|| !-e "$curOutDir/diamond/dia.$term.blast.srt.gz"
		if (!$cD && (!-e "$curOutDir/diamond/dia.$term.blast.gz" && !-e "$curOutDir/diamond/dia.$term.blast.srt.gz" )){$cD = 1; }#system "rm $curOutDir/diamond/dia.$term.blas*.gz";}
		$pD = 1 if (!-e "$curOutDir/diamond/dia.$term.blast.srt.gz.stone");#  <- last version always requires .srt.gz
		#print "$cD, $pD  $curOutDir/diamond/dia.$term.blast.gz\n";
	}
	#die "$cD, $pD\n";
	$cD = 0 if ($pD==0);
	
	if ($MFopt{rewriteAllIfAnyDiamond} && ($cD  || $pD)){ #just delete everything..
		system "rm -r $curOutDir/diamond/" if (-d "$curOutDir/diamond/"); $cD=1; $pD=1;
	}
	if (!$cD && !$pD){
		foreach my $curDB (@alldbs){$progStats{$curDB}{SearchCompl}++;}
	}
	
	return ($cD,$pD);
}

sub prepDiamondDB($ $ $ $){#takes care of copying the respective DB over to scratch
	my ($curDB,$CLrefDBD,$ncore, $searchMode) = @_;
	#searchMode 1: diamond, 2: mmseqs2
	my $diaBin = getProgPaths("diamond");
	my $mmseqs2Bin = getProgPaths("mmseqs2");
	my $dbSuffix = ".db.dmnd";

	my ($DBpath ,$refDB ,$shrtDB) = getSpecificDBpaths($curDB,0);
	system "mkdir -p $CLrefDBD" unless (-d $CLrefDBD);
	my $clnCmd = "";
	if ($MFopt{globalDiamondDependence}->{$curDB} eq "" ){
		my $DBcmd = ""; 
		# check for pr.db.dmnd
		my $ncoreDB=1;
		my $epoTS = 0; my $refDBnew=0;
		my $epoTSdia = 99999999999999;
		
		if ($searchMode == 1 ){
			$epoTSdia = ( stat "$diaBin" )[9] if (-f $diaBin);
		} else {
			$epoTSdia = ( stat "$mmseqs2Bin" )[9];
			$dbSuffix = ".db.mms2";
		}
		$epoTS = ( stat "$DBpath$refDB${dbSuffix}" )[9] if (-e "$DBpath$refDB${dbSuffix}");
		#$timestamp = POSIX::strftime( "%d%m%y", localtime( $epoTS));
		#print "time: $epoTSdia $epoTS \n";
		if (!-e "$DBpath$refDB${dbSuffix}"){
		# age check completely deactivated..
		#if ($epoTS < $epoTSdia){#checks age of binary vs DB creation.. 
			$ncoreDB = $ncore;$refDBnew=1; 
			system "rm -f $DBpath$refDB${dbSuffix}*";
			if ($searchMode == 1){
				$DBcmd .= "$diaBin makedb --in $DBpath$refDB -d $DBpath$refDB${dbSuffix} -p $ncoreDB\n";
			} else {
				$DBcmd .= "$mmseqs2Bin createdb $DBpath$refDB $DBpath$refDB${dbSuffix} --compressed 1\n";
			}
		}
		unless (-e "$DBpath$refDB.length"){
			my $genelengthScript = getProgPaths("genelength_scr");#
			$DBcmd .= "$genelengthScript $DBpath$refDB $DBpath$refDB.length\n";
		}
		#$clnCmd .= "rm -rf $CLrefDBD;" if (length($CLrefDBD)>6);
		#idea here is to copy to central hdd (like /scratch)
		system "rm -f $CLrefDBD/$refDB${dbSuffix}*"  if (($refDBnew || $MFopt{rewriteDiamond} )&& -d $CLrefDBD);
		#print " !-e $CLrefDBD/$refDB.db.dmnd && !-e $CLrefDBD/$refDB.length\n ";
		if ( -e "$CLrefDBD/$refDB${dbSuffix}" && -e "$CLrefDBD/$refDB.length" 
			#&& ($curDB eq "NOG" && !-e "$CLrefDBD/NOG.members.tsv") &&
			#(($curDB ne "KGB" && $curDB ne "KGM" && $curDB ne "KGE")|| -s "$CLrefDBD/genes_ko.list")
			){
			#has to be noted that this doesn't need to happen again
			$MFopt{globalDiamondDependence}->{$curDB}="$shrtDB-1";
		} else {
			$DBcmd .= "mkdir -p $CLrefDBD\n";
			$DBcmd .= "cp $DBpath$refDB${dbSuffix}*  $DBpath$refDB.length  $CLrefDBD\n";
		}
		
		#specialized file copy
		if ($curDB eq "NOG" && !-s "$CLrefDBD/NOG.members.tsv" && !-s "$CLrefDBD/NOG.annotations.tsv"){
			system "rm -f $CLrefDBD/NOG* $CLrefDBD/all_species_data.txt";
			$DBcmd .= "cp $DBpath/all_species_data.txt $DBpath/NOG.members.tsv $DBpath/NOG.annotations.tsv $CLrefDBD\n";
		}
		if ($curDB eq "CZy" && !-s "$CLrefDBD/MohCzy.tax"){
			system "rm -f $CLrefDBD/MohCzy.tax $CLrefDBD/cazy_substrate_info.txt";
			$DBcmd .= "cp $DBpath/MohCzy.tax $DBpath/cazy_substrate_info.txt $CLrefDBD\n";
		}
		if (($curDB eq "KGB" || $curDB eq "KGM" || $curDB eq "KGE") && !-s "$CLrefDBD/genes_ko.list"){
			system "rm -f $CLrefDBD/genes_ko.list $CLrefDBD/kegg.tax.list";
			$DBcmd .= "cp $DBpath/genes_ko.list $DBpath/kegg.tax.list $CLrefDBD\n";
		}
		if ($curDB eq "ABRc" && !-s "$CLrefDBD/card.parsed.f11.tab.map"){ 
			system "rm -f $CLrefDBD/card*";
			$DBcmd .= "cp $DBpath/card*.txt $DBpath/card*.map $CLrefDBD\n";
		}
		if ($curDB eq "PTV" && !-s "$CLrefDBD/PATRIC_VF2.tab"){ 
			system "rm -f $CLrefDBD/PATRIC_VF2.tab";
			$DBcmd .= "cp $DBpath/PATRIC_VF2.tab $CLrefDBD\n";
		}
		if ($curDB eq "VDB" && !-s "$CLrefDBD/VF.tab"){ 
			system "rm -f $CLrefDBD/VF.tab";
			$DBcmd .= "cp $DBpath/VF.tab $CLrefDBD\n";
		}
		if ($curDB eq "PAB" && !-s "$CLrefDBD/all_species_data.txt"){
			#copy NOG taxonomy
			my ($DBpathN ,$refDBN ,$shrtDBN) = getSpecificDBpaths("NOG",0);
			$DBcmd .= "cp $DBpathN/all_species_data.txt $CLrefDBD\n" unless (-e "$CLrefDBD/all_species_data.txt");
		}
		if ($curDB eq "TCDB" && !-s "$CLrefDBD/hir.txt"){ 
			$DBcmd .= "cp $DBpath/TCDBhir.txt  $CLrefDBD\n";
		}

		my $jN = "_DIDB$shrtDB$JNUM"; my $tmpCmd;
#		die "$DBcmd";
		if ($DBcmd ne ""){
			my @preConstr = @{$QSBoptHR->{constraint}};
			#push(@{$QSBoptHR->{constraint}}, "intel");
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			($jN, $tmpCmd) = qsubSystem($logDir."DiamondDBprep$shrtDB.sh",$DBcmd,$ncoreDB,int(100/$ncoreDB)."G",$jN,"","",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$MFopt{globalDiamondDependence}->{$curDB} = $jN;
			@{$QSBoptHR->{constraint}} = @preConstr;
			#die "$MFopt{globalDiamondDependence}->{$curDB}";
		}
		#die $DBcmd."\n";
	}
	return ($refDB,$shrtDB,$clnCmd);
}


sub getRdLibs($ $ $ $){
	my ($ar1,$ar2,$sa1,$mrgHshHR) = @_;
	my @reads1 = @{$ar1}; my @reads2 = @{$ar2}; my @singlRds = @{$sa1};
	my %mrgHsh = %{$mrgHshHR};
	my $mrgMode = 0;
	$mrgMode =1 if (exists ($mrgHsh{mrg}) && $mrgHsh{mrg} ne "");
	my $useLibArrays = 3;
	if ($mrgMode) {$useLibArrays=4;}
	my %ret;
	#print "mrgMode : : $mrgMode $useLibArrays \n@reads1\nYTY\n";
	for (my $kk=0;$kk<$useLibArrays;$kk++){
		my @rds = @singlRds;
		if ($mrgMode){
			if ($kk==1){@rds = $mrgHsh{pair2};}
			if ($kk==2){@rds = $mrgHsh{pair1};}
		} else {
			if ($kk==1){@rds = @reads2;}
			if ($kk==2){@rds = @reads1;}
		}
		if ($kk==3){@rds = $mrgHsh{mrg};}
		$ret{$kk} = \@rds;
	}
	return %ret;
}

sub isLastSampleInAssembly{
	my ($assD,$curD) = @_;
	return 0 if (!-d $assD);
	return 0 if (!-e "$assD/smpls_used.txt");
	open I,"<$assD/smpls_used.txt" or die $!;
	my $lastEntr="";
	while (<I>){
		chomp $_; 
		$lastEntr = $_ if (length($_)>2);
	}
	close I;
	#print "$lastEntr  $curD \n";
	return 1 if ($lastEntr eq $curD);
	return 0;
}


#preparation of secondary mapping, including wildcard resolution, gene calling, index building
sub map2ndPrep{
	unless ($ARGV0 eq "scaffold" || $ARGV0 eq "map2tar" || $ARGV0 eq "map2DB" || $ARGV0 eq "map2GC"){
		#not asked for 2nd map? ok, deactivate all related parameters
		$MFopt{mapModeTogether} = 0; $MFopt{DoMapModeDecoy} = 0; $MFopt{mapModeActive} = 0;
		return;
	}
	if ($ARGV0 eq "scaffold"){
		die"update scaffold\n";
		$scaffTarExternal = $ARGV[1];
		if (!-f $scaffTarExternal){
			die "Could not find scaffold file:\n$scaffTarExternal\n";
		}
		$scaffTarExternalName = $ARGV[2];
		if (@ARGV>3){
			$scaffTarExtLibTar = $ARGV[3];
		}
		return;
			

	}

#in this case primary focus is on mapping and not on assemblies
	if ($ARGV0 eq "map2DB" || $ARGV0 eq "map2GC"){
		$MFopt{mapModeCovDo}=0;$MFopt{DoMapModeDecoy}=0;$MFopt{mapModeTogether}=0;$map2ndMpde=2;
	}
	if ($MFopt{map2Assembly} || $MFopt{DoAssembly}){
		print "Mapping mode: reference ";
		print " Decoy" if ($MFopt{DoMapModeDecoy});
		print " Competitive" if ($MFopt{mapModeTogether}== 1);
		print " combined map, seperate reporting" if ($MFopt{mapModeTogether}== 2);
		print " combined map, combined reporting" if ($MFopt{mapModeTogether}== -1);
		print "\nDeactivating assembly and dependent modules.\n";
		$MFopt{map2Assembly}=0; $MFopt{DoAssembly} =0;
	}
	my $DBsubmCnt=0; my $GENEsubmCnt=0;
	my @refDB1 = split(/,/,$MFopt{refDBall});
	#$MFopt{mapModeTogether} = 0 if (@refDB1 == 1);#could still be a wildcard, wrong assumption
	my @bwt2Name1;
	if (defined $MFopt{bwt2NameAll}){
		@bwt2Name1 = split(/,/,$MFopt{bwt2NameAll}) ;
	} elsif ($ARGV0 eq "map2DB"  ){
		@bwt2Name1 = ("refDB");
	} elsif ($ARGV0 eq "map2GC"){
		@bwt2Name1 = ("GC");
		$map2ndMpde=3;
	} else {
		@bwt2Name1 = ("auto") x scalar(@refDB1);
	}
	my @refDB;my @bwt2Name ;
	my %FNrefDB2ndmap;
	
	for (my $i=0;$i<@refDB1;$i++){
		my @sfiles;
		if ($ARGV0 eq "map2GC"){
			die "can only have one ref to GC: @refDB1\n" unless (@refDB1 == 1);
			@sfiles = ($refDB1[0]."/compl.incompl.95.fna");
			die "Could not find reference in GC dir. Expected: $sfiles[0]\n" unless (-e $sfiles[0]);
		} else {
			@sfiles = glob($refDB1[$i]);
		}
		#die "@sfiles\n$refDB1[$i]\n".@sfiles."\n";
		my $iniBwtNm = $bwt2Name1[$i];
		#die "@sfiles\n$refDB1[$i]\n";
		if (@sfiles>1){
			for (my $j=0;$j<@sfiles;$j++){
				push(@refDB,$sfiles[$j]);
				if ($iniBwtNm eq "auto"){
					#$sfiles[$j] =~ m/ssemblyfind_list_(.*)\.txt\/(.*)\.contigs_/; my $nmnew = $1.$2; $nmnew =~ s/#/_/g;
					$sfiles[$j] =~ m/\/([^\/]+)\.f.*a$/;
					my $nmnew = $1;
					#die "$nmnew\n";
					push(@bwt2Name,$nmnew);

				} else {
					push(@bwt2Name,$bwt2Name1[$i].$j);
				}
			}
		} elsif (@sfiles == 1) {
			push(@refDB,$sfiles[0]);
			if ($iniBwtNm eq "auto"){
				$sfiles[0] =~ m/.*\/([^\/]+)$/;
				$iniBwtNm = $1;
				$iniBwtNm =~ s/\.[^\.]+$//;
			}
			push(@bwt2Name,$iniBwtNm);
		} else {
			die "Could not find file for entry $refDB1[$i]\n";
		}
	}
	#die "@refDB\n@bwt2Name\n";
	my $shrtMapNm = "";	$shrtMapNm = "Comb_" if (@bwt2Name > 1);
	for ( my $i=0;$i< @bwt2Name; $i++){
		my $substrl = length($bwt2Name[$i]);	if (@bwt2Name > 3){$substrl = 8;}	
		if (@bwt2Name > 5){$substrl = 6;}	if (@bwt2Name > 7){$substrl = 3;}
		if ($i==0){
			$shrtMapNm .= substr($bwt2Name[$i],0,$substrl);
		} else {
			$shrtMapNm .= ".". substr($bwt2Name[$i],0,$substrl);
		}
		if ($i>5){$shrtMapNm .= ".X".(@bwt2Name - $i)."X"; last;}
	}
	$map2ndTogRefDB{DB} = "$baseOut/GlbMap/$shrtMapNm/$shrtMapNm.fa";
	#die"$map2ndTogRefDB{DB}\n";
	#decoy mapping setup (only required in map2tar
	$make2ndMapDecoy{Lib} = "";
	#die "decoy mapping not ready for multi fastas\n" if (@refDB > 1);
	if ($MFopt{DoMapModeDecoy} || $MFopt{mapModeTogether}){
		if ($MFopt{mapModeTogether} && ($rewrite || $MFopt{MapRewrite2nd}) ){
			#system("rm -r -f $map2ndTogRefDB{DB}*");#;mkdir -p $bwt2outDl
		}
		for (my $i=0;$i<@refDB; $i++){ #take care of ref DB decoy prep
			#last if (-e $map2ndTogRefDB{DB});
			#print "$refDB[$i]\n";
			#die "XX\n";
			my $aref;
			if ($MFopt{mapModeTogether}){#read into mem & combine
				my $hr = readFasta($refDB[$i],1); 
				$hr = prefixFAhd($hr,$bwt2Name[$i]);#rename header to fasta files...
				my %FN = %{$hr};
				%FNrefDB2ndmap = (%FNrefDB2ndmap , %FN);
				$aref = [keys %FN];
			} else {
				$aref = readFastHD($refDB[$i]);
			}
			if ($MFopt{mapModeTogether}>0 || $MFopt{DoMapModeDecoy}){
				push(@{$make2ndMapDecoy{regions}}, join(" ",@{$aref}) );
				push(@{$make2ndMapDecoy{region_lcs}},  lcp(@{$aref}) );#prefix_find($aref)  );
			}
		}
	}
	#		die "@{$make2ndMapDecoy{region_lcs}}\n";
	
	#build of combined refDB, if competitive mapping..
	my $bwtDBcore = 10; 
	
	if ($MFopt{mapModeTogether}){#the fasta's were already combined into %FNrefDB2ndmap, just built idx now..
		system "mkdir -p $baseOut/GlbMap/LOGandSUB/" unless (-d "$baseOut/GlbMap/LOGandSUB/");
		system "mkdir -p $runTmpDBDirGlobal/$shrtMapNm" unless (-d "$runTmpDBDirGlobal/$shrtMapNm");
		system "mkdir -p $baseOut/GlbMap/$shrtMapNm" unless (-d "$baseOut/GlbMap/$shrtMapNm");
		#system "mkdir -p $map2ndTogRefDB{DB}" unless (-d $map2ndTogRefDB{DB});
		print "$map2ndTogRefDB{DB}\n";
		#die;
		writeFasta(\%FNrefDB2ndmap,"$map2ndTogRefDB{DB}") unless (-e $map2ndTogRefDB{DB} && -s $map2ndTogRefDB{DB});
		if ($MFopt{mapModeTogether}==-1){
			@refDB = ($map2ndTogRefDB{DB});
			@bwt2Name = ($shrtMapNm);
			$MFopt{mapModeTogether} = 0;#deactivate, as from now will be treated as singular ref
		} else {
			my ($cmd,$DBbtRef, $chkFile) = buildMapperIdx($map2ndTogRefDB{DB},$bwtDBcore,$MFopt{largeMapperDB},$MFopt{MapperProg}) ;
			if (!-e $chkFile ){
				my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
				($bwt2ndMapDep,$cmd) = qsubSystem("$baseOut/GlbMap/LOGandSUB/builBwtIdx_comp.sh",$cmd,$bwtDBcore,(int(25/$bwtDBcore)+1) ."G","BWI_compe","","",1,[],$QSBoptHR) ;
				$QSBoptHR->{tmpSpace} =$tmpSHDD;
			}
		}
	}
	
	print "\n=======================\nmap to $MFopt{refDBall}\n with map mode $MFopt{mapModeTogether}\n=======================\n\n";
	$MFopt{mapModeActive} =1; my $cmdBIG = "";my $cmdBIGgene = "";
	#die "@refDB\n";
	
	#build index for each fasta, and predict genes on these
	for (my $i=0;$i<@refDB; $i++){
		#$refDB[$i] =~ m/(.*\/)[^\/]+/;
		#my $refDir = $1;
		my $bwt2outDl = "$baseOut/GlbMap/$bwt2Name[$i]/";
		if ($map2ndMpde == 3){#outdir should be in the GC dir
			$bwt2outDl = "$refDB1[0]/unmappedMap/"
		}
		push @bwt2ndMapNmds , $bwt2Name[$i];
		push(@bwt2outD,$bwt2outDl);
		system "mkdir -p $bwt2outDl" unless (-d $bwt2outDl);
		if ($MFopt{mapModeTogether} >= 0 && ($rewrite || $MFopt{MapRewrite2nd}) ){
			print "Deleting previous mapping DBs..\n" if ($i==0);	$refDB[$i] =~ m/.*\/([^\/]+)$/;
			system("rm -r -f $runTmpDBDirGlobal/$1*");#;mkdir -p $bwt2outDl
		}
		
		$refDB[$i] =~ m/.*\/([^\/]+)$/;
		system "cp $refDB[$i] $bwt2outDl" if ($MFopt{mapModeCovDo}  && !-e "$bwt2outDl/$1");
		#print "\n$refDB[$i]\n";
		#die "$bwt2outDl/$1\n";
		#system "mkdir -p $bwt2outDl/LOGandSUB" unless (-d "$bwt2outDl/LOGandSUB");
		
		
		
		my ($cmd,$DBbtRef,$chkFile) = buildMapperIdx($refDB[$i],$bwtDBcore,$MFopt{largeMapperDB},$MFopt{MapperProg}) ;
		$DBbtRef =~ s/$MFcontstants{bwt2IdxFileSuffix}$//;
		$DBbtRef =~ m/.*\/([^\/]+)$/;
		$DBbtRef = "$runTmpDBDirGlobal/$1";#set up to scratch dir to map onto
		my $idxNFini = 0; $idxNFini = 1 if (!-e $chkFile); #mapperDBbuilt($DBbtRef,$MFopt{MapperProg}); #($MFopt{MapperProg}==3 && !-e "$DBbtRef.pak") || ($MFopt{MapperProg}==1 && !-e "$DBbtRef$MFcontstants{bwt2IdxFileSuffix}.rev.1.bt2");
		$cmd.= "\ncp $refDB[$i]* $runTmpDBDirGlobal\n" if ($idxNFini);# if (!$MFopt{mapModeCovDo} && !-e "$bwt2outDl/$1");
		#print $cmd."\n";
		#die $DBbtRef."\n$runTmpDBDirGlobal/\n";
		if (!$MFopt{mapModeTogether} && $idxNFini){ #not required for these map modi
			#system $cmd 
			#	my ($bwt2ndMapDep2,$cmd2) = qsubSystem($bwt2outDl."/LOGandSUB/builBwtIdx$i.sh",$cmdBIG,$bwtDBcore,(int(20/$bwtDBcore)+1) ."G","BWI".$i,"","",1,[],$QSBoptHR) ;$bwt2ndMapDep .= ";$bwt2ndMapDep2";
			$cmdBIG .= "\n\n#====== $i =======\n".$cmd;
			$DBsubmCnt++;
		}
		
		#$DBbtRefX = $DBbtRef;
		#die "$DBbtRef\n";
		push(@DBbtRefX,$DBbtRef);
		if($MFopt{mapModeCovDo} && $map2ndMpde != 3){ #get the coverage per gene etc; for this I need a gene prediction
												#but not for GC mapping (these are genes already)
			my $gDir = $bwt2outDl."";
			my $nativeGFF = $refDB[$i];$nativeGFF =~ s/\.[^\.]+$/\.gff/;
			my $gffF = "genes.$bwt2Name[$i].gff";
			#die "$nativeGFF\n";
			if (-e "$gDir/$gffF"){
				;
			}elsif (-e $nativeGFF){
				system "cp $nativeGFF $gDir/$gffF";
			} else {
				system "mkdir -p $gDir";
				$logDir = $gDir;
				my $dEGP = $MFopt{DoEukGenePred}; $MFopt{DoEukGenePred} = 0;
				my $tmpDep1 = genePredictions($refDB[$i],$gDir,"",$gDir,"iGP$i","",0);
				$MFopt{DoEukGenePred} = $dEGP;
				$cmdBIGgene .= "#====== $i =======\n".$tmpDep1."cp $gDir/genes.gff $nativeGFF; mv $gDir/genes.gff $gDir/$gffF\n\n\n";
				$GENEsubmCnt++;
				#my ($tmpDep,$tmpCmd) = qsubSystem( $bwt2outDl."/LOGandSUB/cpGenes.sh",  "cp $gDir/genes.gff $nativeGFF; mv $gDir/genes.gff $gDir/$gffF\n",
				#1,"1G","genecop".$i,$tmpDep1,"",1,[],$QSBoptHR);
				#$bwt2ndMapDep .= ";".$tmpDep;
			}
			push(@DBbtRefGFF,$gDir."/$gffF");#"genePred/genes.gff"
		}
		$logDir="";
	} #end for loop building bwtIdx
	#now submit all together as single call..
	my $bwt2outDl = "$baseOut/GlbMap/LOGandSUB/"; system "mkdir -p $bwt2outDl" unless (-d $bwt2outDl);
	#submit mapping index build
	#die "$cmdBIG\n\n";
	if ($DBsubmCnt>0){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		my $mapperMemDB= 20; $mapperMemDB = 40 if ($MFopt{largeMapperDB});
		my ($bwt2ndMapDep2,$cmd2) = qsubSystem($bwt2outDl."/builBwtIdxBIG.sh",$cmdBIG,$bwtDBcore,(int($mapperMemDB/$bwtDBcore)+1) ."G","BWIbig","","",1,[],$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$bwt2ndMapDep .= ";$bwt2ndMapDep2";
	}
	#submit gene predictions
	if ($GENEsubmCnt>0){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		my ($tmpDep,$tmpCmd) = qsubSystem( $bwt2outDl."/GenesPredBIG.sh",  $cmdBIGgene,
			1,"1G","genePred","","",1,[],$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$bwt2ndMapDep .= ";".$tmpDep;
	}
	
	
	if ($map2ndMpde == 3 && !$MFopt{useUnmapped}){
		die "Did you mean to activate \"-mapUnmapped\"?\n";
	}

	#die "@bwt2outD\n";
}


sub runOrthoPlacement(){
	my ($cleanSeqSetHR ,$outD,$tmpP,$jdep) = @_;
	die "runOrthoPlacement no longer active\n";
	my $ar1 = ${$cleanSeqSetHR}{arp1}; my $ar2 = ${$cleanSeqSetHR}{arp2}; my $sa1 = ${$cleanSeqSetHR}{singAr}; 
	my $mrgHshHR = ${$cleanSeqSetHR}{mrgHshHR};
	
	my $sdmBin = getProgPaths("sdm");#"/g/bork3/home/hildebra/dev/C++/sdm/./sdm";
	my $fna2faaBin = getProgPaths("fna2faa");
	system "mkdir -p $outD $tmpP" unless (-d $outD && -d $tmpP);
	my %RdLibs = getRdLibs($ar1,$ar2,$sa1,$mrgHshHR);
	my $numCore=22;
	my $scrP = "";my $cmd="";
	$cmd .= "mkdir -p $tmpP\n";
	#my $tmpF = ${$RdLibs{$kk}}[0];
	system "rm -f $tmpP/6fr.fna";
	for (my $kk=0;$kk<scalar(keys(%RdLibs));$kk++){
		my @rds = @{$RdLibs{$kk}};
		$scrP = $rds[0] if ($scrP eq "");
		foreach my $fnaI (@rds){
			$cmd .= "$sdmBin -i $fnaI -o_fna - | $fna2faaBin -q - >> $tmpP/6fr.fna\n";
		}
	}
	$scrP =~ s/\/[^\/]+$//;
	#die $scrP."\n";
	my $hmmscr = "python /g/bork3/home/hildebra/dev/Perl/SoilHelpers/extract_domains.py";
	my $hmmD = "/g/bork3/home/hildebra/DB/HMMs/FungiAB/";
	my @hmmsDB=("Condensation.hmm","dmat.hmm","AMP-binding.hmm","PKS_KS.hmm","Terpene_synth_C.hmm");
	my @hmmNms = ("Condensation","dmat","AMP","PKS","Terpene");
	my @preMSAs = (); #ref alignments, big DB
	my @refTrees = (); #ref trees, calc with fastree
	my $onceExtr=0;
	my $redo=0;
	for (my $i=0; $i<@hmmsDB; $i++){	
		my $hmmM  = "$hmmD/$hmmsDB[$i]";
		my $hmmName = $hmmsDB[$i]; $hmmName=~s/\.hmm//;
		if (!-e "$outD/$hmmName.fna" || $redo){
			$cmd .= "$hmmscr $tmpP/6fr.fna $hmmM $outD/$hmmName.fna $numCore\n" ;
			$onceExtr=1;
		}
	}
	if($onceExtr ==0 && !$redo){
		$cmd = "";
	}
	$cmd .= "mkdir -p $tmpP\n";
		
	#second new part .. alignment to existing tree
	$onceExtr=0;
	for (my $i=0; $i<@hmmsDB; $i++){	
		my $hmmName = $hmmsDB[$i]; $hmmName=~s/\.hmm//;
		my $hmmN2 = $hmmNms[$i];
		my $tarFile = "$outD/$hmmName.fna";
		if (-e $tarFile && -s "$tarFile" == 0){next;}
		#$tarFile = fixHDs4Phylo($tarFile);
		#fix fasta headers
#		$cmd .= "cat $tarFile | perl -p -e 's/[:|#+]/_/g' > $tarFile\n";

		my $clustaloBin = getProgPaths("clustalo");
		my $raxMLbin = getProgPaths("raxml");
		my $treeDBdir = "/g/bork3/home/hildebra/data/SoilABdoms/Jaime/JaimeT2/";
#		my $refTREE = "$treeDBdir/$hmmN2/final_alg/phylo/FASTTREE_allsites.nwk";
		my $refTREE = "$treeDBdir/$hmmN2/final_alg/phylo/IQtree_fast_allsites.treefile";
		my $refMSA = "$treeDBdir/$hmmN2/final_alg/outMSA.faa";
		
		if (!-e "$outD/RAxML_classification.${hmmName}_place"){
			$onceExtr=1;
			$cmd .= "$clustaloBin --p1 $refMSA -i $tarFile --threads $numCore > $tmpP/$hmmName.msa.faa\n";
			$cmd .= "sed -i 's/[:|#+]/_/g' $tmpP/$hmmName.msa.faa\n";
			$cmd .= "rm -f $tmpP/RAxML*\n";
			$cmd .= "$raxMLbin -f v -s $tmpP/$hmmName.msa.faa  -w $tmpP -t $refTREE -m PROTGAMMALG -n ${hmmName}_place -T$numCore\n";
			#$cmd .= "mv $tmpP/RAxML_fastTreeSH_Support.${hmmName}_place $outD/${hmmName}_place.nwk\n";
			$cmd .= "mv $tmpP/*${hmmName}_place $outD/\n";
		}
	}
	
		#die $cmd;
	
	
	$cmd .= "rm -r $tmpP\n";
	#die $cmd;
	my $jobName = "_ORTH$JNUM"; 
	if ($onceExtr){
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."PABpred.sh",$cmd,$numCore,"4G",$jobName,$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		return $jobName2;
	} else { return "";}
}




sub runDiamond(){
	my ($cleanSeqSetHR,$outD,$CLrefDBD,$tmpP,$jdep,$curDB_o) = @_;
	my $ar1 = ${$cleanSeqSetHR}{arp1}; my $ar2 = ${$cleanSeqSetHR}{arp2}; my $sa1 = ${$cleanSeqSetHR}{singAr}; my $mrgHshHR = ${$cleanSeqSetHR}{mrgHshHR};
	
	#die "runDiamond:: $mrgHshHR\n";
	
	my $secCogBin = getProgPaths("secCogBin_scr");
	my $diaBin = getProgPaths("diamond");
	my $mmseqs2Bin = getProgPaths("mmseqs2");
	my $searchMode = 1;
	
	my %RdLibs = getRdLibs($ar1,$ar2,$sa1,$mrgHshHR);
	
	die "No reads found for sample $outD in runDiamond sub\n" if (scalar(keys(%RdLibs)) == 0);
	
	my $clnCmd="";my $jobN2="";
	#die;
	system "mkdir -p $outD" unless (-d $outD && -d $tmpP);
	#my $diaOfmt = "tab";#old
	#$Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc
	my $diaOfmt = "-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "; #diamond
	
	my $mmsOfmt = "--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"; #mmseqs
	my $ncore = $MFopt{diaCores}; my $sensBlast = "";
	$sensBlast = " --sensitive " if ($MFopt{diaRunSensitive});
	foreach my $curDB (split /,/,$curDB_o){
		$progStats{$curDB}{SearchCompl} = 0 unless (defined($progStats{$curDB}{SearchCompl}));
		$progStats{$curDB}{SearchIncomplete} = 0 unless (defined($progStats{$curDB}{SearchIncomplete}));
		#print "$curDB";
		my ($refDB,$shrtDB,$clnCmd) = prepDiamondDB($curDB,$CLrefDBD,$ncore,$searchMode);
		my $doInterpret = 1;
		#print "$refDB ,$shrtDB,$clnCmd\n";
		#$doInterpret = 0 if ($shrtDB eq "ABR");
		my $getQSeq = 0;
		$getQSeq = 0 if ($curDB eq "PAB");
		my $diaOfmt2 = $diaOfmt;
		$diaOfmt2 .= " qseq" if ($getQSeq); #in case, I want to get the query sequence (matching)
		#run actual diamond
		my @collect = (); my @collectSingl=();
		my $cmd ="mkdir -p $tmpP\n";
		for (my $kk=0;$kk<scalar(keys(%RdLibs));$kk++){
			my $rdsUsed = 0;
			my @rds = @{$RdLibs{$kk}};
			for (my $ii=0;$ii<@rds;$ii++){
				my $query = $rds[$ii];	
				my $outF = "$tmpP/DiaAssignment.sub.$shrtDB.$kk.$ii";
				#my $tmpcnt = `grep -c '^>' $query`; chomp $tmpcnt
				#--comp-based-stats 0
				if ($searchMode==1){
				$cmd .= "$diaBin blastx $diaOfmt2 --masking 0 --comp-based-stats 0 --compress 1 --quiet -t $tmpP --min-orf 25 -d $CLrefDBD$refDB.db -q $query -k 5 -e 1e-4 -o $outF $sensBlast -p $ncore\n"; #
				} elsif ($searchMode==2){
				$cmd .= "$mmseqs2Bin easy-search $query $CLrefDBD$refDB.db.mms2 $outF.gz $tmpP --threads $ncore --max-accept 500 --compressed 1 -s 4 $mmsOfmt \n";
				}
				
				#$cmd .= "$diaBin view -a $outF.tmp -o $outF -f tab\nrm $outF.tmp.daa\n";
				if ($kk==0 || $kk == 3){#single or ext fragments, doesn't need to be sorted
					push(@collectSingl,$outF.".gz");
				} else {
					push(@collect,$outF.".gz");
				}
			}
		}
		
		#die "$cmd\n";
		
		my $out = $outD."dia.$shrtDB.blast";
		my $outgz = "$out.srt.gz";
		#die "$outgz\n";
		#unzip, sort, zip
		$cmd .= "zcat ".join( " ",@collect) ." | sort -t\$'\\t' -k1 -T $tmpP | $pigzBin --stdout -p $ncore > $outgz \n";
		$cmd .= "rm -f ". join( " ",@collect) . "\n";
		if (@collectSingl >= 1){
			#append on gzip, can be done with gzip
			$cmd .= "cat ".join( " ",@collectSingl) ." >> $outgz\nrm -f " . join( " ",@collectSingl) ."\n"; #$out.srt
		}
		$cmd.= "rm -r $tmpP\n";
		#die $cmd."\n";
		my $cmd2 = "$secCogBin -i $outgz -DB $shrtDB -eval $MFopt{diaEVal} -percID $MFopt{DiaPercID} -minAlignLen $MFopt{DiaMinAlignLen} -minFractQueryCov $MFopt{DiaMinFracQueryCov} -mode 0 -LF $CLrefDBD/$refDB.length -reportDomains $getQSeq -DButil $CLrefDBD -tmp $tmpP";
		#$cmd2 .= " " if ($getQSeq);
		if ($curDB eq "ABR"){
			my $KrisABR = getProgPaths("KrisABR_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/ABRblastFilter.pl";
			$cmd2 = "$KrisABR $outgz $outD/ABR/ABR.genes.txt $outD/ABR/ABR.cats.txt\n";
		} elsif ($curDB eq "PAB" && $MFopt{PABtaxChk}){ #NOG assignments
			$cmd2 .= " -NOGtaxChk $outD/dia.NOG.blast.srt ";
		}
		$cmd2 .= "\n";
		$cmd2 .= "rm -f $outgz\n" if ($MFopt{DiaRmRawHits});
		
		#check if the secondary routines for parsing blast out still need to be run
		if ($doInterpret){
			if (!(-e  "$out.gz.stone" || -e  "$out.srt.gz.stone" || -e  "$out.stone")){$doInterpret=1;
			} else {$doInterpret=0;
			}
		}
		my $jobName = $jdep;
		my $globDep = $MFopt{globalDiamondDependence}->{$curDB};
		$globDep = "" if ($MFopt{globalDiamondDependence}->{$curDB} eq "$shrtDB-1");
		
		my $memu = "7G";my $tmpCmd;
		if (!-d $outD || !(-e "$out" || -e "$out.gz"|| -e "$out.srt.gz") ){ #diamond alignments
			$jobName = "_D$shrtDB$JNUM"; 
			my @preConstr = @{$QSBoptHR->{constraint}};
			push(@{$QSBoptHR->{constraint}}, $avx2Constr);
			my $tmpSHDD = $QSBoptHR->{tmpSpace};
			$QSBoptHR->{tmpSpace} = $HDDspace{diamond}; #set option how much tmp space is required, and reset afterwards
			($jobName, $tmpCmd) = qsubSystem($logDir."Diamo$shrtDB.sh",$cmd,$ncore,$memu,$jobName,$jdep.";".$globDep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
			@{$QSBoptHR->{constraint}} = @preConstr;
			$QSBoptHR->{tmpSpace} = $tmpSHDD;
			$doInterpret=1;#run interpret step in any case
			
		} else {
			$jobName = $globDep;
		}
		if ($doInterpret){ #parsing of dia output
			my $jobName2 =  "_DP$shrtDB$JNUM";
			$memu = "30G";
			($jobName, $tmpCmd) = qsubSystem($logDir."Diamo_parse$shrtDB.sh",$cmd2,1,$memu,$jobName2,$jobName,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
			#die "bo";
			$progStats{$curDB}{SearchIncomplete}++;#this needs to be done first..
		} else {
			$progStats{$curDB}{SearchCompl}++;
		}
		if ($jobN2 eq ""){ $jobN2 = $jobName; } else {$jobN2 .= ";".$jobName;}
	}
	return ($jobN2,$clnCmd);
}
sub nopareil(){
	my ($ar1,$outD,$Gdir,$name,$jobd) = @_;
	my $numCore = 4;
	my @re1 = @{$ar1};
	
	my $npBin = getProgPaths("nonpareil");#"/g/bork5/hildebra/bin/nonpareil/nonpareil";

	my $sumOut = "$name.npo";
	my $cmd = "mkdir -p $outD\n";
	$cmd .= "$npBin -s $re1[0] -f fastq -t 20 -m 40000 -b $outD$name -n 10240 -i 0.1 -m 0.2\n";#-t $numCore  -o $sumOut
	$cmd .= "cp $outD/$sumOut $Gdir";
	#R part
	#source('/g/bork5/hildebra/bin/nonpareil/utils/Nonpareil.R');
	#Nonpareil.curve('$outD/$sumOut');
	my $jobName = "_NP$JNUM"; my $tmpCmd;
	if (!-d $outD || !-e "$outD/$sumOut"){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($jobName,$tmpCmd) = qsubSystem($logDir."NonPar.sh",$cmd,1,"42G",$jobName,$jobd,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}
	return $jobName;
}
sub runContigStats{
	my ($path,$jobd,$assD,$subprts,$immSubm, $readL,$readLX, $tmpD,$AssemblyGo,$Nthr,$smpl) = @_;
	my $sepCtsScript = getProgPaths("sepCts_scr");#
	my $ContigStatsDir  = "$path/$dir_ContigStats";
	my $CSfilesComplete = 1;
	#die;
	$CSfilesComplete = 0 if (  (!-s "$ContigStatsDir/Coverage.count_pergene.gz" || !-e "$path/assemblies/metag/assembly.txt") );
	#print "CSfilesComplete $CSfilesComplete $subprts $ContigStatsDir\n!-s $ContigStatsDir/Coverage.count_pergene.gz || !-e $path/assemblies/metag/assembly.txt\n";
	$CSfilesComplete = 0  if ($MFopt{kmerPerGene} && $AssemblyGo && !-s "$assD/ContigStats/scaff.pergene.4kmer.pm5.gz" );
	$CSfilesComplete = 0  if ($MFopt{mapSupport2Assembly} && $map{$smpl}{"SupportReads"} ne "" &&  !-s "$ContigStatsDir/Cov.sup.count_pergene.gz" );
	$CSfilesComplete = 0  if ($subprts =~ m/F/ && !-s "$assD/ContigStats/FMG/FMGids.txt" );
	$CSfilesComplete = 0  if ($subprts =~ m/G/ && !-s "$assD/ContigStats/GTDBmg/marker_genes_meta.tsv" );
	#die "$CSfilesComplete";
	return ("","",0) if ($CSfilesComplete);
	
	print "Running Contig Stats on assembly ($immSubm)\n";
	$Nthr =1 unless ($subprts =~m/[FGEm]/);
	
	
	my $jobName = "_CS$JNUM"; $QSBoptHR->{LocationCheckStrg}=""; my $tmpCmd="";
	my $jobDep = "";
	#system "mkdir -p $cwd" unless ($cwd eq "" || -d $cwd);
	my $cmd = "";
	#$cmd .= "mkdir -p $cwd\n" unless ($cwd eq "" );
	$cmd .= "$sepCtsScript -inD $path -assD $assD -subparts $subprts -readLength $readL -readLengthSup $readLX -tmpD $tmpD -threads $Nthr";
	#$jobName = "_CS$JNUM"; 
	($jobDep,$tmpCmd) = qsubSystem($logDir."ContigStats.sh",$cmd,$Nthr,int(50/$Nthr)."G",$jobName,$jobd,"",$immSubm,[],$QSBoptHR);
	$tmpCmd = "" if ($immSubm);
	#die "$cmd\n$logDir.ContigStats.sh,$cmd,$Nthr,int(50/$Nthr).G,$jobName,$jobd,$cwd,$immSubm\n";
	return ($jobDep,$tmpCmd, 1);
}
sub calcCoverage{
	my ($cov,$gff,$RL,$cstNme,$jobd,$dirsHr) = @_;
	my $readCov_Bin =getProgPaths("readCov");
	my $jobName = ""; $QSBoptHR->{LocationCheckStrg}=""; my $tmpCmd="";
	my $cmd = "$readCov_Bin $cov $gff $RL";
	my $qdir = $logDir; $qdir = ${$dirsHr}{qsubDir} if (exists( ${$dirsHr}{qsubDir} ));
	my $ret = $cmd;
	if (!-s $cov.".pergene" || !-s $cov.".percontig" || !-s $cov.".median.percontig" ){
		$jobName = "_COV$JNUM";
		$jobName = "$cstNme"."_$JNUM" if ($cstNme ne "");
		#die "$cmd\n";
		if (${$dirsHr}{submit}){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			($jobName,$tmpCmd) = qsubSystem($qdir."COV$cstNme.sh",$cmd,1,"20G",$jobName,$jobd,"",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
		} 
	} else {
		$jobName = $jobd;
	}
	return ($jobName,$cmd);
}

sub checkDrives{
	my ($aref) = @_;
	my @locs = @{$aref};
	my $retStr = "\n####### BEGIN file location check ######\n";
	foreach my $llo (@locs){$retStr.="ls -l $llo > /dev/null \n";}
	$retStr.=" \nsleep 3\n"; my $cnt=0;
	foreach my $llo (@locs){ $cnt++;
		$retStr.="if [ ! -d \"$llo\" ]; then echo \'Location $cnt does not exist\'; exit 5; fi\n";
	}
	$retStr .= "####### END file location check ######\n";
	return $retStr;
}

# $sdmjN = cleanInput($cfp1ar,$cfp2ar,$sdmjN,$smplTmpDir);
 sub cleanInput($ $ $ $){
	my ($hr,$sdmjN,$saveD) = @_;
	my %seqSet = %{$hr};
	my @c1 = @{$seqSet{"pa1"}}; my @c2 = @{$seqSet{"pa2"}};
	my $cmd = "";
	for (my $i=0;$i<@c1;$i++){
		if ($c1[$i] =~ m/$saveD/){
			$cmd .= "rm -f $c1[$i] $c2[$i]\n" if (-e $c1[$i] || -e $c2[$i]);
		} else {
			#die "$c1[$i] =~ m/$saveD/\n";
		}
	}
	my @cs = @{$seqSet{"pas"}};
	for (my $i=0;$i<@cs;$i++){
		if ($c1[$i] =~ m/$saveD/){
			$cmd .= "rm -f $cs[$i] \n" if (-e $cs[$i] );
		} 
	}

	#die $cmd."\n";
	my $jobName = $sdmjN;
	if ($cmd ne ""){
		print "Removing raw input fastqs..\n";
		system $cmd;
		#$jobName = "_PC$JNUM"; my $tmpCmd;
		#my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		#($jobName, $tmpCmd) = qsubSystem($logDir."ClnUnzip.sh",$cmd,1,"1G",$jobName,$sdmjN,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		#$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}
	return $jobName;
 }
 
 
#Gap Filler to refine scaffolds
sub GapFillCtgs{
	my($ar1,$ar2,$scaffolds,$GFdir_a,$dep,$xtrTag) = @_; #.= "_GFI1";
	my $GFbin = getProgPaths("gapfiller");#"perl /g/bork5/hildebra/bin/GapFiller/GapFiller_n.pl";
	my @pa1 = @{$ar1}; my @pa2 = @{$ar2};
	my @inserts;
	my $prefi = "GF";
	my $numCore = 16;
	mkdir($GFdir_a.$prefi);
	my $log = $GFdir_a.$prefi."/GapFiller.log";
	#system("mkdir -p $GFdir_a$prefi");
	my $GFlib = ($GFdir_a."GFlib.opt");
	my @libFiles;
	for (my $i=0;$i<@pa1;$i++){
		push(@libFiles, ($pa1[$i].",".$pa2[$i]));
		push(@inserts,450);#default value for our hiSeq runs..
	}
	createGapFillopt($GFlib,\@libFiles,\@inserts);
#GF round 1
	my $cmd = "";
	$cmd .= "mkdir -p $GFdir_a$prefi\n";
	$cmd .= $GFbin . " -l $GFlib -s $scaffolds -m 75 -o 2 -r 0.7 -d 70 -t 10 -g 1 -T $numCore -b ".$prefi." -D $GFdir_a > $log \n";
	$cmd .= "rm -r $GFdir_a$prefi/alignoutput $GFdir_a$prefi/intermediate_results $GFdir_a$prefi/reads\n";
	my $finalGFfile = $GFdir_a."$prefi/$prefi.gapfilled.final.fa";
	#die $cmd."\n";
	if ($cmd ne ""){
		my $tmpCmd;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($dep, $tmpCmd) = qsubSystem($logDir."GapFill_ext_$xtrTag.sh",$cmd,$numCore,int(80/$numCore)."G","_GFE$JNUM",$dep,"",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}
	return $dep;
}

sub movePreAssmData{
	my ($metagD, $mvD,$mapD, $tmpD, $CSdir, $smplID ) = @_;
	#very thorough checks that everything is correctly prepped
	my $mvSTO = "$mvD/moved.sto";
#	system "rm -fr $metagD;\n" if ($AssemblyGo && -e $mvSTO);
	if ( -e $mvSTO){
		print "Moved assmbl already.. $mvSTO \n";
		return;
	}
	print "Preparing assembly from preassembly..";
	die "Coverage does not exist.. can't move preAssembly\n" if (!-e "$CSdir/Coverage.percontig.gz");
	die "Not a preassebmly?: $metagD\n" unless (-e "$metagD/$STOpreAssmblDone");
	die "Couldn't find ContigStats in $metagD\n" unless (-d $CSdir);
	die "Couldn't find Assembly $metagD/scaffolds.fasta.filt\n" unless (-e "$metagD/scaffolds.fasta.filt");
	my $cmd = "";
	#my $newCovFile = "$tmpD/$smplID.coverage.perCtg";
	$cmd .= "mkdir -p $tmpD;cp -r $metagD/scaffolds.fasta.filt $metagD/$STOpreAssmblDone $CSdir/Coverage* $tmpD;\n";
	#$cmd .= "cp $CSdir/Coverage.median.percontig $newCovFile\n";
	$cmd .= "rm -fr $mapD;\n" ;
	$cmd .= "cp -r $metagD/AssemblyStats.txt $logDir/preAssmStat.txt\n";
	#if ($AssemblyGo){$cmd .= "rm -fr $metagD;\n" ; print "removed preASsmbl dir";}
	$cmd .= "mkdir -p $mvD;\n";
	#$cmd .= "cp $tmpD/$STOpreAssmblDone $metagD;\n";
	$cmd .= "cp -r $tmpD/* $mvD/;\n";
	$cmd .= "rm -fr $tmpD\n";
	$cmd .= "touch $mvSTO\n";
	#die $cmd;
	systemW $cmd;
	print " Done \n";
	#return $newCovFile;
}

#used in hybrid assemblies
sub prepPreAssmbl{
	my ($metagD, $mvD,$mapD, $tmpD , $CSdir, $curSmpl, $cAssGrp, $assmDone, $ePreAssmblPck) = @_;
	#die "$mvD\n";
	my $doPreAssmFlag = 0;
	return ($doPreAssmFlag,0,$ePreAssmblPck) if ($assmDone);
	my $eCOV = 0; $eCOV = 1 if (-e "$CSdir/Coverage.percontig.gz");
	#print "$eCOV $CSdir/Coverage.percontig\n";
	if ($MFopt{DoAssembly} == 5 && $map{$curSmpl}{"SupportReads"} =~ m/PB:/ ){ #condition: right assembly mode and actually secondary support reads
		$doPreAssmFlag = 1 ;
		if (!$eCOV && !$ePreAssmblPck){
			#print "preAssmbl: nothing done yet.. \n";#$mvD\n$metagD";
			return ($doPreAssmFlag, 0, $ePreAssmblPck );
		}
	} else {
		return (0,0,0); 
	}
	#die "$mvD\n";
	#my $eCOV = 0; $eCOV =1 if ( -e "$CSdir/Coverage.percontig");
	$AsGrps{$cAssGrp}{CntPreAss} = 0 unless (exists($AsGrps{$cAssGrp}{CntPreAss}));
	
	
	if (($ePreAssmblPck || $eCOV) && -e "$metagD/$STOpreAssmblDone" ){
		#die "preAssmX: $PostAssemblyGo $doPreAssmFlag     $AsGrps{$cAssGrp}{CntPreAss} >= $AsGrps{$cAssGrp}{CntAimAss}\n";
		$doPreAssmFlag = 0 ;#no prep needed any longer.. files will/are saved already!
		#all ready for second assembly step!
		movePreAssmData($metagD, $mvD,$mapD, $tmpD , $CSdir, $curSmpl) ;
		$ePreAssmblPck = 1;
	}  
	#die "XAS\n";
	if ($ePreAssmblPck){
		$AsGrps{$cAssGrp}{CntPreAss} ++ ;
		$doPreAssmFlag = 0 ; #everyone else needs to keep
		push(@{$AsGrps{$cAssGrp}{preAsmblDir}}, $mvD);
	}
	my $PostAssemblyGo = 0;
	$PostAssemblyGo = 1 if (!$doPreAssmFlag && ($AsGrps{$cAssGrp}{CntPreAss} >= $AsGrps{$cAssGrp}{CntAimAss}) ); #has already seen enough complete preAssmblies
	$doPreAssmFlag = 1 if (!$PostAssemblyGo); 
	#print "-e $CSdir/Coverage.percontig   $metagD/$STOpreAssmblDone\n" ;
	#print "preAssm:  $doPreAssmFlag     $AsGrps{$cAssGrp}{CntPreAss} >= $AsGrps{$cAssGrp}{CntAimAss} :: $ePreAssmblPck $PostAssemblyGo\n";
	#die "$doPreAssmFlag\n";
	return ($doPreAssmFlag,$PostAssemblyGo,$ePreAssmblPck);
}

#scaffolding via mate pairs
sub scaffoldCtgs{
	my ($AsgHR,$ASG, $xar1,$xar2, $refCtgs, $tmpD1,$outD,$dep,$Ncore,$smplName,$spadesRef,$xtraTag) = @_;
	my ($ar1,$ar2,$ars,$liar,$rear) = getRawSeqsAssmGrp($AsgHR,$ASG,0);
	my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
	my $besstBin = getProgPaths("BESST");#"/g/bork3/home/hildebra/bin/BESST/./runBESST";
	my @pa1 = @{$ar1}; my @pa2 = @{$ar2}; my @libs = @{$liar};
	my @xpa1 = @{$xar1}; my @xpa2 = @{$xar2};
	return ("",$dep) if (@pa1 ==0 );
	my $tmpD = $tmpD1."/scaff/";
	my $bwtIdx = $refCtgs.$MFcontstants{bwt2IdxFileSuffix}; my $cmdDB="";my $chkFile = "";
	my $clnCmd = ""; my $spadesDir = ""; my $spadFakeDir = "";
	if ($spadesRef ){#&& -e ("$refCtgs/scaffolds.fasta")){#this is supposed to be the spades dir
		my $oldDir = "spades_ori";
		$spadFakeDir = $refCtgs;
		$clnCmd .= "mkdir -p $refCtgs/../$oldDir/; mv $refCtgs/* $refCtgs/../$oldDir/; mv $refCtgs/../$oldDir $refCtgs\n";
		$spadesDir = "$refCtgs$oldDir"; 
		$clnCmd .= "cp $spadesDir/smpls_used.txt $refCtgs\n";
		$refCtgs .= "$oldDir/scaffolds.fasta";
		$clnCmd .="\ngzip $spadesDir/*";
		$clnCmd .="\ngunzip $refCtgs";
		($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$refCtgs",$Ncore,0,$MFopt{MapperProg});#$Ncore);
	} elsif (!-e $bwtIdx){
		($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$refCtgs",$Ncore,0,$MFopt{MapperProg});
	}
	
	
	#my @rd1; my @rd2;
	my $algCmd = "$clnCmd\n$cmdDB\n";
	$algCmd .= "mkdir -p $tmpD\n";
	my @bams; my $cnt=0; my @insSiz; my @orientations;
	for (my $i=0;$i<@pa1;$i++){
		next unless ($libs[$i] =~ m/mate/i);
		#push @rd1,$rds[$i];push @rd2,$rds[$i+1];
		my $tmpOut = "$tmpD/tmpMateAlign$cnt.bam";
		my $tmpBAM = "$tmpD/tmpMateAlign$cnt.srt.bam";
		$algCmd .= "$bwt2Bin --no-unal --end-to-end -p $Ncore -x $bwtIdx -X $MFconfig{mateInsertLength} -1 $pa1[$i] -2 $pa2[$i] | $smtBin view -b -F 4 - > $tmpOut\n";
		$algCmd .= "$smtBin sort -@ $Ncore -T kk -O bam -o $tmpBAM $tmpOut; $smtBin index $tmpBAM\n";
		#$algCmd .= "$novosrtBin --ram 50G -o $tmpBAM -i $tmpOut \n";
		$algCmd .= "rm $tmpOut\n\n";
		push(@bams,$tmpBAM); push(@insSiz,10000);push(@orientations,"fr");
		$cnt++;
		#print "sc  mat\n";
	}
	for (my $i=0;$i<@xpa1;$i++){#assumes paired end reads
		#push @rd1,$rds[$i];push @rd2,$rds[$i+1];
		my $tmpOut = "$tmpD/tmpMateAlign$cnt.bam";
		my $tmpBAM = "$tmpD/tmpMateAlign$cnt.srt.bam";
		$algCmd .= "$bwt2Bin --no-unal --end-to-end -p $Ncore -x $bwtIdx  -1 $xpa1[$i] -2 $xpa2[$i] | $smtBin view -b -F 4 - > $tmpOut\n";
		$algCmd .= "$smtBin sort -@ $Ncore -T kk -O bam -o $tmpBAM $tmpOut; $smtBin index $tmpBAM\n";
		#$algCmd .= "$novosrtBin --ram 50G -o $tmpBAM -i $tmpOut \n";
		$algCmd .= "rm $tmpOut\n\n";
		push(@bams,$tmpBAM);push(@insSiz,500);push(@orientations,"fr");
		$cnt++;
		#print "sc  mat\n";
	}
	
	#die "\n\n\n$algCmd\n\n\n";
	return ("",$dep) if (@bams ==0 );
	#create bowtie2 mapping

	my $zcmd = "-z 5000";  #-z 10000";
	 my $ori = "--orientation ".join(" ",@orientations);
	#for (my $i=1;$i<@bams;$i++){ $ori .= " fr"}#$zcmd .= " 10000";
	my $cmd = "";
	my $bams = join(" ",@bams);
	system "mkdir -p $outD" unless (-d $outD);
	$cmd .= "mkdir -p $outD\n";
	$cmd .= "$besstBin $zcmd $ori -f $bams -o $outD -c $refCtgs\n"; #-q $Ncore <- unstable?
	my $jobName = $dep;
	$cmd .= "rm -r $tmpD\n";
	#cleanup2, fake spades result folder
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	if ($spadesDir ne ""){
		my $assStatScr = getProgPaths("assStat_scr");#"perl assemblathon_stats.pl";
		my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
		$cmd .= "mv $outD/BESST_output/pass1/Scaffolds_pass1.fa $spadFakeDir/scaffolds.fasta\n";
		$cmd .= "$renameCtgScr $spadFakeDir/scaffolds.fasta $smplName\n";
		$cmd .= "$sizFiltScr $spadFakeDir/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
		$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $spadFakeDir/scaffolds.fasta > $spadFakeDir/AssemblyStats.500.txt\n";
		$cmd .= "$assStatScr $spadFakeDir/scaffolds.fasta > $spadFakeDir/AssemblyStats.ini.txt\n";
		$cmd .= "$assStatScr $spadFakeDir/scaffolds.fasta.filt > $spadFakeDir/AssemblyStats.txt\n\n";
		my ($cmdX,$bwtIdxX,$chkFileX) = buildMapperIdx("$spadFakeDir/scaffolds.fasta.filt",$Ncore,0,$MFopt{MapperProg});
		$cmd .= $cmdX."\n";

	}
	my $newScaffFNA = "$outD/BESST_output/pass2/Scaffolds_pass2.fa";
	$cmd .= "$renameCtgScr $newScaffFNA $smplName\n";
	$cmd .= "\ntouch $outD/scaffDone.sto\n" ;
	$cmd = "" if (-e "$outD/scaffDone.sto");
	
#die "scaff cmd $cmd\n";
	if ($cmd ne ""){
		my $tmpCmd;
		($dep, $tmpCmd) = qsubSystem($logDir."BesstScaff_ext_$xtraTag.sh",$algCmd.$cmd,$Ncore,int(50/$Ncore)."G","_BBE$JNUM$xtraTag",$dep,"",1,[],$QSBoptHR);
	}
	return ($newScaffFNA,$dep);
}

#preprocess mate pairs (use nxtrim on them and communicate results)
 sub check_mates($ $ $ $ $){
	my ($ar,$ifastasPre,$mateD,$doMateCln,$dep) = @_;
	my @mat = @{$ar};
	my $cmd = "";
	my @sarPre; my $mateC=0;
	my @mates;
	my $nxtrimBin = getProgPaths("nxtrim");#"/g/bork3/home/hildebra/bin/NxTrim/./nxtrim";

	foreach my $matp (@mat){
		system "mkdir -p $mateD" unless (-d $mateD);
		my @mateX = split /,/,$matp;
		push(@sarPre,"$mateD/mate${mateC}.se.fastq.gz");
		push(@mates,"$mateD/mate.${mateC}_R1.unknown.fastq.gz","$mateD/mate.${mateC}_R2.unknown.fastq.gz","$mateD/mate.${mateC}_R1.mp.fastq.gz","$mateD/mate.${mateC}_R2.mp.fastq.gz");
		next if ( -e "$mateD/matesDone.sto");
		#--rf keeps reads in rf; --joinreads joins pe
		$cmd .= "$nxtrimBin --ignorePF --separate -1 $mateX[0] -2 $mateX[1] -O $mateD/mate.$mateC\n";
		#$nxtrimBin --stdout-mp -1 $rd1 -2 $rd2 | $bwaBin mem $refCtgs -p - > out.sam
		$ifastasPre .= ";$mateD/mate.${mateC}_R1.pe.fastq.gz,$mateD/mate.${mateC}_R1.pe.fastq.gz";
		$mateC++;
	}
	#die "@mates SDS\n";
	$cmd .= "touch $mateD/matesDone.sto\n";
	#die $cmd;
	$cmd = "" if ($doMateCln == 0 || $mateC==0);
	my $jobName = $dep;
	if ($cmd ne ""){
		my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."mateClean.sh",$cmd,1,"1G","_NX$JNUM",$dep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
	}
	return ($ifastasPre,\@sarPre,\@mates,$jobName) ;
 }

 
 #preprocess mate pairs (use nxtrim on them and communicate results)
 sub check_matesL($ $ $ $){
	my ($pa1,$pa2,$mateD,$doMateCln) = @_;
	my %ret;
	my $cmd = "";
	my $nxtrimBin = getProgPaths("nxtrim");#"/g/bork3/home/hildebra/bin/NxTrim/./nxtrim";

	$cmd .= "rm -rf $mateD; mkdir -p $mateD\n";
	 my $mateC=0;
	#foreach my $matp (@mat){
	if ($mateC > 0){die "mate pairs only supports single library\n";}
	system "mkdir -p $mateD" unless (-d $mateD);
	#my @mateX = split /,/,$matp;
	#--rf keeps reads in rf; --joinreads joins pe
	$cmd .= "$nxtrimBin --ignorePF --separate -1 $pa1 -2 $pa1 -O $mateD/mate.$mateC\n";
	$cmd .= "rm -f $pa1 $pa2\n";
	#$nxtrimBin --stdout-mp -1 $rd1 -2 $rd2 | $bwaBin mem $refCtgs -p - > out.sam
	$ret{pe1} = "$mateD/mate.${mateC}_R1.pe.fastq.gz"; $ret{pe2} = "$mateD/mate.${mateC}_R1.pe.fastq.gz";
	$ret{se} =  "$mateD/mate${mateC}.se.fastq.gz";
	$ret{un1} = "$mateD/mate.${mateC}_R1.unknown.fastq.gz"; $ret{un2} = "$mateD/mate.${mateC}_R2.unknown.fastq.gz";
	$ret{mp1} = "$mateD/mate.${mateC}_R1.mp.fastq.gz"; $ret{mp2} = "$mateD/mate.${mateC}_R2.mp.fastq.gz";
	$mateC++;
	#}
	#die "@mates SDS\n";
	my $locStone = "$mateD/matesDone.sto";
	$cmd .= "touch $locStone\n";
	#die $cmd;
	$cmd = "" if ($doMateCln == 0 || $mateC==0);
	
	return (\%ret,$cmd,$locStone) ;
 }

 
 #determines what read types are present and starts cleaning of mates, if required
 sub get_ifa_mifa($ $ $ $ $ $ $){
	my ($ifastas,$ifastasS,$libInfoAr, $mateD, $doMateCln, $jdep, $singlReadMode) = @_;
	my @allFastas = split /;/,$ifastas;
	my @libInfo = @{$libInfoAr}; my @miSeqFastas = (); my $mcnt=0;
	my @mates;
	#die "@libInfo\n";
	foreach (@libInfo){
		if ($_ =~ m/.*miseq.*/i) {
			push (@miSeqFastas, $allFastas[$mcnt]);
			splice(@allFastas, $mcnt, 1);
		} elsif ($_ =~ m/.*mate.*/i){ #remove from process
			push (@mates, $allFastas[$mcnt]);
			splice(@allFastas,$mcnt,1);
		}
		$mcnt++;
	}
	$ifastas = join(";",@allFastas);
	my $mi_ifastas = join(";",@miSeqFastas);
	my $singleIfas = $ifastasS; my $matRef = [];
	#die "$ifastas\n$singleIfas\n\n";
	#moved to seedUnzip2tmp
	#($ifastas,$singleAddAr,$matRef,$jdep) = check_mates(\@mates,$ifastas,$mateD,$doMateCln,$jdep);

	return ($ifastas, $mi_ifastas, $singleIfas, $matRef, $jdep);
 }
 sub get_sdm_outf($ $ $ $ $ $){
	my ($ifastas, $mi_ifastas,$finD,$singlReadMode,$pairedReadMode, $useXtras) = @_;
	my @ret1;my @ret2;my @sret;
	my $fEnd = "fq";
 	if ($MFopt{gzipSDMOut}){
		$fEnd = "fq.gz";
	}
	my $baseFname = "filtered";
	$baseFname = "filtered.suppl" if ($useXtras);
	
	if ( $pairedReadMode){
		push(@ret1,$finD."$baseFname.1.$fEnd");push( @ret2, ($finD."$baseFname.2.$fEnd"));push(@sret, ($finD."$baseFname.singl.$fEnd"));
	} 
	if ($singlReadMode && !$pairedReadMode){ #only single reads avaialble, different file ending..
		push(@sret, $finD."$baseFname.s.$fEnd");
	}
	if ($mi_ifastas ne ""){
		push(@ret1,$finD."${baseFname}_mi.1.$fEnd");push( @ret2, ($finD."${baseFname}_mi.2.$fEnd"));push(@sret, ($finD."${baseFname}_mi.singl.$fEnd"));
	}

	return (\@ret1,\@ret2,\@sret);
}
 sub check_sdm_loc(){
	die "defunct functon  check_sdm_loc!\n";
	my ($ar1,$ar2,$libInfoAr,$sdmO) = @_;
 	my $ifastasPre = ${$ar1}[0].",".${$ar2}[0];
	for (my $i=1;$i<@{$ar1};$i++){if ($i>0){$ifastasPre .= ";".${$ar1}[$i].",".${$ar2}[$i];}}
	my ($ifastas, $mi_ifastas, $singlAddAr, $matAr, $jdep) = get_ifa_mifa($ifastasPre,"",$libInfoAr,"",0,"","???");
	my ($ar1x,$ar2x,$sar) = get_sdm_outf($ifastas, $mi_ifastas,$sdmO,"???","",0);
	my @ret1=@{$ar1x};my @ret2=@{$ar2x};my @sret=@{$sar};
	my $presence=1;
	foreach my $loc (@ret2,@ret1){	if (!-e $loc || -z $loc){$presence=0;}	}
	#for (my $i=0;$i<;$i++){	if ( !-e $ret2[0]|| -z $ret1[0]){$presence=0;}	}
	if (!-e $sret[0]){$presence=0;}
	my $stone = "$sdmO/filterDone.stone";
	if (!-e $stone){$presence=0;}
	if (!$presence){return 0;}
	my $assInputFlaw = 0;
	if (-e $logDir."spaderun.sh.otxt"){open I,"<$logDir/spaderun.sh.otxt" or die "Can't open old assembly logfile $logDir\n"; my $str = join("", <I>); close I;
		if ($str =~ /Assertion `seq_\.size\(\) == qual_\.size\(\)' failed\./){$assInputFlaw=1;}	
		if ($str =~ /paired_readers\.hpp.*Unequal number of read-pairs detected in the following files:/){$assInputFlaw=1;}	
	}

	if ($presence==0 || $assInputFlaw==1){
		return 0;
	} else {
		return 1;
	}
}

#($mergRdsHsh,$mergJbN) = mergeReads($arp1,$arp2,$sdmjN,$smplTmpDir."merge_clean/");
sub mergeReads(){
	my ($cleanSeqSetHR,$jdep,$outdir,$doMerge,$runThis) = @_;
	my $arp1 = ${$cleanSeqSetHR}{arp1}; my $arp2 = ${$cleanSeqSetHR}{arp2}; 

	my $flashBin = getProgPaths("flash");
	my $numCores = 8;
	my $outT = "sdmCln";
	my %ret = (mrg => "", pair1 => "", pair2 => "");
	if (!exists(${$cleanSeqSetHR}{mrgHshHR} )) {${$cleanSeqSetHR}{mrgHshHR} = {};}
	#print "XSADS\n!$doMerge || !$MFconfig{readsRpairs} || !$runThis\n";
	if (!$doMerge || !$MFconfig{readsRpairs} || @{$arp2} == 0 || !$runThis){return ($cleanSeqSetHR,"");}
	#print "XSADS1\n";
	#if (@{$arp1} > 1 ){die "Array with reads provided to merging routine is too large!\n@{$arp1}\n";}
	my $mergCmd  = "";
	for (my $i=0; $i<@{$arp1};$i++){
		my $outTL = $outT;
		$outTL .= ".$i" if ($i > 0);
		$mergCmd .= "$flashBin -M 250 -z -o $outT -d $outdir -t $numCores ${$arp1}[0] ${$arp2}[0]\n";
		if ($i > 0){
			$mergCmd .= "cat $outTL.extendedFrags.fastq.gz >> $outT.extendedFrags.fastq.gz;cat $outTL.notCombined_2.fastq.gz >> $outT.notCombined_2.fastq.gz; ";
			$mergCmd .= "cat $outTL.notCombined_1.fastq.gz >> $outT.notCombined_1.fastq.gz;\n";
		}
	}
	my $stone = "$outdir/$outT.sto";
	$mergCmd .="touch $stone\n ";
	my $jobName = "";
	if (-e $stone && -e "$outdir/$outT.extendedFrags.fastq" && !-e "$outdir/$outT.extendedFrags.fastq.gz"){
		#zip
		$mergCmd = "$pigzBin -f -p $numCores $outdir/$outT.extendedFrags.fastq $outdir/$outT.notCombined_1.fastq $outdir/$outT.notCombined_2.fastq\n";
		system "rm $stone";
	}
	if (!-e $stone){
		$jobName = "_FL$JNUM"; my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."flashMrg.sh",$mergCmd,$numCores,"3G",$jobName,$jdep,"",1,[],$QSBoptHR);
	}
	
	#die "$logDir/flashMrg.sh\n$mergCmd\n";
	$ret{mrg} = "$outdir/$outT.extendedFrags.fastq.gz";
	$ret{pair1} = "$outdir/$outT.notCombined_1.fastq.gz";
	$ret{pair2} = "$outdir/$outT.notCombined_2.fastq.gz";
	
	${$cleanSeqSetHR}{mrgHshHR} = \%ret;
	
	
	return ($cleanSeqSetHR,$jobName);
}
 
sub adaptSDMopt{
#adaptSDMopt($baseSDMopt,$globalLogDir,$samplReadLength);
	my ($baseSF,$oDir,$RL,$RT) = @_;
	my $newSDMf = $oDir."/sdmo_${RL}_${RT}.txt";
	my $nRL = $RL - 10 - int($RL/15);
	open I,"<$baseSF" or die "Can't open sdm opt in:\n $baseSF\n"; my $str = join("", <I>); ;close I;
	#print "\n\n$RT\n$baseSF\n";
	if ($RT ne "proto" && $RT ne "PB" && $RT ne "ONT" && $RL != 0){
		#$str =~ s/minSeqLength\t\d+/minSeqLength\t$nRL/;
		if ($RL < 50){$str =~ s/maxAccumulatedError\t.*\n/maxAccumulatedError\t0.5\n/;
		} elsif ($RL < 90){$str =~ s/maxAccumulatedError\t\d+\.\d+/maxAccumulatedError\t1.2/; #really shitty GAII platform
		} elsif ($RL < 200){$str =~ s/maxAccumulatedError\t\d+\.\d+/maxAccumulatedError\t2.5/;} #really shitty GAII platform
		if ($RL < 90){$str =~ s/TrimWindowWidth	\d+/TrimWindowWidth	8/;$str =~ s/TrimWindowThreshhold	\d+/TrimWindowThreshhold	16/;	
		$str =~ s/maxAmbiguousNT	\d+/maxAmbiguousNT	1/;
		} else {$str =~ s/TrimWindowWidth	\d+/TrimWindowWidth	18/;$str =~ s/TrimWindowThreshhold	\d+/TrimWindowThreshhold	20/;
		}
		$str =~ s/maxAmbiguousNT	\d+/maxAmbiguousNT	2/;
	}
	if ($MFopt{sdmProbabilisticFilter} ==0 ){
		$str =~ s/BinErrorModelAlpha\t.*\n/BinErrorModelAlpha\t-1\n/;
	}
	foreach my $so (keys %{$MFopt{sdm_opt}}){$str =~ s/$so	[^\n]*\n/$so	$MFopt{sdm_opt}->{$so}\n/;}
	#die $str."\n$baseSF\n";
	open O,">$newSDMf" or die "Can't open new sdm opt out:\n $newSDMf\n"; print O $str; close O;
	#print $newSDMf."\n";
	return $newSDMf;
}

sub getReadTechInMap{
	my $curReadTec = "hiSeq"; #default	
	if (exists($map{$curSmpl}{SeqTech})){ $curReadTec = $map{$curSmpl}{SeqTech};}
	return $curReadTec;
}

sub getReadTechInMapSingl($){
	my ($curReadTec )= @_;
	if (exists($map{$curSmpl}{SeqTechSingl})){ $curReadTec = $map{$curSmpl}{SeqTechSingl};}
	return $curReadTec;
}


sub sdmOptSet{ 
	my ($curSmpl,$samplReadLength, $curReadTec, $curSTech) = @_;
	
	if ($MFopt{sdmOpt} ne ""){
		if (! -f $MFopt{sdmOpt}){die "-customSDMopt must point to file! (currently: $MFopt{sdmOpt})\n";}
		if (! -s $MFopt{sdmOpt}){die "-customSDMopt must point to non-empty file! (currently: $MFopt{sdmOpt})\n";}
		
		return ($MFopt{sdmOpt},$MFopt{sdmOpt});
	}
	my $curSDMopt = $baseSDMopt; 
	#my $iqualOff = 33; #62 for 1st illu
	
	#print "XXXXXXXXXXX $curReadTec XXXXXXXXXX\n";
		
	if ($curReadTec eq "GAII_solexa" || $curReadTec eq "GAII"|| $curReadTec eq "hiSeq"){
		#$iqualOff = 59; #really that old??
	} elsif ($curReadTec eq "miSeq"){ $curSDMopt = $baseSDMoptMiSeq; 
	} elsif ($curReadTec eq "proto"){ $curSDMopt = getProgPaths("baseSDMoptProto"); 
	} elsif ($curReadTec eq "PB"){ $curSDMopt = getProgPaths("baseSDMoptPacBio"); 
	} elsif ($curReadTec eq "ONT"){ $curSDMopt = getProgPaths("baseSDMoptONT"); 
	}
	
	if (!exists($sampleSDMs{$curReadTec}{$samplReadLength})){
		#die "wrong";
		$curSDMopt = adaptSDMopt($curSDMopt,$globalLogDir,$samplReadLength,$curReadTec);
	}
	my $curSDMoptSingl = $baseSDMopt;
	#$curSDMopt = $sampleSDMs{$curReadTec}{$samplReadLength};
	
	#singletons..
	if ($curSTech eq ""){$curSDMoptSingl=$curSDMopt;
	} elsif ($curSTech eq "454"){$curSDMoptSingl = getProgPaths("baseSDMopt454"); 
	} elsif ($curSTech eq "miSeq"){ $curSDMoptSingl = $baseSDMoptMiSeq; 
	} elsif ($curSTech eq "proto"){ $curSDMoptSingl = getProgPaths("baseSDMoptProto"); 
	} elsif ($curSTech eq "PB"){ $curSDMoptSingl = getProgPaths("baseSDMoptPacBio"); 
	} elsif ($curReadTec eq "ONT"){ $curSDMopt = getProgPaths("baseSDMoptONT"); 
	}
	
	if ($curSTech eq "PB" && $samplReadLength < 1000){
		print "WARNING: it seems sequencing technology is PacBio (\"PB\"), but read length is very short: $samplReadLength\nConsider adjusting via \"-inputReadLength\" or  \"-inputReadLengthSuppl\"\n";
	}
	if ($curSTech eq "ONT" && $samplReadLength < 1000){
		print "WARNING: it seems sequencing technology is Oxford Nanopore (\"ONT\"), but read length is very short: $samplReadLength\nConsider adjusting via \"-inputReadLength\" or  \"-inputReadLengthSuppl\"\n";
	}
	#print "$curSDMopt,$curSDMoptSingl\n";
	return ($curSDMopt,$curSDMoptSingl);
}
sub sdmClean(){
	my ($curOutDir,$seqSetHR,$cleanSeqSetHR,$mateD,$finD,$jobd,$runThis, $useXtras ) = @_;
	#create ifasta path
	#die "cllll\n";
	my %seqSet = %{$seqSetHR};
	#sdm options..
	#my $sdmO = $seqSet{"curSDMopt"}; my $sdmS = $seqSet{"curSDMoptSingl"};
	
	my ($ar1,$ar2,$ars,$libInfoAr,$seqTec) = ($seqSet{"pa1"},$seqSet{"pa2"},$seqSet{"pas"},$seqSet{"libInfo"},$seqSet{"seqTech"});
	my $samplReadLength = $seqSet{"samplReadLength"};
	#my ($ar1X,$ar2X,$arsX,$libInfoArX) = ($seqSet{"paX1"},$seqSet{"paX2"},$seqSet{"paXs"},$seqSet{"libInfoX"})
	if ($useXtras){
		#if useXtras flag set, only consider "xtra" defined input reads (eg scaffoling reads, 3rd gen in supplement of main reads etc)
		($ar1,$ar2,$ars,$libInfoAr,$seqTec) = ($seqSet{"paX1"},$seqSet{"paX2"},$seqSet{"paXs"},$seqSet{"libInfoX"},$seqSet{"seqTechX"});
		$samplReadLength = $seqSet{"samplReadLengthX"};
		if (!@{$ar1} && !@{$ars}) { #nope, no additional reads are requested..
			#return [],[],[],[],"";
			return ($cleanSeqSetHR,$jobd);
		}
		print "cleaning support reads..\n";
	}
	my $sdmBin = getProgPaths("sdm");# sdm program from LotuS2 pipeline
	my $comprCores=$MFopt{sdmCores};
	if (@{$ar1} == 0 && @{$ars}==0){die "Empty array to sdmClean given\n";}
	my $ifastasPre="";my $ifastasPreS=""; my $singlReadMode=0; my $pairReadMode=0;
	if (@{$ar1} != 0){
		$ifastasPre = ${$ar1}[0].",".${$ar2}[0]; $pairReadMode=1;
		for (my $i=1;$i<@{$ar1};$i++){if ($i>0){$ifastasPre .= ";".${$ar1}[$i].",".${$ar2}[$i];}}
		#die "${$ar1}[0]\n";
	} 
	if (@{$ars} > 0)	{
		$ifastasPreS = ${$ars}[0];
		for (my $i=1;$i<@{$ars};$i++){if ($i>0){$ifastasPreS .= ";".${$ars}[$i];}}
		$singlReadMode=1;
		
	}
	
	my $curReadTec = getReadTechInMap();
	my $curSTech = getReadTechInMapSingl($curReadTec);
	if ($seqTec ne ""){#override defaults from map
		$curReadTec = $seqTec; $curSTech = $seqTec;
	} 
	checkSeqTech($curReadTec, "MG-TK.pl::sdmClean");
	my $is3rdGen = is3rdGenSeqTech($curReadTec);

	
	my ($sdmO,$sdmS) = sdmOptSet($curSmpl,$samplReadLength, $curReadTec, $curSTech );
	#print "curReadTec :: $curReadTec $is3rdGen\n$sdmO $sdmS\n";
	
	if (!$singlReadMode && !$pairReadMode){
		die "Couldn't find any valid input files (sdmClean)\n";
	}
	#die "$pairReadMode\n$singlReadMode\n";
	my $stone = "$finD/filterDone.stone";
	if (!-e "$curOutDir/input_fil.txt" ){
		open O,">$curOutDir/input_fil.txt"; print O $ifastasPre; close O;
	}

	#system("mkdir -p $finD");
	my $cmd = "";
	#$cmd .= "mkdir -p $tmpD\n" unless ($tmpD eq "");
	
	my ($ifastas, $mi_ifastas, $singlIfas, $matAr, $jdep) = get_ifa_mifa($ifastasPre,$ifastasPreS,$libInfoAr,$mateD,0,$jobd,$singlReadMode);
	$jobd = $jdep; #replaces old dep
	#die "$ifastas\n";
	my ($ar1x,$ar2x,$sar) = get_sdm_outf($ifastas, $mi_ifastas,$finD,$singlReadMode,$pairReadMode, $useXtras);
	
	#push(@{$sar},@{$singlAddAr});
	my @ret1=@{$ar1x};my @ret2=@{$ar2x};my @sret=@{$sar};
	$cmd .= "mkdir -p $finD\n" unless ($finD eq "");
	$cmd .= "rm -f $stone  $sret[0]\n"; #sret needs to be removed, as I concat later onto it..
	$cmd .= "rm -f $sret[1]\n" if (@sret > 1);
	$cmd .= "mkdir -p $logDir/sdm/\n";
	
	
	#die "$mi_ifastas\n$ifastas XX\n@ret1\n@sret\n";
	#$baseSDMoptMiSeq;
	my $ofiles = "";
	my $mi_ofiles = "";
	#system("mkdir -p $logDir/sdm/") unless (-d "$logDir/sdm/");
	my $useLocalTmp = 1;
	#if ($tmpD eq ""){$useLocalTmp = 0;$tmpD = $finD;}
	#my $nsdmPair = 2;
	my $sdm_def= " -ignore_IO_errors 1 -i_qual_offset auto -binomialFilterBothPairs 1 -threads $MFopt{sdmCores} ";
	$sdm_def .= " -XfirstReads $MFconfig{XfirstReads} " if ($MFconfig{XfirstReads} > 0);
	$sdm_def .= " -illuminaClip 1 " if ($MFopt{trimAdapters} && $curReadTec ne "PB" && $curReadTec ne "ONT" );
	my $sdm_Extr = "";
	$sdm_Extr .= "-logLvsQ 1 " if ($MFopt{SDMlogQualvsLen});
	# $ret{$curSmp}{cut5pR2}
	my $sdm_cut = "";
	$sdm_cut .= "-5PR1cut $map{$curSmpl}{cut5pR1} " if ($map{$curSmpl}{cut5pR1} > 0); $sdm_cut .= "-5PR2cut $map{$curSmpl}{cut5pR2} " if ($map{$curSmpl}{cut5pR2} > 0);
	$sdm_cut .= "-XfirstReadsRead $map{$curSmpl}{firstXrdsRd} " if ($map{$curSmpl}{firstXrdsRd} > 0);
	$sdm_cut .= "-XfirstReadsWritten $map{$curSmpl}{firstXrdsWr} " if ($map{$curSmpl}{firstXrdsWr} > 0);
	
	my $catCmd = "cat ";
	#$catCmd = "$pigzBin -p $comprCores -c " if ($MFopt{gzipSDMOut});
	if ($pairReadMode){
		$ofiles = $ret1[0].",".$ret2[0];# $finD."/filtered.1.fq,".$finD."/filtered.2.fq";
		$cmd .= "$sdmBin -i \"$ifastas\" -o_fastq $ofiles -options $sdmO -paired 2 $sdm_Extr -log $logDir/sdm/filter.log $sdm_def $sdm_cut\n\n";
		my $tmpTar = "$ret1[0]";
		if ($MFopt{gzipSDMOut}){$tmpTar =~ s/\.\d\.fq\.gz/\.\*\.singl\.fq\.gz/g ;
		} else {$tmpTar =~ s/\.\d\.fq/\.\*\.singl\.fq/g ;}
		#ls /path/to/your/files* 1> /dev/null 2>&1; 
		my $cmdSA = "if ls $tmpTar 1> /dev/null 2>&1; then $catCmd $tmpTar ";
		$cmd .= $cmdSA . ">>  $sret[0]; \n" ;
		$cmd .= "rm -f $tmpTar\nfi\n";
	}
	if ($singlReadMode){
		my $tmpO = $sret[0];$tmpO =~ s/(.*)\//$1\/tmp\./;
		#quick fix for sdm problems with single fq out
		if (1){	$catCmd = "$pigzBin -f -p $comprCores -c "; $tmpO =~ s/\.fq\.gz/\.fq/ ;}
		#die $tmpO."\n";
		$cmd .= "$sdmBin -i \"$singlIfas\" -o_fastq $tmpO -options $sdmS $sdm_Extr -paired 1  -log $logDir/sdm/filterS.log $sdm_def $sdm_cut\n\n";
		if ($pairReadMode){
			$cmd .= "$catCmd $tmpO >>  $sret[0]\nrm $tmpO\n" ;
		} else {
			if ($tmpO =~ m/\.gz$/){
			$cmd .= "mv $tmpO $sret[0]\n";
			} else { 
				$cmd .= "$catCmd $tmpO >>  $sret[0]\nrm $tmpO\n" ;
			}
		}		
		
	}
	if ($mi_ifastas ne ""){ #miSeq specific filtering
		$mi_ofiles =  $ret1[1].",".$ret2[1];
		$cmd .= " $sdmBin -i \"$mi_ifastas\" -o_fastq $mi_ofiles $sdm_Extr -options $baseSDMoptMiSeq -paired 2  -log $logDir/sdm/filter_xtra.log $sdm_def $sdm_cut\n\n";
		if (!$singlReadMode){
			my $tmpTar = "$ret1[1]";
			if ($MFopt{gzipSDMOut}){$tmpTar =~ s/\.\d\.fq\.gz/\.\*\.singl\.fq\.gz/g ;
			} else {$tmpTar =~ s/\.\d\.fq/\.\*\.singl\.fq/g ;}
			my $cmdSA = "if ls $tmpTar 1> /dev/null 2>&1; then $catCmd $tmpTar ";
			$cmd .= $cmdSA . ">>  $sret[1]; rm -f $tmpTar; fi\n" ;
		}
	}
	#die $cmd."\n";
#	if ($useLocalTmp){
#		if ($MFconfig{unpackZip} || $MFopt{gzipSDMOut}){#zip output
#			foreach my $of (split(/,/,$ofiles)){
#				$of =~ m/\/([^\/]+$)/; my $fname = $1;
#				$cmd .= "$pigzBin -p $comprCores -c $of > $finD/$fname.gz\n";
#			}
#			if (!$singlReadMode){
#				$sret[-1] =~ m/\/([^\/]+$)/; my $fname = $1;
#				#$cmd .= "$pigzBin -p $comprCores -c ".$sret[-1]." > $finD/$fname.gz\n" ;
#			}
#		} else {
#			$cmd .= "mv -f ".join(" ",(split(/,/,$ofiles),split(/,/,$mi_ofiles)))." $finD\n";
#			$cmd .= "rm -rf $tmpD\n";
#		}
#	} elsif ($MFconfig{unpackZip} || $MFopt{gzipSDMOut}){#zip output
#		foreach my $of (split(/,/,$ofiles)){
#			$of =~ m/\/([^\/]+$)/; my $fname = $1;
#			$cmd .= "$pigzBin -p $comprCores $of \n";
#		}
#		if (!$singlReadMode){
#			$sret[-1] =~ m/\/([^\/]+$)/; my $fname = $1;
#			#$cmd .= "$pigzBin -p $comprCores ".$sret[-1]." \n" ;
#		}
#	}
	if (!$singlReadMode){
		$sret[-1] =~ m/\/([^\/]+$)/; my $fname = $1;
		#$cmd .= "$pigzBin -p $comprCores ".$sret[-1]." \n" ;
	}
	$cmd .= "touch $stone\n";
	#die "$cmd\n";
	my $jobName = "";
	my $presence=1; my $gzPres=1;
	if (!$MFconfig{unpackZip}){
		foreach my $loc (@ret2,@ret1){	if (!-e $loc || -z $loc){$presence=0;} my $loc2 = $loc;$loc2 =~ s/\.gz$//;	if (!-e $loc2){$gzPres=0;}}
		if (!-e $sret[0]){$presence=0;}
	}
	if (!-e $stone){$presence=0;}
	
	#recovery part, that uses previous sdm filtered data..
	if (!$presence && $gzPres && -e $stone){
		#exists, but not gzipped..
		$cmd = "";
		foreach my $loc (@ret2,@ret1,@sret){
			my $loc2 = $loc;
				$loc2 =~ s/\.gz$//;
			$cmd .= "$pigzBin -f -p $comprCores $loc2\n";
		}
	}
	#die "$presence\n";
	
	#common error: spade input failed
	my $assInputFlaw = 0;
	my $qsubFile = $logDir."sdmReadCleaner.sh";
	if ($useXtras){
		${$cleanSeqSetHR}{arpX1} = \@ret1;${$cleanSeqSetHR}{arpX2} = \@ret2;
		${$cleanSeqSetHR}{singArX} = \@sret;${$cleanSeqSetHR}{matArX} = $matAr;
		${$cleanSeqSetHR}{readTecX} = $curReadTec; ${$cleanSeqSetHR}{is3rdGenX} = $is3rdGen;
		$qsubFile = $logDir."sdmReadCleanerSuppl.sh";
	} else {#default..
		${$cleanSeqSetHR}{arp1} = \@ret1;${$cleanSeqSetHR}{arp2} = \@ret2;
		${$cleanSeqSetHR}{singAr} = \@sret;${$cleanSeqSetHR}{matAr} = $matAr;
		${$cleanSeqSetHR}{readTec} = $curReadTec;  ${$cleanSeqSetHR}{is3rdGen} = $is3rdGen;
	}
				#print "@{${$cleanSeqSetHR}{arp1}} YY \n";

	if ( ($presence==0 || $assInputFlaw==1 )&& $runThis){
		#die "yes\n";
		$jobName = "_SDM${useXtras}_$JNUM"; my $tmpCmd;
		my $preHDDspace=$QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = 0;
		($jobName, $tmpCmd) = qsubSystem($qsubFile,$cmd,$MFopt{sdmCores},"15G",$jobName,$jobd.";".$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		$QSBoptHR->{tmpSpace} = $preHDDspace;
	}
	#die "$presence presi\n@ret1\n";
	#print $cmd;
	#return (\@ret1,\@ret2,\@sret,$matAr,$jobName);
	#$cleanSeqSet{arp1},$cleanSeqSet{arp2},$cleanSeqSet{singAr},$cleanSeqSet{matAr}
	return ($cleanSeqSetHR,$jobName);
}
sub mocat_reorder(){
	my($ar1,$ar2,$singlAr,$inJob) = @_;
	
	#my @pairs = split(";",$ifastas);
	my @ret1 = @{$ar1}; my @ret2 = @{$ar2}; 
	#print "@ret1\n@ret2\n";
	#foreach (@pairs){		my @spl = split /,/;		push(@ret1,$spl[0]);push(@ret2,$spl[1]);	}
	#my $jobName = $inJob;
	return (\@ret1,\@ret2,$singlAr,$inJob);
}
sub SEEECER(){
	my ($p1ar,$p2ar,$tmpD) = @_;
	die("SEECER deactive\n");
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	if ( $p1[0] =~ m/.*\.gz/){
		system("gunzip -c ".join(" ",@p1)." > $tmpD/pair.1.fastq");	system("gunzip -c ".join(" ",@p2)." > $tmpD/pair.2.fastq");
		@p1 = ("$tmpD/pair.1.fastq"); @p2 = ("$tmpD/pair.2.fastq");
	} elsif ( @p1 > 1 ){
		system("cat ".join(" ",@p1)." > $tmpD/pair.1.fastq");	system("cat ".join(" ",@p2)." > $tmpD/pair.2.fastq");
		@p1 = ("$tmpD/pair.1.fastq"); @p2 = ("$tmpD/pair.2.fastq");
	}
	my $SEEbin = "bash /g/bork5/hildebra/bin/SEECER-0.1.3/SEECER/bin/run_seecer.sh";
	system("mkdir -p $tmpD/tmpS");
	my $cmd = $SEEbin . " -t $tmpD/tmpS $p1[0] $p2[0]";
	#qsubSystem($logDir."SEECERCleaner.sh",$cmd,1,"30G",1);
	my @ret1= ($tmpD."pair.1.fastq_corrected.fa");my @ret2= ($tmpD."pair.2.fastq_corrected.fa");
	return (\@ret1,\@ret2);
}


#calculates md5 sums needed to upload filtered raw fastq's to EBI/SRA
sub unploadRawFilePostprocess{
	#die "XXn\n";
	if (@EBIjobs == 0){return;}
	if ($MFconfig{uploadRawRds} eq ""){return ;}
	if (-e "$MFconfig{uploadRawRds}/R1.md5"){return;}
	print "Postprocess upload postprocess\n";
	#md5sum --version
	my $cmd = "md5sum $MFconfig{uploadRawRds}/*.R2.fq.gz > $MFconfig{uploadRawRds}/R2.md5\n";
	$cmd .= "md5sum $MFconfig{uploadRawRds}/*.R1.fq.gz > $MFconfig{uploadRawRds}/R1.md5\n";
	my ($jobN, $tmpCmd) = qsubSystem("$globalLogDir/postEBI.sh",$cmd,1,"20G","_PP",join(";",@EBIjobs),"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
}

#cleans raw fastqs fastq's (human DNA) to prepare upload to EBI/SRA
sub uploadRawFilePrep{
	if ($MFconfig{uploadRawRds} eq ""){return "" ;}
	my ($inD,$tmpD,$cfp1ar,$cfp2ar,$cfpsar,$smplID,$libInfoAr,$totalInputSizeMB) = @_;
	my @pa1 = @{$cfp1ar}; my @pa2 = @{$cfp2ar}; my @sa = @{$cfpsar};
	my @libInfo = @{$libInfoAr};

	my $unsplBin = getProgPaths("unsplitKrak_scr"); #unsplitting merged reads..
	my $krk2Bin = getProgPaths("kraken2");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
	#flag whether to use kraken v1 or v2
	my $krkBin = "";
	my $krak1= 1;	if ($krk2Bin ne ""){ $krak1=0; }
	if ($krak1){
		$krkBin = getProgPaths("kraken");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
	}

	#prepare databases
	my $DBdir = $krakenDBDirGlobal."/";
	my @DBname = ("hum1stTry"); 
	if (@filterHostDB > 0){
		@DBname = (); my $cnt=0;
		foreach my $fhdb (@filterHostDB){
			if (!-d $fhdb){die "-filterHostKrakDB is not a dir: \"$fhdb\"\n";}
			$DBname[$cnt] = $fhdb; $DBdir=""; $cnt++;
		}
	}


	my $outD  ="$MFconfig{uploadRawRds}"; 
	my $numThr = 4;
	system "mkdir -p $outD/tmp/ " unless (-d "$outD/tmp");
	my $rd;my @rds;
	my $cmd = "";#"rm -rf $outD/tmp/;mkdir -p $outD/tmp/\n";
	for (my $i=0;$i<3;$i++){
		next if ($i==1); #read2 will be dealt with read1
		if ($i==0){$rd="1";@rds=@pa1;}
		#if ($i==1){$rd="2";@rds=@pa2;}
		my @rds2= @pa2;
		if ($i==2){$rd="single";@rds=@sa;}
		my $idx=0;
		foreach (my $x=0;$x<@rds;$x++){
			my $f =$rds[$x];
			my $rd2=$rds2[$x];
			my $xtra="";
			if ($libInfo[$idx] =~ m/.*mate.*/i){$xtra = "mate.";}
			if ($libInfo[$idx] =~ m/.*miseq.*/i){$xtra = "miSeq.";}
			my $of = "$tmpD/$smplID.$xtra$idx.R$rd.fq";  
			my $of2 = "$tmpD/$smplID.$xtra$idx.R2.fq";  
			my $tmpF = "$tmpD/$smplID.$xtra$idx.R1R2.fq";  
			my $ff = "$outD/$smplID.$xtra$idx.R$rd.fq";  
			my $ff2 = "$outD/$smplID.$xtra$idx.R2.fq";  
			my $ofT = "$tmpD/tmp/$smplID.$xtra$idx.TEMP.R$rd.fq";
			my $ofT2 = "$tmpD/tmp/$smplID.$xtra$idx.TEMP.R2.fq";
			if ($f =~ m/\.gz$/){$of .= ".gz";$ff .= ".gz";$ofT .= ".gz";$of2 .= ".gz";$ff2 .= ".gz";$ofT2 .= ".gz";}
			$idx++;
			#print "$ff\n";
			next if (-e $ff);
			$cmd .= "rm -fr $tmpD;\nmkdir -p $tmpD/tmp/ $outD\n";
			#$cmd .= "rm -f $ofT;cp $inD$f $ofT\n" unless (-e $of);
			$cmd .= "rm -f $ofT;ln -s $inD$f $ofT\n" unless (-e $of);
			if ($i==2){
				$cmd .= krakHSapSingl("$ofT",$of,$numThr)."\n" unless (-e $of.".gz");
				$ofT2="";$of2="";
			}
			if ($i==0){
				my $gzFlag = "";$gzFlag =  "--gzip-compressed" if ($ofT =~ m/\.gz$/);
				my $gzEnd = ""; $gzEnd = ".gz" if ($gzFlag ne "");

				$cmd .= "rm -f $ofT2;cp $inD$rd2 $ofT2\n" unless (-e $of);
				if ($krak1){
					$cmd .= "$krkBin --preload --threads $numThr --paired --fastq-input $gzFlag --unclassified-out $tmpF --db $DBdir$DBname[0]  $ofT $ofT2 > /dev/null\n";
					#overwrites input files
					$of =~ s/\.gz$//;$of2 =~ s/\.gz$//;
					$cmd .= "$unsplBin $tmpF $of $of2\nrm -f $tmpF*\n";
					if ($gzFlag ne ""){
						$cmd .= "$pigzBin -f -p $numThr $of $of2\n";
						$of .= ".gz";$of2 .= ".gz";
					}
				} else {
					
					$tmpF ="$tmpD/krak.tmp#.fq"; my $tmpF1 = "$tmpD/krak.tmp_1.fq";my $tmpF2 = "$tmpD/krak.tmp_2.fq";
					$cmd .= "$krk2Bin --threads $numThr $gzFlag --paired --unclassified-out $tmpF --db $DBdir$DBname[0] --output - $MFopt{filterHostKr2QuickMode} --confidence $MFopt{krakHostConf} $ofT $ofT2 \n";
					$cmd .= "$pigzBin -f -p $numThr $tmpF1 $tmpF2\n" unless ($gzFlag eq "");
					$cmd .= "mv $tmpF1$gzEnd $of; mv $tmpF2$gzEnd $of2\n";
				}
				
				#die "$cmd";
			}
			$cmd .= "rm -f $ofT $ofT2\n";
			$cmd .= "mv $of $of2 $outD\n";
		}
	}
	$cmd .= "rm -rf $tmpD\n" if ($cmd ne "");
	
	#$cmd = "rm -rf /local/hildebra/MF/\n";
	
	#die "$cmd\n";
	my $retJob = "";
	if ($cmd ne ""){
		#systemW $cmd;
		my $preHDDspace=$QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = int($totalInputSizeMB/1024*6)+30  ."G";
#		$QSBoptHR->{tmpSpace} = $HDDspace{prepPub}; #increase local space..  # use $totalInputSizeMB ??
		#print "$QSBoptHR->{tmpSpace}\n";
		my ($jobN, $tmpCmd) = qsubSystem("$logDir/prepEBI.sh",$cmd,$numThr,int(50/$numThr) . "G","_PP$JNUM",$krakDeps,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} = $preHDDspace;
		$retJob = $jobN;
	}
	return $retJob;
}
#function to copy files into mocat compatible format .. due to mocat inflexibility for input *sic*
sub mocatFileCpy($$$$$){
	if ($MFconfig{mocatLinkDir} eq ""){return 0 ;}
	my ($inD,$cfp1ar,$cfp2ar,$cfpsar,$smplID) = @_;
	my @pa1 = @{$cfp1ar}; my @pa2 = @{$cfp2ar}; my @sa = @{$cfpsar};
	my $outD  ="$MFconfig{mocatLinkDir}/$smplID"; 
	my $cmd = "";
	$cmd = "mkdir -p $outD\n" unless (-d $outD);
	my $idx=0;
	foreach my $f (@pa1){
		my $of = "$outD/raw.$idx.1.fq"; $of .= ".gz" if ($f =~ m/\.gz$/);
		$cmd .= "ln -s $inD$f $of\n" unless (-e $of);
		$idx++;
	}
	$idx=0;
	foreach my $f (@pa2){
		my $of = "$outD/raw.$idx.2.fq"; $of .= ".gz" if ($f =~ m/\.gz$/);
		
		$cmd .= "ln -s $inD$f $of\n"  unless (-e $of);
		$idx++;
	}
	$idx=0;
	foreach my $f (@sa){
		my $of  = "$outD/raw.$idx.single.fq"; $of .= ".gz" if ($f =~ m/\.gz$/);
		$cmd .= "ln -s $inD$f $of\n"  unless (-e $of);
		$idx++;
	}
	
	#die "$cmd\n";
	if ($cmd ne ""){
		systemW $cmd;
		#my ($jobN, $tmpCmd) = qsubSystem("$logDir/cp2mocat.sh",$cmd,1,"1G","_CPM$JNUM",$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
	}
	return 1;
}



sub complexGunzCpMv($ $ $ $ $ $){
	my ($fastap,$in,$scrathD,$finDest,$ncore, $allowLinks) = @_;
	my $unzipcmd = "";
	my $lowEffort=1;
	
	my $out = $in;
	if (0&&$in =~ m/\.gz$/){ #deactivated.. takes too much space on file sys
		$unzipcmd .= "$pigzBin  -f -p $ncore -d -c $fastap$in";
		$out =~ s/\.gz$//;
		$unzipcmd .= " > $finDest$in \n";
		$lowEffort=0;
	}elsif ($in =~ m/\.bz2/){
		my $bzip2B  = getProgPaths("bzip2");
		
		$out =~ s/\.bz2$/\.gz/;
		$unzipcmd .= "$bzip2B --decompress --stdout $fastap$in | $pigzBin  -f -p $ncore -c -f > $finDest$out";
		$lowEffort=0;
	}elsif ($allowLinks){
		$unzipcmd .= "ln -s $fastap$in $finDest/$out;\n";
		
	} else {
		$unzipcmd .= "cp $fastap$in $finDest;\n";
		$lowEffort=0;
	}
	if (0 && $scrathD ne $finDest){
		$unzipcmd .= "mv -f $scrathD$in $finDest\n";
		$lowEffort=0;
		die "$scrathD ne $finDest";
	}
	my $newRDf = $finDest."$out";
	
	return ($unzipcmd,$newRDf,$lowEffort);
}

sub outfiles_trimall($ $){
	my ($opath,$fil) = @_;
	my $OFp1 = "$opath/$fil"; 
	my $OFu1 = "$opath/$fil";	
	if ($OFu1 =~ m/\.gz$/){
		$OFu1 =~ s/(\.[^\.]+\.gz)$/\.sing$1/ ;#$OFp1 =~ s/\.gz$// 
	} elsif ($OFu1 =~ m/\.f[^\.]*q$/){
		$OFu1 =~ s/(\.f[^\.]*q)$/\.sing$1/ ;
	} else {
		die "Unknown file ending: $fil\n";
	}
	#die "$OFp1\n$OFu1\n";
	return ($OFp1,$OFu1);
}
sub outfiles_Bam($ $){
	my ($opath,$fil) = @_;
	my $OFu1 = "$opath/$fil"; 
	if ($OFu1 =~ m/\.bam$/){
		$OFu1 =~ s/(\.bam)$/\.unbam\.fq\.gz/ ;
	} else {
		die "Unknown file ending: $fil\n";
	}
	#die "$fil\n$OFu1\n";
	return $OFu1;
}


sub valid_files{
	my ($path, $paIR) = @_;
	my @paI = @{$paIR};
	my @paO;
	foreach my $ff (@paI){
		if (-e "$path/$ff" && !-d "$path/$ff"){
			push (@paO, $ff) ;
			print "$path/$ff\n";
		}
	}
	return \@paO;
}

sub seedUnzip2tmp{
	my ($fastp,$curSmpl,$jDepe,$tmpPath,$finDest, $WT,
		$libNum,$calcUnzp,$finalMapDir,$porechopFlag,$inputRawFile) = @_;
	my $mocatImport = $MFconfig{importMocat};
	my $himipeSeqAd = getProgPaths("illuminaTS3pe"); #for trimomatic
	my $himiseSeqAd = getProgPaths("illuminaTS3se");
	my $trimJar = getProgPaths("trimomatic");

	my $smplPrefix = $map{$curSmpl}{prefix};
	$smplPrefix =~ s/\./\\\./g;   #make sure dots are recognized as such

	my $xtrMapStr = $map{$curSmpl}{SupportReads};
	my $doMateCln = 1; #nxtrim  mate pairs ?
	my $numCore=$MFopt{unzipCores};
	my $rawReads=""; my $mmpu = "";
	my @fastp2 = (""); my $xtraRdsTech = "";
	$xtrMapStr = "" if (!defined $xtrMapStr) ;
	my $totalInputSizeMB=10000; #default to something sensible

	if ($fastp eq "") { print "No primary dir.. \n"; }
	#main source for input files..
	my @pa1 = (); my @pa2 = (); my @pas = (); 
	#support reads.. 
	my @paX1= (); my @paXs= (); my @paX2= ();
	#information on libraries (name of library)
	my @libInfo= (); my @libInfoX= ();
	#information on seq tech
	my $seqTech = ""; my $seqTechX = "";
	#store bam formated input..
	my @paBam; my @paBamX;
	
		#die "Mapping $pa1 to ref\n";
	if ( $xtrMapStr ne ""){
		#die "XX$xtrMapStr\n";
		my @spl = split(/:/,   $map{$samples[$JNUM]}{"SupportReads"}   );
		if (@spl > 1){
			$xtraRdsTech = $spl[0];#.$libNum;
			@fastp2 = split(/,/,$spl[1]);
			checkSeqTech($xtraRdsTech,"MATFILER.pl::Support Reads");
			$seqTechX = $xtraRdsTech;
		}
		#die $fastp2."\n".$xtraRdsTech."\n";
	}
	my $is3rdGen = is3rdGenSeqTech($seqTech);
	my $is3rdGenX = is3rdGenSeqTech($seqTechX);
	
	
	#create empty return object
	my %seqSet = (pa1 => \@pa1, pa2 => \@pa2, pas => \@pas, seqTech => $seqTech, is3rdGen => $is3rdGen,
			paX1 => \@paX1, paX2 => \@paX2, paXs => \@paXs, seqTechX => $seqTechX, is3rdGenX =>  $is3rdGenX,
			libInfo => \@libInfo, libInfoX => \@libInfoX,
			totalInputSizeMB => -1,
			rawReads => $rawReads,
			mmpu => $mmpu, WT => $WT,
			);

	
	#check if unmapped reads requested..
	if ($MFopt{useUnmapped}){
		$seqSet{pa1} = ("$finalMapDir/unaligned/unal.1.fq.gz");
		$seqSet{pa2} = ("$finalMapDir/unaligned/unal.2.fq.gz");
		$seqSet{pas} = ("$finalMapDir/unaligned/unal.fq.gz");
		$seqSet{totalInputSizeMB} = filsizeMB((@pa1,@pa2,@pas));
		$seqSet{libInfo} = ("unmapped");
#		return ("",\@pa1,\@pa2, \@pas, 0, "", "",\@libInfo, "",$totalInputSizeMB);
		return ("", \%seqSet);

	}
	

	if ($fastp ne ""){
		#die "here\n";
		if ( !-d $fastp ){
			my $msg = "Infile dir not existing: $fastp\n";
			die $msg if ($MFconfig{abortOnEmptyInput});
			print $msg; $seqSet{totalInputSizeMB}=0;
			#return ("EMPTY_DO_NEXT",\@pa1,\@pa2, \@pas, 0, "", "",\@libInfo, "",$totalInputSizeMB);
			return ("EMPTY_DO_NEXT", \%seqSet);

		}
		#print "$fastp\n@pa1\n$smplPrefix\n";print "MOCA is $mocatImport $MFconfig{readsRpairs}\n";
		#if ($mocatImport==0){
		#print "$fastp\n";
		#die "$smplPrefix$rawFileSrchStr2\n";
		opendir(DIR, $fastp) or die "Infile dir not existing: $fastp\n";	
		my @DIRf = readdir(DIR) ;
		close(DIR);
		if ($MFconfig{rawFileSrchStrSingl} ne ""){
			@pas = sort ( grep { /$smplPrefix($MFconfig{rawFileSrchStrSingl})/  && -e "$fastp/$_"} @DIRf);	
		}
		if ($MFconfig{readsRpairs}!=0){
			@pa2 = sort ( grep { /$smplPrefix($MFconfig{rawFileSrchStr2})/ && -e "$fastp/$_" } @DIRf );	
			@pa1 = sort ( grep { /$smplPrefix($MFconfig{rawFileSrchStr1})/  && -e "$fastp/$_"  } @DIRf );	#rewinddir(DIR);
			#die "readsRpairs::@pa1\n$smplPrefix$rawFileSrchStr1\n";
		}
		
		#sort out multi assignments of read files (grep not uniquely defined)
		if ($MFconfig{rawFileSrchStrSingl} ne ""){
			my %h;
			my $rdPairsBf=scalar @pa1;
			if ($MFconfig{prefSinglFQgreps}){
				@h{(@pas)} = undef;
				@pa1 = grep {not exists $h{$_}} @pa1;
				@pa2 = grep {not exists $h{$_}} @pa2;
			} else {
				@h{(@pa1,@pa2)} = undef;
				@pas = grep {not exists $h{$_}} @pas;
			}
			my $rdPairsAft=scalar @pa1;
			if ($rdPairsAft != $rdPairsBf){
				print "corrected from $rdPairsBf read pairs to $rdPairsAft read pairs from input dir\n";
				print "s: @pas\n1: @pa1\n2: @pa2\n";
			}
		}
		
		if ($MFconfig{rawFileBamSrchSing} ne "" ){
			@paBam = sort ( grep { /$smplPrefix($MFconfig{rawFileBamSrchSing})/  && -e "$fastp/$_"  } @DIRf );	
		}
		
	}
	#die "1:@pa1\nS:@pas\n$MFconfig{rawFileSrchStrSingl}\n$MFconfig{readsRpairs}\n\n";
	
	#import xtra long rds
	if (@pa1 != @pa2 && $MFconfig{readsRpairs}!=0){
		print "For dir $fastp, unequal fastq files exist:\nP1\n".join("\n",@pa1)."\nP2:\n".join("\n",@pa2)."\n";@pa1=(); @pa2=();
	}
	$totalInputSizeMB = filsizeMB($fastp,(@pa1,@pa2,@pas,@paBam));
	$inputFileSizeMB{$curSmpl} = $totalInputSizeMB;
	#die "@pa1\n$fastp\n$totalInputSizeMB\n";
	#check for date of files (special filter)
	if ($MFconfig{doDateFileCheck}){
		my @pa1t = @pa1; my @pa2t = @pa2; undef @pa1; undef @pa2;
		for (my $i=0; $i<@pa1t;$i++){
			my $date = POSIX::strftime( "%m", localtime( ( stat "$fastp/$pa1t[$i]" )[9] ) );
			#print $date;
			if ($date < 3){push(@pa1, $pa1t[$i]);push(@pa2, $pa2t[$i]);}
		}
		if (@pa1 != 1){die "still too many file: $fastp\n @pa1\n@pa2\n";}
	}
	
	#create libinfo array 
	@libInfo = ("lib$libNum") x int @pa1;
	if (@pa1==0 && @pas!=0 || @pas > @pa1 #case of just single read in metag
			 ){ #or actually always..
		@libInfo = ("libS$libNum") x int @pas;
	} elsif(@pas > 0){
		@libInfo = (@libInfo, ("libS$libNum") x int @pas);
	} elsif(@paBam > 0){
		@libInfo = (@libInfo, ("libS$libNum") x int @paBam);
	}
	$locStats{hasPaired} = 1 if (@pa1 > 0);
	$locStats{hasSingle} = 1 if (@pas > 0);
	if ($fastp2[0] ne ""){
		#die "$fastp2\n\n";
		if (-d $fastp2[0]){
			opendir(DIR, $fastp2[0]) or die "Fasta indir not found $fastp2[0]\n";	
			my @DIRf = readdir(DIR) ;
			close(DIR);
			@paX2 = sort ( grep { /$MFconfig{rawFileSrchStrXtra2}/ && -e "$fastp2[0]/$_" } @DIRf );	#rewinddir DIR;
			@paX1 = sort ( grep { /$MFconfig{rawFileSrchStrXtra1}/  && -e "$fastp2[0]/$_"} @DIRf );	#close(DIR);
			if (@paX2 != @paX1){die "For dir $fastp, unequal fastq files exist:P1\n".join("\n",@paX1)."\nP2:\n".join("\n",@paX2)."\n";}
			if (@paX2 == 0){die "Can't find files with pattern $MFconfig{rawFileSrchStrXtra2} in $fastp2[0]\n";}
			#die $paX1[0]."\n";
			push @pa1 , @paX1; push @pa2 , @paX2; @libInfoX = $xtraRdsTech x @paX1; push @libInfo, @libInfoX;
		} elsif (-e $fastp2[0] ){
			if ($fastp2[0]=~ m/,/){die "MG-TK.pl::support reads have comma: $fastp2[0]!\nExiting..\n";}
			if ($fastp2[0] !~ m/\.bam$/ && $fastp2[0] !~ m/\.fastq\.gz?$/ && $fastp2[0] !~ m/\.fq\.gz$/){die "MG-TK.pl::support reads not bam or fq.gz formatted: $fastp2[0]!\nExiting..\n";}
			push @libInfoX, $xtraRdsTech;
			if ($fastp2[0] =~ m/\.bam$/){
				push @paBamX, @fastp2;
			} else {
				push @paXs, @fastp2;
			}
		} else {
			die "MG-TK.pl:: Can't find support reads @fastp2!\nExiting..\n";
		}
	}
		#die "@libInfo\n$libInfo[1]\n";

	#die "@pa1";
	#die "For dir $fastp, unual fastq files exist:P1\n".join("\n",@pa1)."\nP2:\n".join("\n",@pa2)."\n";
	
	#check if raw file is a symlink and if this is valid & create raw read link (HD times):
	for (my $i = 0; $i<@pa1; $i++){
		my $pp = $fastp;
		#$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
		if (-l $pp.$pa1[$i]){
			if (!-f abs_path($pp.$pa1[$i])){die "File $pa1[$i] is not file.\n";}
		}
		if (-l $pp.$pa2[$i]){
			if (!-f abs_path($pp.$pa2[$i])){die "File $pa2[$i] does not exist.\n";}
		}
		if ($i==0){
			$rawReads="$pp$pa1[$i],$pp$pa2[$i]";
		} else {
			$rawReads.=";".$pp.$pa1[$i].",".$pp.$pa2[$i];
		}
#		print "$pp$pa1[$i]\n";
		my $realP = `readlink -f $pp$pa1[$i]`;
		if (defined $realP && $realP =~ m/\/(MMPU[^\/]+)\//){
			$mmpu = $1;
		}
	}
	#die $rawReads."\n";
	#DEBUG fix to reduce file sizes
	#@fastap2 = ($fastap2[0]); @fastap1 = ($fastap1[0]);
	
	#could not find any files.. report this to user
	if (@pa1 == 0 && @pas ==0 && @paBam == 0){
		my $msg = "Can;t find files in $fastp\nUsing search pattern: $smplPrefix$MFconfig{rawFileSrchStr1}  $smplPrefix$MFconfig{rawFileSrchStr2}\n$smplPrefix$MFconfig{rawFileSrchStrSingl}\n";
		die $msg if ($MFconfig{abortOnEmptyInput});
		print $msg;$totalInputSizeMB=0;$inputFileSizeMB{$curSmpl} =0;
		#return ("EMPTY_DO_NEXT",\@pa1,\@pa2, \@pas, 0, "", "",\@libInfo, "",$totalInputSizeMB);
		$seqSet{pa1} = \@pa1;$seqSet{pa2} = \@pa2;$seqSet{pas} = \@pas;
		$seqSet{libInfo} = \@libInfo;
		return ("EMPTY_DO_NEXT", \%seqSet);

	}
	
	#set up stone to mark end of run
	my $finishStone = "$finDest/rawRds/done.sto";
	my $trimoStone = "";
	my $porechStone = "";
	if ($porechopFlag){
		$porechStone = "$finDest/rawRds/poreChopped.stone";
	}
	if (!$porechopFlag ){
		$MFconfig{filterFromSource} = 0;#too dangerous..
	}
	
	#report on what was found so far..
	if ($inputFileSizeMB{$curSmpl} > 0){
		print "Input size raw (Mb): " . int($inputFileSizeMB{$curSmpl}) ;
		print " Fastq pairs: " . scalar(@pa1) if (@pa1);
		print " Fastq Singls: " . scalar (@pas) if (@pas);
		print " Fastq Supports Singls: " . scalar (@paXs) if (@paXs);
		print " Bam Singls: " . scalar (@paBam) if (@paBam);
		print " Bam Supports Singls: " . scalar (@paBamX) if (@paBamX);
		
		
		print "\n" ;
	}

	
	
	#relinking in mocat file structure, if requested ## not used any longer
	my $mocatFCDone = mocatFileCpy($fastp,\@pa1,\@pa2,\@pas,$curSmpl);
	
	#prepare files to be uploaded to EBI etc, if requested
	my $uplDone = uploadRawFilePrep($fastp,$tmpPath,\@pa1,\@pa2,\@pas,$curSmpl,\@libInfo,$totalInputSizeMB);
	push (@EBIjobs, $uplDone);

	if ($mocatFCDone || $uplDone ne ""){ $totalInputSizeMB=0;$inputFileSizeMB{$curSmpl}=0;
		#return ("EMPTY_DO_NEXT",\@pa1,\@pa2, \@pas, 0, "", "",\@libInfo, "",$totalInputSizeMB);
		$seqSet{pa1} = \@pa1;$seqSet{pa2} = \@pa2;$seqSet{pas} = \@pas;
		$seqSet{libInfo} = \@libInfo;
		return ("EMPTY_DO_NEXT", \%seqSet);

	}
	die "tmpPath empty: $tmpPath" if ($tmpPath eq "");
	$tmpPath.="/rawRds/";
	my $unzipcmd = "";
	#$unzipcmd .= "set -e\n"; #sleep $WT\n
	$unzipcmd .= "rm -rf $finishStone $trimoStone $porechStone $finDest/rawRds/\nmkdir -p $finDest/rawRds/;\nsleep 1;\n";
	my $unzipcmdTMP = "rm -r -f $tmpPath;\nmkdir -p $tmpPath;\n";
	#make sure input is unzipped <- deprecated, in newer MF versions input is .gz
	my $testf1 = "";my $testf2 = "";
	my $lowEffort =-1; my $allowLinks=0;
	$allowLinks = 1 if (@pa1 == 1);
	$allowLinks = 1 if (@pa1 == 0 && @pas == 1 && !$porechopFlag);
	my $illCLip = $himipeSeqAd;
	if ($map{$curSmpl}{clip} ne ""){$illCLip = $map{$curSmpl}{clip};}
	#die "Can't find illumina trimming file: $illCLip\n" if ($useTrimomatic && !-e $illCLip);
	
	# check if only stone exists, but not the files any longer..
	if (-d "$finDest/rawRds"){
		opendir(my $dh, "$finDest/rawRds") or die "Could not open $finDest/rawRds for reading: $!\n";
		my @files = grep { -f "$finDest/rawRds/$_" } readdir($dh);
		closedir($dh);
		if (scalar(@files) == 1 && -e $finishStone) {
			system "rm -rf $finDest/rawRds" ;
		}
	}

	
	
	#first conversion of bam to fastqs::
	#this is currently only working with unpaired reads!
	for (my $i=0; $i<@paBam; $i++){
		my $pp = $fastp;
		my $smtBin = getProgPaths("samtools");
		$pas[$i] = outfiles_Bam("$finDest/rawRds/",$paBam[$i]);
		$unzipcmd .= "\necho \"Converting bam $i to fastq\"\n";
		$unzipcmd .= "$smtBin fastq -@ $numCore -t $pp/$paBam[$i] -0 $pas[$i];\n"; #| $pigzBin -p $numCore -c >
		$lowEffort = 0;
	}
	#and also take care of support reads in bam format
	for (my $i=0; $i<@paBamX; $i++){
		my $pp = $fastp;
		my $smtBin = getProgPaths("samtools");
		my $BamF = $paBamX[$i]; $BamF =~ s/.*\//\//;
		my $supportDir = "$finDest/rawRds/Support/";
		system "mkdir -p $supportDir" if ($i==0 && !-d $supportDir);
		$paXs[$i] = outfiles_Bam($supportDir,$BamF);
		$unzipcmd .= "echo \"Converting support bam $i to fastq\"\n";
		$unzipcmd .= "mkdir -p $supportDir;\n";
		$unzipcmd .= "$smtBin fastq -@ $numCore -t $paBamX[$i] -0 $paXs[$i];\n"; # | $pigzBin -p $numCore -c > 
		$lowEffort = 0;
	}
	
	
	
	
	
	for (my $i=0; $i<@pa1; $i++){
		#print $pa1[$i]."\n";
		my $pp = $fastp."/";
		if (@paBam){
			#do nothing, all done in paBam step
			;
		} elsif ($MFconfig{filterFromSource}){
			$lowEffort =1 if ($lowEffort != 0);
			$pa1[$i] = $pp.$pa1[$i]; $pa2[$i] = $pp.$pa2[$i];
		} elsif ($porechopFlag){
#			$unzipcmd .= "$porechBin \n";
			die "porechop is not implemented for read pairs!\n";
		} elsif (0){#$useTrimomatic){
			die "trimomatic no longer supported.. function is now implemented in sdm!\n";
			#$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
			#trimomatic instead of unzip
			my ($OFp1,$OFu1) = outfiles_trimall("$finDest/rawRds/",$pa1[$i]);
			my ($OFp2,$OFu2) = outfiles_trimall("$finDest/rawRds/",$pa2[$i]);
			$unzipcmd .= "java -jar $trimJar PE -threads $numCore $pp/$pa1[$i] $pp/$pa2[$i] $OFp1 $OFu1 $OFp2 $OFu2 ILLUMINACLIP:$illCLip:2:30:10\n";
			#for now: discard of singletons
			$unzipcmd .= "rm -f $OFu1 $OFu2\n";
			$pa1[$i] = $OFp1; $pa2[$i] = $OFp2;
			$lowEffort=0;
		} else {
		#old style
			my ($tmpCmd,$newF,$LEloc) = complexGunzCpMv($pp,$pa1[$i],$tmpPath,$finDest."/rawRds/",$numCore,$allowLinks);
			$unzipcmd .= $tmpCmd."\n";
			$pa1[$i] = $newF;
			$lowEffort = 0 if ($LEloc==0);
			($tmpCmd,$newF,$LEloc) = complexGunzCpMv($pp,$pa2[$i],$tmpPath,$finDest."/rawRds/",$numCore,$allowLinks);
			$unzipcmd .= $tmpCmd."\n";
			$pa2[$i] = $newF;
			
			$lowEffort = 0 if ($LEloc==0);
			#$lowEffort = 1 if ($allowLinks && $lowEffort == -1);
		}
	}
	#$unzipcmd .= "chmod +w ".join(" ",@pa1)."\n" if (!$MFconfig{filterFromSource} && @pa1>0);
	#$unzipcmd .= "chmod +w ".join(" ",@pa2)."\n" if (!$MFconfig{filterFromSource} && @pa2>0);
	#die "$unzipcmd\n";
	#die "@pa1\n";
	#for porechop this might be tons of files..
	
	for (my $i=0; $i<@pas; $i++){
		next if (scalar(@paBam));
		my $porechopped = "$finDest/rawRds/$pas[$i]"; $porechopped .= ".gz" unless ($porechopped =~ m/\.gz$/);
		my $pp = $fastp;
		if ($i==0 && ($porechopFlag && $is3rdGen) && !$allowLinks){$unzipcmd .=  "\nrm -f $porechopped\ntouch $porechopped\n\n";}
		#print "$libInfo[$i] eq $xtraRdsTech\n";
		#$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
		if ($MFconfig{filterFromSource}){
			$pas[$i] = $pp.$pas[$i]; 
		} elsif ($porechopFlag){
			#porechop is running really slow and instable, probably better to get fast5 and use modern basecaller, that will do this automatically..
			my $porechBin = getProgPaths("porechop");
			$unzipcmd .= "$porechBin -i $pp/$pas[$i] -t $numCore  --adapter_threshold 90 |gzip -c >> $porechopped\n";
			if (@pa1 > 0 ){die "no paired end reads can be given together with porechopped long reads!\n";}
		} elsif ($is3rdGen){
			if ($pas[$i] =~ m/\.gz$/){
				if (@pas > 1 || !$allowLinks){
					$unzipcmd .= "cat $pp/$pas[$i] >> $porechopped\n";
				} else {
					$unzipcmd .= "ln -s  $pp/$pas[$i] $porechopped\n";
					$lowEffort = 1 if ($lowEffort != 0);
				}
			} else {
				$unzipcmd .= "$pigzBin -f -p $numCore -c $pp/$pas[$i] >> $porechopped\n";
			}
		} else {
			my ($tmpCmd,$newF) = complexGunzCpMv($pp,$pas[$i],$tmpPath,$finDest."/rawRds/",$numCore,$allowLinks);
			$unzipcmd .= $tmpCmd."\n";
			$pas[$i] = $newF;
			if ($MFconfig{splitFastaInput} != 0){ #in case of input assemblies (MG-RAST.. arghh!!)
				my $sizSplitScr = getProgPaths("sizSplit_scr");
				$unzipcmd .= "\n$sizSplitScr $pas[$i] $MFconfig{splitFastaInput}\n";
			}
			$lowEffort = 1 if ($allowLinks &&  $lowEffort != 0);
		}
		if ($porechopFlag && $is3rdGen){
			$unzipcmd .= "touch $porechStone\n" if ($porechopFlag);
			$pas[$i] = ($porechopped);
		}
	}
	$lowEffort = 0 if ($lowEffort == -1); #FALLBACK option
	#die "$lowEffort\n";

	$unzipcmd .= "touch $finishStone\n";
	#if ($useTrimomatic && !$porechopFlag){$unzipcmd .= "touch $trimoStone\n" ;}
	my $jobN = "";
	my $jobNUZ = $jobN;
	#die;
	#die "unipss $unzipcmd\n";
	#print "  HH ".-s $testf2 < -s $testf1." FF \n";
	my $tmpCmd;
	#die "$unzipcmd\n$calcUnzp\n";
	if ($calcUnzp && !-e $finishStone && !$MFconfig{filterFromSource}){ #submit & check for files
		my $presence=1;
		for (my $i=0;$i<@pa1;$i++){	if (!-e $pa1[0] || -z $pa1[0]){$presence=0; }	}#die "$pa1[0]\n";
		for (my $i=0;$i<@pa2;$i++){	if ( !-e $pa2[0] || -z $pa2[0]){$presence=0;}	}
		#die "$presence presi\n";
		if (!$presence || !-e $finishStone  ){#|| ($useTrimomatic && !-e $trimoStone) ){
			if ($lowEffort==1){
				#die "$unzipcmd\n";
				systemW $unzipcmd;
			} elsif ($lowEffort==0) {
				$jobN = "_UZ$JNUM"; 
				$unzipcmd = $unzipcmdTMP . $unzipcmd ;

				$unzipcmd = "" if ($presence && -e $finishStone && -e $trimoStone);
				my $tmpSHDD = $QSBoptHR->{tmpSpace};
				$QSBoptHR->{tmpSpace} = $HDDspace{kraken}; #set option how much tmp space is required, and reset afterwards
				($jobN, $tmpCmd) = qsubSystem($logDir."UNZP.sh",$unzipcmd,$numCore,"20G",$jobN,$jDepe,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
				$QSBoptHR->{tmpSpace} = $tmpSHDD;
			} else {
				die "MG-TK.pl:: unzipcmd not excecutable due to \$lowEffort not correctly set ($lowEffort)!\n";
			}
			$WT += 30;
			#print " FDFS ";
		}
	#	die($jobN);
		#### 2 : remove human contamination
		#DB just needs to be loaded way too often.. do after sdm to have single files
		#$jobN = removeHostSeqs(\@pa1,\@pa2, \@pas,$tmpPath,$jobN);
	}
	#die "$unzipcmd\n";
	
		#check already here for mate pair support reads, deactivate fastp2
	my ($mateCmd,$mateSto) = ("","/sh");
#	print "@libInfo\n";
	my $ii=0; my @matePrps;
	while($ii<@libInfo){
		#print $ii." \n";
		#next;
		unless ($libInfo[$ii] =~ m/.*mate.*/i){$ii++;next;} #remove from process
		my $href;
		($href,$mateCmd,$mateSto) = check_matesL($pa1[$ii],$pa2[$ii],$finDest."mateCln/",$doMateCln);
		#remove this ori file from raw reads
		splice(@pa1,$ii,1);splice(@pa2,$ii,1);splice(@libInfo,$ii,1);
		push(@matePrps,$href);
	}
	if ( ($mateSto ne "/sh" && !-e $mateSto) && $calcUnzp){
		$mateCmd = "" if ( -e $mateSto);
		($jobN, $tmpCmd) = qsubSystem($logDir."MATE.sh",$mateCmd,1,"20G","_MT$JNUM",$jobN,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		#### 3 : if mate pairs, process these now (and remove corresponding raw fiiles
	}
	foreach my $hrr (@matePrps){
		my %nateFiles = %{$hrr};
		push(@pa1,$nateFiles{pe1});push(@pa2,$nateFiles{pe2});push(@pas,$nateFiles{se}); push(@libInfo,"pe4mt$ii");
		push(@pa1,$nateFiles{mp1});push(@pa2,$nateFiles{mp2});push(@libInfo,"mate$ii");
		push(@pa1,$nateFiles{un1});push(@pa2,$nateFiles{un2});push(@libInfo,"mate_unkn$ii");
	}
#set global var
	$inputFileSizeMB{$curSmpl} = $totalInputSizeMB;
	#die "HJASD:@pa1\n@pas\n";
	%seqSet = (pa1 => \@pa1, pa2 => \@pa2, pas => \@pas, seqTech => $seqTech, is3rdGen =>$is3rdGen,
			paX1 => \@paX1, paX2 => \@paX2, paXs => \@paXs, seqTechX => $seqTechX, is3rdGenX => $is3rdGenX,
			libInfo => \@libInfo, libInfoX => \@libInfoX,
			totalInputSizeMB => $totalInputSizeMB,
			rawReads => $rawReads,
			mmpu => $mmpu, WT => $WT,
			);
			
	if (!-e $inputRawFile || -s $inputRawFile == 0){
		open O,">$inputRawFile"; print O $seqSet{"rawReads"}; close O;
	}

	#return ($jobN,\@pa1,\@pa2, \@pas, $WT, $rawReads, $mmpu,\@libInfo, $totalInputSizeMB);
	return ($jobN, \%seqSet);
}


sub mOTU2Mapping{
	my ($cleanSeqSetHR,$tmpD,$finOutD,$smp,$Ncore,$deps) = @_;
	my $inF1a = ${$cleanSeqSetHR}{arp1}; my $inF2a = ${$cleanSeqSetHR}{arp2}; my $inFSa = ${$cleanSeqSetHR}{singAr}; #my $mergRdsHr = ${$cleanSeqSetHR}{mrgHshHR};

	
	
	#die "yes\n\n";
	my $DBdir = getProgPaths("motus2_DB",0);
	my $m2Bin = getProgPaths("motus2");#"python "."$m2Glb/$m2sub/motus";

	my $stone = $finOutD."$smp.Motu2.sto";
	my $DBstr = ""; $DBstr = "-db $DBdir " unless ($DBdir eq "");
	if (!$MFopt{DoMOTU2} ||  -e $stone){return;}
	my @car1 = @{$inF1a}; my @car2 = @{$inF2a}; my @sar = @{$inFSa};
	my $inF1 = join(",",@car1); my $inF2 = join(",",@car2); my $inFS = join(",",@sar); 
	system "mkdir -p $finOutD\n" unless (-d $finOutD);
	my $cmd = "mkdir -p $tmpD\n";
	#$cmd .= "which bwa\n";
	$cmd .= "$m2Bin profile ";
	$cmd .= "-f $inF1 -r $inF2 " if (@car1 > 0 );
	$cmd .= "-s $inFS " if (@sar > 0);
	$cmd .= " -n $smp $DBstr -q -c -t $Ncore -q | gzip -c > $finOutD/$smp.motu2.tab.gz\n";
	$cmd .= "if [ -s $finOutD/$smp.motu2.tab.gz ] ; then touch  $stone ; fi\n";
	#$cmd .= "touch  $stone\n";
	my $jobN = "mOT$JNUM";
	#die $cmd."\n";
	my ($jobN2,$tmpCmd) = qsubSystem($logDir."mOTU2_prof.sh",$cmd,$Ncore,"3G",$jobN,$deps,"",1,[],$QSBoptHR);
	$jobN  = $jobN2;
	return $jobN;
}


sub prepMOTU2(){
	my $orimotu2 = getProgPaths("motus2_DB",0);
	return "" if (!$MFopt{DoMOTU2} || $orimotu2 eq "");
	#my $m2Bin = getProgPaths("motus2");#"python "."$m2Glb/$m2sub/motus";
	my $m2sub="";
	if ($orimotu2 =~ m/\/([^\/]+)\/?$/){
		$m2sub = $1;
	} elsif (! -d $orimotu2) {
		die "$orimotu2 seems not to be a dir\n";
	}
	my $m2Glb = $runTmpDirGlobal."DB/";
	my $cmd = "";
	#print "$m2Glb/$m2sub && !-e $m2Glb/$m2sub/motu2";
	if (!-d "$m2Glb/$m2sub"){#&& !-e "$m2Glb/$m2sub/motu2"){
		$cmd .= "mkdir -p $m2Glb\n";
		$cmd .= "cp -r $orimotu2 $m2Glb\n";
	}
	#$cmd .= "exit(0);\n";
	#die $cmd;
	my $tmpCmd=""; my $jobN="";
	if ($cmd ne ""){
		$jobN= "_m2DB";
		($jobN, $tmpCmd) = qsubSystem("$baseOut/LOGandSUB/motu2DB.sh",$cmd,1,"1G",$jobN,"","",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
	}
	
	return ($jobN);
}
sub prepKraken(){
	#my ($DBdir) = @_;
	my %DBname;
	if ($MFopt{humanFilter} && $MFopt{filterHostDB1} eq ""){
		$DBname{"hum1stTry"} = 1;
	}
	#die "$MFopt{humanFilter} && $MFopt{filterHostDB1}\n";
	$DBname{"minikraken_2015"} = 1 if ($MFopt{DoEukGenePred});
	$DBname{$MFopt{globalKraTaxkDB}} = 1 if ($MFopt{globalKraTaxkDB} ne "");
	my $cmd = "";my $jobN= "";
	return $jobN if (keys %DBname == 0);
	my $oriKrakDir = getProgPaths("Kraken_path_DB");
	#die "krakper\n";
	foreach my $kk (keys %DBname ){
		if (!-d "$oriKrakDir$kk"){die "can't find kraken db $oriKrakDir$kk\n";}
		if (!-d "$krakenDBDirGlobal/$kk" && !-e "$krakenDBDirGlobal/$kk/cpFin.stone" ){
			$cmd =  "cp -r $oriKrakDir$kk $krakenDBDirGlobal/\n";
			$cmd .= "touch $krakenDBDirGlobal/$kk/cpFin.stone\n";
		}
	}
	my $tmpCmd="";
	if ($cmd ne ""){
		$jobN= "_KRDB";
		($jobN, $tmpCmd) = qsubSystem("$baseOut/LOGandSUB/krkDB.sh",$cmd,1,"1G",$jobN,"","",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		#print "Fix KRDB 2416\n";
	}
	#die $cmd;
	return $jobN;
}

sub removeHostSeqs($ $ $ $){
	my ($cleanSeqSetHR,$tmpD,$jDep,$checkIfExists) = @_;
	return $jDep unless ($MFopt{humanFilter});
	my @pa1 = @{${$cleanSeqSetHR}{arp1}}; my @pa2 = @{${$cleanSeqSetHR}{arp2}}; my @pas = @{${$cleanSeqSetHR}{singAr}};
	my $outputExists=1;
	my $fileDir = "";

	if ($fileDir eq ""){
		if (@pa1 > 0){
			$pa1[0] =~ m/(.*\/)[^\/]+$/; $fileDir= $1 ;
		} else {#has to be in pas then
			$pas[0] =~ m/(.*\/)[^\/]+$/; $fileDir= $1 ;
		}
	}
	my $hostRMVer=$MFopt{humanFilter}; #0: no, 1:kraken2, 2: kraken1, 3:hostile
	#files don't need to be checked.. it's in any case just sdm files..
	$outputExists = 0 if (!-e "$fileDir/krak.stone");
	#die "$outputExists && $checkIfExists\n";
	return $jDep if ($outputExists && $checkIfExists);
	
	#flag whether to use kraken v1 or v2
	my $krk2Bin = getProgPaths("kraken2");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
	my $krkBin = "";
	my $hostileBin = getProgPaths("hostile");
	my $hostileDB = getProgPaths("hostileDB");
	#my $krak1= 1;	if ($krk2Bin ne ""){ $krak1=0; }
	if ($hostRMVer==2){
		$krkBin = getProgPaths("kraken");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
	}

	#my ($DBdir,$DBname) = @_;
	
	my $DBdir = $krakenDBDirGlobal."/";
	my $unsplBin = getProgPaths("unsplitKrak_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/unsplit_krak.pl ";
	my @DBname = ("hum1stTry"); my $numThr = 4;
	if (@filterHostDB > 0){
		@DBname = (); my $cnt=0;
		foreach my $fhdb (@filterHostDB){
			if (!-d $fhdb){die "-filterHostKrakDB is not a dir: \"$fhdb\"\n";}
			$DBname[$cnt] = $fhdb; $DBdir=""; $cnt++;
		}
	}
	my $cmd = "\n\nmkdir -p $tmpD\n\n"; my $tmpF ="$tmpD/krak.tmp.fq";
	
	#loop around different DBs..
	
	for (my $j=0;$j<@DBname ; $j++){
		for (my $i=0;$i<@pa1;$i++){
			my $r1 = $pa1[$i]; my $r2 = $pa2[$i];
			#if ($fileDir eq ""){
			#	$pa1[$i] =~ m/(.*\/)[^\/]+$/; $fileDir= $1 ;
			#	$cmd .= "\n\nrm -f $fileDir/krak.stone\n\n";
			#}
			my $gzFlag = "";$gzFlag =  "--gzip-compressed" if ($pa1[$i] =~ m/\.gz$/);
			my $gzEnd = ""; $gzEnd = ".gz" if ($gzFlag ne "");

			if ($hostRMVer==2){
				$cmd .= "$krkBin --preload --threads $numThr --fastq-input $gzFlag --unclassified-out $tmpF --db $DBdir$DBname[$j]  $r1 $r2 > /dev/null\n";
				#overwrites input files
				$r1 =~ s/\.gz$//; $r2 =~ s/\.gz$//;
				$cmd .= "$unsplBin $tmpF $r1 $r2\nrm -f $tmpF*\n";
				if ($gzFlag ne ""){
					$cmd .= "$pigzBin -f -p $numThr $r1 $r2\n";
				}
			} elsif($hostRMVer==1) { #kraken2
				$tmpF ="$tmpD/krak.tmp#.fq"; my $tmpF1 = "$tmpD/krak.tmp_1.fq";my $tmpF2 = "$tmpD/krak.tmp_2.fq";
				$cmd .= "$krk2Bin --threads $numThr $gzFlag --paired --unclassified-out $tmpF --db $DBdir$DBname[$j] --output - $MFopt{filterHostKr2QuickMode} --confidence $MFopt{krakHostConf} $r1 $r2 \n";
				$cmd .= "$pigzBin -f -p $numThr $tmpF1 $tmpF2\n" unless ($gzFlag eq "");
				$cmd .= "rm -f $r1 $r2; \n";
				$cmd .= "mv $tmpF1$gzEnd $r1; mv $tmpF2$gzEnd $r2\n";
			} elsif($hostRMVer==3) { #hostile
				system "rm $hostileDB/human-t2t-hla.mmi" if (-z "$hostileDB/human-t2t-hla.mmi");
			
				$cmd .= "export HOSTILE_CACHE_DIR=$hostileDB\n";
				$cmd .= "$hostileBin clean --fastq1 $r1 --fastq2 $r2 --index $hostileDB/$MFopt{hostileIndex} --aligner auto --output $tmpD/ --threads $numThr --airplane --force\n";
				my $newR1 = $r1; my $newR2 = $r2; 
				$newR1 =~ s/.*\///; $newR2 =~ s/.*\///; #remove path
				$newR1 =~ s/\.fq([\.gz])?/\.clean_1\.fastq$1/;$newR2 =~ s/\.fq([\.gz])?/\.clean_2\.fastq$1/;
				$cmd .= "rm $r1 $r2\nmv $tmpD/$newR1 $r1;\nmv $tmpD/$newR2 $r2\n";
				
				#die $cmd."\n";
			}
		}
		for (my $i=0;$i<@pas;$i++){
			my $rs = $pas[$i]; 
			my $gzFlag = "";$gzFlag =  "--gzip-compressed" if ($pas[$i] =~ m/\.gz$/);
			if ($hostRMVer==2){
				$cmd .= "$krkBin --preload --threads $numThr --fastq-input $gzFlag --unclassified-out $tmpF --db $DBdir$DBname[$j]  $rs > /dev/null\n";
			} elsif ($hostRMVer==1) {
				$cmd .= "$krk2Bin --threads $numThr $gzFlag --unclassified-out $tmpF --db $DBdir$DBname[$j] --output - $MFopt{filterHostKr2QuickMode} --confidence $MFopt{krakHostConf} $rs \n";
				if ($gzFlag eq ""){$cmd .= "rm -f $rs; mv $tmpF $rs\n";
				} else{$cmd .= "rm -f $rs; $pigzBin -f -p $numThr -c $tmpF > $rs\nrm -f $tmpF\n";}
			} elsif($hostRMVer==3) {#hostile
				system "rm $hostileDB/human-t2t-hla.mmi" if (-z "$hostileDB/human-t2t-hla.mmi");
				$cmd .= "export HOSTILE_CACHE_DIR=$hostileDB\n";
				$cmd .= "$hostileBin clean --fastq1 $rs --index $MFopt{hostileIndex} --aligner auto --output $tmpD/ --threads $numThr --airplane --force\n";
				my $newR1 = $rs;
				$newR1 =~ s/.*\///;$newR1 =~ s/\.fq([\.gz])?/\.clean\.fastq$1/;
				$cmd .= "rm $rs \nmv $tmpD/$newR1 $rs;\n";
			}
			#overwrites input files
		}
	}
	$cmd .= "\n\n";
	$cmd .= "touch $fileDir/krak.stone\n";
	my $jobN = ""; my $tmpCmd="";
	unless (-e "$fileDir/krak.stone"){
		#die "$fileDir\n";
		my $tmpSHDD = $QSBoptHR->{tmpSpace};
		my $reqSpace = $HDDspace{kraken};
		if ($inputFileSizeMB{$curSmpl} > 10000){
			#convert input to a) Gb and b) *6 to account for extra space needed by kraken
			$reqSpace = int($inputFileSizeMB{$curSmpl}*6/1024)+30  ."G";
		} 
		#die "krakHS : $reqSpace $inputFileSizeMB{$curSmpl} \n\n";
		$QSBoptHR->{tmpSpace} = $reqSpace; #set option how much tmp space is required, and reset afterwards
		#### 1 : UNZIP
		$jobN = "_KR$JNUM";
		($jobN, $tmpCmd) = qsubSystem($logDir."KrakHS.sh",$cmd,$numThr,"16G",$jobN,$jDep.";$krakDeps","",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
			$QSBoptHR->{tmpSpace} = $tmpSHDD;
	}
	return $jobN;
}

sub krakHSapSingl($ $){ #just on single file..
	my ($inF, $outF,$numThr) = @_;
	return "" if ( -e $outF);
	my $krk2Bin = getProgPaths("kraken2");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
	#my ($DBdir,$DBname) = @_;
	my $DBdir = $krakenDBDirGlobal."/";
	#my $unsplBin = getProgPaths("unsplitKrak_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/unsplit_krak.pl ";
	my $DBname = "hum1stTry";# my $numThr = 6;
	my $cmd  = ""; my $doGZIP=0;
	if ($outF =~ m/\.gz$/){ $doGZIP = 1;
		$outF =~ s/\.gz//;
	}
	my $gzFlag = "";$gzFlag =  "--gzip-compressed" if ($inF =~ m/\.gz$/);
	#flag whether to use kraken v1 or v2
	my $krak1= 1;	if ($krk2Bin ne ""){ $krak1=0; }
	my $krkBin = "";
	$krkBin = getProgPaths("kraken") if ($krak1);
	if ($krak1){
		$cmd .= "$krkBin --preload --threads $numThr --fastq-input $gzFlag --unclassified-out $outF --db $DBdir$DBname  $inF > /dev/null\n";
	} else {
		$cmd .= "$krk2Bin --threads $numThr $gzFlag --unclassified-out $outF --db $DBdir$DBname --output - $MFopt{filterHostKr2QuickMode} --confidence $MFopt{krakHostConf} $inF \n";
	}
	if ($doGZIP){
		$cmd .= "$pigzBin -f -p $numThr $outF \n";
	}
	return $cmd;
}

sub genoSize(){
	my ($cleanSeqSetHR,$oD,$jdep) = @_;
	my $arp1 = ${$cleanSeqSetHR}{arp1}; my $arp2 = ${$cleanSeqSetHR}{arp2}; my $ars = ${$cleanSeqSetHR}{singAr}; my $mergRdsHr = ${$cleanSeqSetHR}{mrgHshHR};

	my @pa1 = @{$arp1}; my @pa2 = @{$arp2}; my @pas = @{$ars};
	my $cmd = "mkdir -p $oD\n";
	my $microCensBin = getProgPaths("microCens");
	for (my $i=0;$i<@pa1;$i++){
		$cmd .= "$microCensBin -t 2 $pa1[$i],$pa2[$i] $oD/MC.$i.result\n";
	}
#	die $cmd;
	my $jobName = "_GS$JNUM"; my $tmpCmd;
	
	($jobName,$tmpCmd) = qsubSystem($logDir."MicroCens.sh",$cmd,2,"40G",$jobName,$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
	return $jobName;
}
sub krakenTaxEst(){
	my ($cleanSeqSetHR, $outD, $tmpD,$name,$jobd) = @_;
		my $arp1 = ${$cleanSeqSetHR}{arp1}; my $arp2 = ${$cleanSeqSetHR}{arp2}; my $ars = ${$cleanSeqSetHR}{singAr}; 

	my @pa1 = @{$arp1}; my @pa2 = @{$arp2}; my @pas = @{$ars};
	#$outD.= "$MFopt{globalKraTaxkDB}/";
	my $krakStone = "$outD/krakDone.sto";
	return $jobd if (-d $outD && -e $krakStone);
	my $numCore = $MFopt{krakenCores};
	my @thrs = (0.01,0.02,0.04,0.06,0.1,0.2,0.3);
	my $cmd = "mkdir -p $outD\nrm -rf $tmpD\nmkdir -p $tmpD\n";
	my $it =0;
	my $curDB = "$krakenDBDirGlobal/$MFopt{globalKraTaxkDB}";

	my $krkBin = "";my $krak1 = 1;
	$krkBin = getProgPaths("kraken") if ($krak1);

	#die $curDB."\n";
	#paired read tax assign
	for (my $i=0;$i<@pa1;$i++){
		my $r1 = $pa1[$i]; my $r2 = $pa2[$i];
		$pa1[$i] =~ m/(.*\/)[^\/]+$/; 
		$cmd .= "$krkBin --paired --preload --threads $numCore --fastq-input  --db $curDB  $r1 $r2 >$tmpD/rawKrak.$it.out\n";
		for (my $j=0;$j< @thrs;$j++){
			$cmd .= "$krkBin-filter --db $curDB  --threshold $thrs[$j] $tmpD/rawKrak.$it.out | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak_$thrs[$j]"."_$it.out\n";
		}
		$it++;
	}
	$cmd .= "\n\n";
	#single read tax assign
	for (my $i=0;$i<@pas;$i++){
		my $rs = $pas[$i]; 
		$cmd .= "$krkBin --preload --threads $numCore --fastq-input  --db $curDB  $rs >$tmpD/rawKrak.$it.out\n";
		#$cmd .= " | tee ";#$krkBin-filter --db $curDB  --threshold 0.01 | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak$it.out\n";
		for (my $j=0;$j< @thrs;$j++){
			$cmd .= "$krkBin-filter --db $curDB  --threshold $thrs[$j] $tmpD/rawKrak.$it.out | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak_$thrs[$j]"."_$it.out\n";
		}
		#overwrites input files
		$it++;
	}
	
	#TODO: 1: make table; 2: copy to outD
	$cmd .= "\n\n";
	my $krakCnts1 = getProgPaths("krakCnts_scr");
	for (my $j=0;$j< @thrs;$j++){
		$cmd .= "cat $tmpD/krak_$thrs[$j]"."_*.out > $tmpD/allkrak$thrs[$j].out\n";
		$cmd .= "$krakCnts1 $tmpD/allkrak$thrs[$j].out $outD/krak.$thrs[$j].cnt.tax\n";
	}

	$cmd .= "touch $krakStone\n";
	$cmd .= "rm -fr $tmpD";

	my $jobName = "_KT$JNUM"; my $tmpCmd;
	if (!-d $outD || !-e $krakStone){
		($jobName,$tmpCmd) = qsubSystem($logDir."KrkTax.sh",$cmd,$numCore,"20G",$jobName,$jobd.";$krakDeps","",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
	}
	return $jobName;
}

sub check_map_done{
	my ($doCram, $finalD, $baseN, $mappDir) = @_;
	
	#my $aa = (stat "$finalD/$baseN-smd.bam")[7];
	#die "$mappDir/$baseN-smd.bam "."\n";

	my $retVal = 1;
	if ($doCram && (!-e "$finalD/$baseN-smd.cram.sto" )){# && !-e "$mappDir/$baseN-smd.cram.sto"  ) ){#&& !$params{bamIsNew} ) ){
		$retVal = 0;
	} elsif (!$doCram && !-e "$finalD/$baseN-smd.bam" ){#&& !$params{bamIsNew}){
		$retVal = 0;
	}
	
	return $retVal;
}
sub check_depth_done{
	my ($doCram, $finalD, $baseN, $mappDir) = @_;
	my $retVal = 1;
	if ( !-e "$mappDir/$baseN-smd.bam.coverage.gz" ){#&& (-e "$mappDir/$baseN-smd.bam" || -e "$mappDir/$baseN-smd.cram") ){# && !$params{bamIsNew}){ #already stored in mapping dir, still needs to be copied
		#die "$mappDir/$baseN-smd.bam.coverage.gz";
		$retVal = 0;
	}
	return $retVal;
}

#create string for mapper to register libraries
sub getRgStr{ 
	my ($smpl,$libsOri,$libsOriX,$usePairs,$mapper) = @_;
	my $rgStr ="noReg";
	if ($mapper > 1){ #bwa/minimap2 have same format..
		$rgStr = "'\@RG\\tID:$1\\tSM:$smpl\\tPL:ILLUMINA";$rgStr .= "\\tLB:lib1'";
	}
	if ($mapper==1){ #bowtie2
		$rgStr = "--rg SM:$smpl --rg PL:ILLUMINA "; #PU:lib1
		if ($usePairs){
			$rgStr .= "--rg LB:$libsOri ";
			$rgStr .= " -X $MFconfig{mateInsertLength} " if ($libsOri =~ m/mate/);
		} else {
			$rgStr .= "--rg LB:$libsOriX ";
			$rgStr .= " -X $MFconfig{mateInsertLength} " if ($libsOriX =~ m/mate/);
		}
	}
	return $rgStr;
}



sub getAlgnCmdBase{
	my ($MapperProg,$NcoreL, $readTec,$mapModeTogether,$unaligned) = @_;
	my $algCmdBase = "";
	if ($MapperProg==1){ #bowtie2
		my $bwt2Bin = getProgPaths("bwt2");
		$algCmdBase = "$bwt2Bin  --end-to-end -p $NcoreL  "; #--no-unal
		if ($mapModeTogether == 2 || $mapModeTogether == -1){$algCmdBase .= "-a ";}
	} elsif ($MapperProg==2){
		my $bwaBin = getProgPaths("bwa");
		$algCmdBase = "$bwaBin mem -t $NcoreL ";
	} elsif ($MapperProg==3){#minimap2 
		my $mini2Bin = getProgPaths("minimap2");
		$algCmdBase = "$mini2Bin -2 -t $NcoreL --secondary=no ";
		if ($readTec eq "ONT"){
			$algCmdBase .= " -x map-ont" ; #use nanopore optimzed for now..
		} elsif ($readTec eq "PB"){
			$algCmdBase .= " -x map-pb" ; #use nanopore optimzed for now..
		} else {
			print"Warning: Minimap2 used for short reads; not recommended\n";
			$algCmdBase .= " -x sr" ;
		}
		#$algCmdBase .= " --sam-hit-only " unless ($unaligned); #deactivated as important for counting
	} elsif ($MapperProg==4){ #kma 
		my $kmaBin = getProgPaths("kma");
		my $consID = 0.95;my $minPhred = 15; my $minMapQ = 20; my $minQueryCov = 0.2; #-bc $minPhred -tmp $nodeTmp/${baseN}.kmatmp/
		$algCmdBase = "$kmaBin -nc -na -nf -sam 4 -apm p -mrc $minQueryCov -mq $minMapQ -bcd 1 -ID $consID -ref_fsa -t $NcoreL   ";
		if ($readTec eq "ONT"){$algCmdBase .= " -bcNano -ont " ;
		} elsif ($readTec eq "PB"){$algCmdBase .= " -mint3  " ;
		} else {$algCmdBase .= " -mint2 " ;
		}
	} elsif ($MapperProg==5){ #strobealign
		my $stroBin = getProgPaths("strobealign");
		$algCmdBase = "$stroBin -t $NcoreL --no-progress ";
		#$algCmdBase .= " -U " unless ($unaligned);

	} else {
		die "getAlgnCmdBase:: unknown mapper \"$MapperProg\"\nAborting..\n";
	}
	return $algCmdBase;
}



sub announce_MGTK{
	
	
#	print "888b     d888        d8888 88888888888     d8888 8888888888 8888888 888      8888888888 8888888b.  \n";#
#	print "8888b   d8888       d88888     888        d88888 888          888   888      888        888   Y88b \n";
#	print "88888b.d88888      d88P888     888       d88P888 888          888   888      888        888    888 \n";
#	print "888Y88888P888     d88P 888     888      d88P 888 8888888      888   888      8888888    888   d88P \n";
#	print "888 Y888P 888    d88P  888     888     d88P  888 888          888   888      888        8888888P\"  \n";
#	print "888  Y8P  888   d88P   888     888    d88P   888 888          888   888      888        888 T88b   \n";
#	print "888   \"   888  d8888888888     888   d8888888888 888          888   888      888        888  T88b  \n";
#	print "888       888 d88P     888     888  d88P     888 888        8888888 88888888 8888888888 888   T88b \n";

#	print "/------------------------------------------------------------------------\\\n";
#	print "|  _______ _______ _______ _______ _______ _____        _______  ______  |\n";
#	print "|  |  |  | |_____|    |    |_____| |______   |   |      |______ |_____/  |\n";
#	print "|  |  |  | |     |    |    |     | |       __|__ |_____ |______ |    \\_  |\n";
#	print "|                                                                        |\n";
	print "/--------------------------------------------\\\n";
	print "|  ███╗   ███╗ ██████╗    ████████╗██╗  ██╗  |\n";
	print "|  ████╗ ████║██╔════╝    ╚══██╔══╝██║ ██╔╝  |\n";
	print "|  ██╔████╔██║██║  ███╗█████╗██║   █████╔╝   |\n";
	print "|  ██║╚██╔╝██║██║   ██║╚════╝██║   ██╔═██╗   |\n";
	print "|  ██║ ╚═╝ ██║╚██████╔╝      ██║   ██║  ██╗  |\n";
	print "|  ╚═╝     ╚═╝ ╚═════╝       ╚═╝   ╚═╝  ╚═╝  |\n";
	print "\\--------------------------------------------/\n";

	print "This is MG-TK v$MATFILER_ver\n";
}




sub getMapProgNm{
	my ($MapperProg) = @_;
	my $mapProgNm = "undefined mapper";
	if ($MapperProg==1){$mapProgNm = "bowtie2";
	}elsif ($MapperProg==2){$mapProgNm = "bwa";
	}elsif ($MapperProg==3){$mapProgNm = "minimap2";
	}elsif ($MapperProg==4){$mapProgNm = "kma";
	}elsif ($MapperProg==5){$mapProgNm = "strobealign";
	} else {print "coult not recognize mapper $MapperProg!!\n";
	}
	return $mapProgNm
}

#			my %postTreat = (MapperProg=>$MapperProg, readTec=>$readTec, NcoreL => $NcoreL,nodeTmp =>$nodeTmp,tmpOut21=>$tmpOut21, 
#										subBamsAR => \@subBams,unaligned =>$unaligned,baseN => $baseN,xtraSamSteps1 => $xtraSamSteps1,
#										decoyModeActive => $decoyModeActive, map2ndTogether => $map2ndTogether,regsAR => \@regs, 
#										doCram => $doCram, finalDSar => \@finalDS, outNmsAR => \@outNms, mappDirAR => \@mappDir, 
#										reg_lcsAR => \@reg_lcs)
sub alignPostTreat{
	my ($postTreatHR, $i, $kk) = @_;
	my %postTreat = %{$postTreatHR};
	my $bamHdFilt_scr = getProgPaths("bamHdFilt_scr");
	
	my $algCmd = "";
	my ($MapperProg, $readTec, $NcoreL,$nodeTmp,$tmpOut21, $subBamsAR,$unaligned,$baseN) = ($postTreat{MapperProg}, $postTreat{readTec}, $postTreat{NcoreL},$postTreat{nodeTmp}, $postTreat{tmpOut21}, $postTreat{subBamsAR},$postTreat{unaligned},$postTreat{baseN});
	my $xtraSamSteps1 = $postTreat{xtraSamSteps1}; 
	
	my $bamfilter = getProgPaths("bamFilter_scr");
		#my $filterStep="";
	if ( ($MapperProg==3 || $MapperProg==5 ) && ($readTec eq "ONT" || $readTec eq "PB")){ #low id long reads...
		$algCmd .= " | $bamfilter $MFopt{bamfilterPB}  ";
	} else { #illumina parameters
		$algCmd .= " | $bamfilter $MFopt{bamfilterIll} ";
	}


	
	my $iTO = "$tmpOut21.$i.$kk";
	if ($unaligned ne "" ){ #basic sort of unmapped reads via samtools (general purpose step)
		$algCmd .= " | $smtBin view -b1 -@ $NcoreL > $iTO.t\n";
		$algCmd .= "$smtBin view -u -h $iTO.t | $xtraSamSteps1 $smtBin view -b1 -@ $NcoreL -F 4 -  > $iTO\n";
		#sort out unaligned reads
		$algCmd .= "$smtBin view -u -h -@ $NcoreL -f 4 $iTO.t | $smtBin fastq -1 $nodeTmp/$baseN.$i.1.fq.gz -2 $nodeTmp/$baseN.$i.2.fq.gz -s $nodeTmp/$baseN.$i.s.fq.gz - \n";
		#and copy them already to final destination.. no reason to keep them around..
		$algCmd .= "cat $nodeTmp/$baseN.$i.1.fq.gz >> $unaligned/unal.1.fq.gz;\ncat $nodeTmp/$baseN.$i.2.fq.gz >> $unaligned/unal.2.fq.gz;\n cat $nodeTmp/$baseN.$i.s.fq.gz >> $unaligned/unal.fq.gz;\n";
		#and remove all the temp files..
		$algCmd .= "rm -f $nodeTmp/$baseN.$i.*fq.gz $iTO.t\n";
	} else {
		$algCmd .= " | $xtraSamSteps1  $smtBin view -b1 -@ $NcoreL -F 4 - > $iTO\n";
	}
	
	if ($postTreat{decoyModeActive} || $postTreat{map2ndTogether}>0){#remove unnecessary reads & filter specific reads for each ref genome into separate bam
		$algCmd .= "$smtBin index $iTO\n";
		for (my $k=0;$k<@{$postTreat{regsAR}};$k++){
			#		print "$k\t$finalDS[$k]\n";
			if(check_map_done(${$postTreat{doCram}}, ${$postTreat{finalDSar}}[$k], ${$postTreat{outNmsAR}}[$k], ${$postTreat{mappDirAR}}[$k])){$$subBamsAR[$k]="";next;}
			$algCmd .= "\n\nset +e \n" if ($k==0);
			$algCmd .= "#  %%%%%%%%%%%%%%%% $k %%%%%%%%%%%%%%%% \n";
			$algCmd .= "$bamHdFilt_scr $iTO ${$postTreat{reg_lcsAR}}[$k] 0 > $iTO.decoy.sam.$k\n";
			$algCmd .= "sleep 0.05\n";
			$algCmd.= "  $smtBin view -@ $NcoreL $iTO ${$postTreat{regsAR}}[$k] >> $iTO.decoy.sam.$k\n";
			$algCmd.= "  $smtBin view -b -h -@ $NcoreL $iTO.decoy.sam.$k > $iTO.decoy.$k\nrm -f $iTO.decoy.sam.$k \n";
			$$subBamsAR[$k] .= " $iTO.decoy.$k";  
		}
		
		$algCmd.= "rm -f $iTO\n";# mv $tmpOut21.decoy.$i $tmpOut21.$i\n";
		$algCmd .= "\n\nset -e \n";

	} else {
		$$subBamsAR[0] .= " $iTO"; 
	}
	
	#die "$algCmd\n";
	return  ($algCmd,$subBamsAR);
}


sub fastCovCalc{
	my ($AsgHR,$ASG) = @_;
	my $cmd = "";
	#my ($par1,$par2,$parS,$liar,$rear) = getRawSeqsAssmGrp($AsgHR,$ASG,0);

	#$cmd .= getAlgnCmdBase(1,$NcoreL, $readTec,1);
	#$algCmd .= "$algCmdBase -x $bwtIdxs[$kk] -1 ".join(",",@{$par1}) ." -2 ".join(",",@{$par2});

}




sub mapReadsToRef{
	my ($dirsHr, $AsgHR, $jDepe) = @_;
	#my ($par1,$par2,$parS,$liar,$rear) = getRawSeqsAssmGrp($AsgHR,$ASG,$supportRds);
	#wrong: needs to be reads of the sample, not assembly group!
	my $outName = $dirsHr->{smplName};my $ASG = $dirsHr->{assGrp};
	my $is2ndMap = $dirsHr->{is2ndMap};
	my $doCram = $dirsHr->{cramAlig};
	my $immediateSubm = $dirsHr->{submNow};
	my $unaligned = $dirsHr->{unalDir}; #store unaligned reads somewhere ? leave "" to deactivate
	my $Ncore = $dirsHr->{mapCores};
	my $supportRds = $dirsHr->{mapSupport};
	my $supTag = ""; if ($supportRds){$supTag = ".sup";}
	#die "$supTag $supportRds $dirsHr->{mapSupport} \n";

	my $REF = $dirsHr->{sbj}; #target to map onto, can by ","-spearated list
	my ($par1,$par2,$parS,$liar,$rear) = getRawSeqsAssmGrp($AsgHR,$ASG,$supportRds,$outName);
	my @libsOri = @{$liar};
	#simple rule for mapper program: for now set to bowtie2

	#die "$REF\n";
			# ($dirsHr,$outName, $is2ndMap, $par1,$par2,$Ncore,#hm,m,x,x,x,x,x
	#	$REF,$jDepe, $unaligned,$immediateSubm,#m,x,x,x,m
	#	$doCram,$smpl,$libAR) = @_;#x,x,x
	#get mapper progs..
	#my $bwaBin = "";
	# if ($MapperProg==2);

		#get essential dirs...
	my $mappDirPre = ${$dirsHr}{glbMapDir};	my $nodeTmp = ${$dirsHr}{nodeTmp}; #m,s
	$nodeTmp.="_map${supTag}/";
	my $qdir = $logDir; $qdir = ${$dirsHr}{qsubDir} if (exists( ${$dirsHr}{qsubDir} ));
	my $readTec = ${$dirsHr}{readTec};
	if ($readTec eq "" && @libsOri){$readTec = $libsOri[0];}
	my $tmpOut = ${$dirsHr}{glbTmp};	my $finalD = ${$dirsHr}{outDir}; #s,m
	my $mapperProgLoc = decideMapper($MFopt{MapperProg},$readTec);
	#die "$finalD\n";
	my @finalDS = split /,/,$finalD;
	my @mappDir = split /,/,$mappDirPre;

	#my $bwtIdx = $REF."$MFcontstants{bwt2IdxFileSuffix}";
	my @bwtIdxs;
	#already exists
	#print "IS${unaligned}SI\n";
	my $baseN = "$outName$supTag";#$RNAME."_".$QNAME;
	
	my @outNms = split /,/,$outName;
	if ($outName ne "" ){
		$baseN = $outNms[0].$supTag;
	} else {
		@outNms = ($baseN);
	}
	if ($supportRds){for (my $ii=0;$ii<@outNms; $ii++){$outNms[$ii] = $outNms[$ii].$supTag;}}
	#die "$baseN @outNms\n";
	my @tmpOut22 = ($tmpOut."/$baseN.iniAlignment.bam");
	my @tmpOutxtra = ($nodeTmp."/$baseN.iniAlignment.xtra");
	#global value overwrites local value
	if ($doCram){$doCram = $MFopt{doBam2Cram};}
	my @pa1 = @{$par1}; my @pa2 = @{$par2}; my @paS = @{$parS};

	#calculate total input size (to get handle on req disk space
	

	
	my %params;
	my $bamFresh = 0; #is the bam newly being created?
	my $decoyModeActive=0; #decoy mapping
	my $map2ndTogether = $MFopt{mapModeTogether}; #map competetively among all reference genomes provided might change mapping result, if other genome set is used)
	$decoyModeActive=1 if ( $MFopt{DoMapModeDecoy} && exists($make2ndMapDecoy{Lib}) && -e $make2ndMapDecoy{Lib});
	my $isSorted = 0;		$isSorted=1 if (@pa1==1 && $MFopt{DoMapModeDecoy} && $decoyModeActive);
#	die $REF;
	if (!$decoyModeActive){
		@bwtIdxs = split /,/,$REF;
		for (my $kk=0;$kk<@bwtIdxs;$kk++){
			$bwtIdxs[$kk] .= $MFcontstants{bwt2IdxFileSuffix};
		}
	}
	#die "$isSorted\n";
	my $anyUsedPairs= 1; $anyUsedPairs = 0 if (scalar @pa1 == 0);
	$params{sortedbam}=$isSorted; $params{bamIsNew} = $bamFresh; $params{is2ndMap} = $is2ndMap;
	$params{immediateSubm} =  $immediateSubm; $params{usePairs} = $anyUsedPairs;

	my $outputExistsNEx = 0;
	for (my $k=0;$k<@outNms;$k++){
		$tmpOut22[$k] = $tmpOut."/$outNms[$k].iniAlignment.bam"; 
		$tmpOutxtra[$k] = $tmpOut."/$outNms[$k].iniAlignment.xtra"; 
		#print "$tmpOut22[$k]\n";
	}
	for (my $k=0;$k<@outNms;$k++){
		if ($MFopt{MapRewrite2nd}){$outputExistsNEx++; system "rm -fr $tmpOut22[$k] $finalDS[$k]/$outNms[$k]-smd* $mappDir[$k]/$outNms[$k]-smd*"; next;}
		my $outstat = check_map_done($doCram, $finalDS[$k], $outNms[$k], $mappDir[$k]);
		next if ($outstat);#-e "$finalDS[$k]/$outNms[$k]-smd.bam.coverage.gz" );#|| -e "$mappDir[$k]/$outNms[$k]-smd.bam.coverage.gz");
		#print "-e $finalDS[$k]/$outNms[$k]-smd.bam.coverage.gz || -e $mappDir/$outNms[$k]-smd.bam.coverage.gz";
		if (-e $tmpOut22[$k] && -s $tmpOut22[$k] < 100){unlink $tmpOut22[$k];}
		if (!-e $tmpOut22[$k]){ $outputExistsNEx++; }
	}
	#die "@outNms \n$outputExistsNEx\n";
	if ($outputExistsNEx == 0){
		#die;
		return ("","",\%params);
	} 
	
	my $NcoreL = int($Ncore);#correction for threads
	#my $nxtCRAM = "$tmpOut/$baseN-smd.cram";
	my $tmpOut21 = $nodeTmp."/$baseN.iniAlignment.bam";
	#my $sortTMP = $nodeTmp."/$baseN.srt";supportRds
	#print "$nxtBAM\n";
	#die("too far\n");
	
	my $retCmds="";my $tmpCmd=""; my $xtraSamSteps1="";
	#my $nxtCRAM = "$mappDir/$baseN-smd.cram";

	system("mkdir -p ".join(" ",@mappDir)." " .join(" ",split(/,/,$tmpOut)));
	#move ref DB & unzip
	#system("rm -r $tmpOut\n mkdir -p $tmpOut\n cp $REF $tmpOut");
	#$REF = $tmpOut.basename($REF);
	my $unzipcmd = "";
	if ($REF =~ m/\.gz$/){$unzipcmd .= "gunzip $REF\n";$REF =~ s/\.gz$//;}
	my $jobN = "";
	my $tmpUna = $tmpOut."/unalTMP/";
	#Error: No EOF block on /g/scb/bork/hildebra/SNP/MeHiAss/MH0411//tmp//mapping/Alignment.bam, possibly truncated file.
	my $bashN = "";	if ($is2ndMap){$bashN = "$outName"; $bashN =~ s/,/./g; if (length($bashN)>50){$bashN=substr($bashN,0,40)."_etc";}}
	if (0&& -e $qdir.$bashN."bwtMap2.sh.etxt"){open I,"<$qdir/".$bashN."bwtMap2.sh.etxt" or die "Can't open old bowtie_2 logfile $qdir\n"; my $str = join("", <I>); close I;
		if ($str =~ /Error: No EOF block on (.*), possibly truncated file\./){system("rm -f $1 ".join (" ",@tmpOut22) );}	
		close I;	
	}
	#depending on setup, mappDir == tmpOut
	my $algCmd = "";
	$algCmd .= "\nrm -rf $tmpOut $nodeTmp\nmkdir -p $tmpOut\nmkdir -p $nodeTmp\n"; 
	#die "$algCmd\n@mappDir\n";
	#my $algCmd = "rm -rf ".join(" ",split(/,/,$mappDir))."\nmkdir -p ".join(" ",split(/,/,$mappDir))." $tmpOut $nodeTmp\n"; 
	$algCmd .= "mkdir -p $qdir/mapStats/\n";
	my $statsF = $qdir."mapStats/".$bashN."mapStats.txt";
	
	my @regs;#subset of DB seqs to filter for 
	my @reg_lcs;
	$REF =~ m/([^\/]+)$/; my $REFnm = $1; $REFnm =~ s/\.f.*a$//; 
	#die "$REF\n$REFnm\n";
	if ($decoyModeActive || $map2ndTogether){
		#first create a new ref DB, including the targets and the assemblies from this dir
		my $bwtIdx = "$nodeTmp/$REFnm.decoyDB.fna";
		if ($map2ndTogether){
			$bwtIdx = "$map2ndTogRefDB{DB}";
		} else {
			my $decoyDBscr = getProgPaths("decoyDB_scr"); 
			$algCmd.= "\n$decoyDBscr $REF $make2ndMapDecoy{Lib} $bwtIdx $NcoreL $outName $finalD\n";
		}
		#die "$algCmd\n";
		$bwtIdx .= $MFcontstants{bwt2IdxFileSuffix};
		
		$xtraSamSteps1 = "$smtBin sort -@ $NcoreL -T $nodeTmp./$baseN.srt - |";
		@regs = @{$make2ndMapDecoy{regions}};
		@reg_lcs =  @{$make2ndMapDecoy{region_lcs}};
		@bwtIdxs = ($bwtIdx);
	} 
	my $bwt2DBsuf = "bt2"; $bwt2DBsuf = "bt2l" if ($MFopt{largeMapperDB}); 
	$algCmd .= "if [ ! -e $bwtIdxs[0].1.$bwt2DBsuf ] ;then\n	echo \"Could not find assembly bowtie2 index: $bwtIdxs[0].1.$bwt2DBsuf\";\n	exit 23;\nfi\n" if ($mapperProgLoc == 1); #needs to exit on error

	#die "@bwtIdxs\n";
	
	my $cntAli=0; my $totlRefs = scalar(@bwtIdxs);
	my $numLib = scalar @pa1 + scalar @paS; 
	#if there's too many refgenomes, copy reads onto tmp dir
	if ($totlRefs > 5){
		$algCmd .= "\n\n#copying read files to local tmp\n";
		for (my $ii=0;$ii<@pa1;$ii++){
			$algCmd .= "cp $pa1[$ii] $nodeTmp\n";$pa1[$ii]=~ m/\/([^\/]+$)/;$pa1[$ii] = $nodeTmp."/$1";
		}
		for (my $ii=0;$ii<@pa2;$ii++){
			$algCmd .= "cp $pa2[$ii] $nodeTmp\n";$pa2[$ii]=~ m/\/([^\/]+$)/;$pa2[$ii] = $nodeTmp."/$1";
		}
		for (my $ii=0;$ii<@paS;$ii++){
			$algCmd .= "cp $paS[$ii] $nodeTmp\n";$paS[$ii]=~ m/\/([^\/]+$)/;$paS[$ii] = $nodeTmp."/$1";
		}
	}
	
#	if ($numLib==0){
#		$numLib = scalar @paS; 
#	}
	#die "@bwtIdxs\n"; 
	#organize mapping against 1 or more refs
	for (my $kk=0;$kk<@bwtIdxs;$kk++){  #iterator over different genomes
		#die  "$finalDS[$kk]/$outNms[$kk]-smd.bam\nXXX\n$mappDir[$kk]\n";
		if ($decoyModeActive || $map2ndTogether>0){
			my $allDone=1;
			for (my $k=0;$k<@outNms;$k++){
				#print "$finalDS[$k] $outNms[$k] $mappDir[$k]\n";
				$allDone =0 unless(check_map_done($doCram, $finalDS[$k], $outNms[$k], $mappDir[$kk]));
			}
			#die;
			last if ($allDone);
		} elsif (-e $tmpOut22[$kk] || check_map_done($doCram, $finalDS[$kk], $outNms[$kk], $mappDir[$kk])){ #this check is for non-decoy mode
			#(-e "$finalDS[$kk]/$outNms[$kk]-smd.bam" && -e "$finalDS[$kk]/$outNms[$kk]-smd.bam.coverage.gz") ){
			next;
		} #|| -e "$mappDir[$kk]/$outNms[$kk]-smd.bam.coverage.gz"
		$cntAli++;
		$bamFresh=1;
		#print "ali\n";
		$algCmd .= "#--------------- $kk ---------------\nrm -rf $mappDir[$kk]\nmkdir -p $mappDir[$kk];\n"; 
		$algCmd .= "rm -rf $unaligned $tmpUna;\nmkdir -p $unaligned $tmpUna;\n" if ($unaligned ne "");
		
		my $algCmdBase = getAlgnCmdBase($mapperProgLoc,$NcoreL, $readTec,$MFopt{mapModeTogether},$unaligned);
#die "algCmdBase::$algCmdBase\n";
		my @subBams; 
		my @accR1=(); my @accR2=(); my @accRS=();
		for (my $i=0; $i< $numLib; $i++){
			my $usePairs=1;
			$anyUsedPairs =1 if ($usePairs);
			if ($i >= scalar @pa1){$usePairs=0;}
			my $iS = $i - scalar @pa1;
			
			#test if more reads can be accummulated in the next round...
			#print "@libsOri  : $i $numLib ". scalar @pa1 ."\n";
			if ($mapperProgLoc==1 ){
				#accumulate reads..
				if ($usePairs){ push(@accR1,$pa1[$i]); push(@accR2, $pa2[$i]); 
				} else { push(@accRS,$paS[$iS]); }
				if ( ($i+1) < $numLib && ((($i+1) >= scalar @pa1 && $usePairs==0) || (($i+1) < scalar(@pa1) && $usePairs==1))
						&& $libsOri[$i] eq $libsOri[$i+1]){ #same reads, same lib in next round, all set!
					next;
				}
			} else { #other mappers can't use multiple input files..
				if ($usePairs){ @accR1 = ($pa1[$i]); @accR2 = ($pa2[$i]); 
				} else { @accRS  = ($paS[$iS]); }
			}
			
			#$pa1[$i] =~ m/\/([^\/]+)\.f.*q$/;
			my $rgID = "$outName";
			my $rgStr = getRgStr($outName,$libsOri[$i],$libsOri[$iS],$usePairs,$mapperProgLoc);
			#die "$rgStr\n";
			if ($mapperProgLoc==1){ #bowtie2
				if ($usePairs){
					$algCmd .= "$algCmdBase -x $bwtIdxs[$kk] -1 ".join(",",@accR1) ." -2 ".join(",",@accR2);
				} else {
					$algCmd .="$algCmdBase -x $bwtIdxs[$kk] -U ".join(",",@accRS);#$paS[$iS];
				}
				$algCmd .= " --rg-id $rgID $rgStr ";
			} elsif ($mapperProgLoc==2){ #bwa
				die "single end mapping not implemented for bwa\n" if (!$usePairs);
				$algCmd .= $algCmdBase." -R $rgStr $REF " . join(",",@accR1). " " . join(",",@accR2); ##$pa1[$i]." ".$pa2[$i];
			} elsif ($mapperProgLoc==3){ #minimap2
				$algCmd .= $algCmdBase." -R $rgStr -a $REF$MFcontstants{mini2IdxFileSuffix} ";
				if ($usePairs){ $algCmd .= join(",",@accR1) . " " . join(",",@accR2); #$pa1[$i]." ".$pa2[$i]
				} else { $algCmd .= join(" ",@accRS)
				}
			}elsif ($mapperProgLoc==4){ #kma
				$algCmd .= "$algCmdBase -t_db $REF$MFcontstants{kmaIdxFileSuffix}  ";
				if ($usePairs) {$algCmd .= " -ipe " . join(",",@accR1). " " . join(",",@accR2) ;}# $pa1[$i] $pa2[$i] " 
				if ($iS >= (scalar @pa1)){$algCmd .= " -i $paS[$iS] " ;}
				$algCmd .= " -o $tmpOutxtra[$i] ";
			}elsif ($mapperProgLoc==5){ #strobealign
				$algCmd .= "$algCmdBase --rg=$rgStr $REF  ";
				if ($usePairs) {$algCmd .=  join(",",@accR1). " " . join(",",@accR2) ;}# $pa1[$i] $pa2[$i] " 
				if ($iS >= (scalar @pa1)){$algCmd .= " $paS[$iS] " ;}
			}
			@accR1=(); @accR2=(); @accRS=(); #and empty what was already mapped

			#samtools filter and other postfiltering..
			my %postTreat = (MapperProg=>$mapperProgLoc, readTec=>$readTec, NcoreL => $NcoreL,nodeTmp =>$nodeTmp,tmpOut21=>$tmpOut21, 
										subBamsAR => \@subBams,unaligned =>$unaligned,baseN => $baseN,xtraSamSteps1 => $xtraSamSteps1,
										decoyModeActive => $decoyModeActive, map2ndTogether => $map2ndTogether,regsAR => \@regs, 
										doCram => $doCram, finalDSar => \@finalDS, outNmsAR => \@outNms, mappDirAR => \@mappDir, 
										reg_lcsAR => \@reg_lcs);
			
			my ($algCmdPost,$subBamAR) = alignPostTreat(\%postTreat, $i, $kk);
			$algCmd .= $algCmdPost;
			@subBams = @{$subBamAR};
			
		}

		my $filterStep="";

		#bottleneck step: bam filtering
		#my $bamfilter = getProgPaths("bamFilter_scr");
		#my $filterStep="";
		#if ( ($mapperProgLoc==3 || $mapperProgLoc==5 ) && ($readTec eq "ONT" || $readTec eq "PB")){ #low id long reads...
		#	$filterStep = " | $smtBin view -@ $NcoreL - | $bamfilter $MFopt{bamfilterPB} | $smtBin view -@ $NcoreL -b1 - ";
		#} else { #illumina parameters
	#		$filterStep = " | $smtBin view -@ $NcoreL - | $bamfilter $MFopt{bamfilterIll} | $smtBin view -@ $NcoreL -b1 - ";
	#	}

		for (my $k=0;$k<@subBams;$k++){
			#$tmpOut22 = $tmpOut."/$outNms[$k].iniAlignment.bam";
			#print $tmpOut22."\n";
			#if (-e $tmpOut22){ next;}
			my $tarBam = $tmpOut22[$kk];
			$tarBam = $tmpOut22[$k] if ($decoyModeActive || $map2ndTogether>0);
			next if ($subBams[$k] eq "");#case that decoy map has already created parts of the mappings..
			if (1){   #always active #($numLib > 1){
				$algCmd .= "\n$smtBin cat ".$subBams[$k]." $filterStep > $tarBam\n"; #$k here, because this refers to @regs
				$algCmd .= "\nrm -f ".$subBams[$k]."\n";
			} else { #nothing else, compression is now in the bowtie2 phase
				#$algCmd .= "\n$smtBin view -bS -F 4 -@ $Ncore $subBams[0] > $tmpOut22\n";
				
				$algCmd .=  "\nmv $subBams[$k] $tarBam\n"; #kk here, because this is the alignment to each single bwtIdx
			}
		}

	#die "$algCmd\n";
	}
	#die "$algCmd\n";
	my $mapProgNm = getMapProgNm($mapperProgLoc);

	foreach my $mpdSS (@mappDir) {
		if (!-e "$mpdSS/map.sto"){$bamFresh=1;}
		#$algCmd .= "echo '$mapProgNm' > $mpdSS/map.sto\n";
	}

	#system($algCmd." -U $p2 -S $tmpOut22\n");
	#die "$algCmd\n";
	my $nodeCln = "\nrm -rf $nodeTmp\n";
	#print "$cntAli new alignments\n" if ($cntAli > 0);
#	die $logDir.$bashN."bwtMap.sh";
	if ( $bamFresh  ){
		#die "$unzipcmd.$algCmd.$nodeCln";
		my $baseMapHDD = $HDDspace{mapping} ;  $baseMapHDD =~ s/G$//;
		my $preHDDspace=$QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = int((200+$inputFileSizeMB{$outName})*$baseMapHDD/1024)+15  ."G";
		#print "$QSBoptHR->{tmpSpace}\n";
		$jobN = "_MAP$JNUM$supTag.$outNms[0]"; $bamFresh = 1; 
		($jobN,$tmpCmd) = qsubSystem($qdir.$bashN."map$supTag.sh",
				$unzipcmd.$algCmd.$nodeCln,  #$unalignCmd
				$Ncore,int($MFopt{MapperMemory}/$Ncore+1)."G",$jobN,$jDepe,"",$immediateSubm,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		$retCmds .= $tmpCmd;
		$QSBoptHR->{tmpSpace}=$preHDDspace;
	#die "\nBAM  $bamFresh\n$immediateSubm\n";
	}

	#die "TOOfar\n";
	#xtraSamSteps1 = isSorted
	$params{sortedbam}=$isSorted; $params{bamIsNew} = $bamFresh; $params{is2ndMap} = $is2ndMap;
	$params{immediateSubm} = $immediateSubm; $params{usePairs} = $anyUsedPairs;
	#die "$jobN\n$immediateSubm\n";
	return($jobN,$retCmds,\%params);
}	
	










sub bamDepth{
	my ($dirsHr, $jDep,$mapparhr) = @_;
	#die "bamdep\n";
	my $readCov_Bin = getProgPaths("readCov");
	my $outName = $dirsHr->{smplName};
	my $doCram =  $dirsHr->{cramAlig};
	my $mappDir = ${$dirsHr}{glbMapDir};	my $nodeTmp = ${$dirsHr}{nodeTmp}."_bamDep/$outName/";
	my $tmpOut = ${$dirsHr}{glbTmp};	my $finalD = ${$dirsHr}{outDir};
	my $qdir = $logDir; $qdir = ${$dirsHr}{qsubDir} if (exists( ${$dirsHr}{qsubDir} ));
	my $supportRds = $dirsHr->{mapSupport};
	my $supTag = ""; if ($supportRds){$supTag = ".sup";}
	my $REF = $dirsHr->{sbj}; #target to map onto, can by ","-spearated list	
	#my $bedCovBin = getProgPaths("bedCov");#"/g/bork5/hildebra/bin/bedtools2-2.21.0/bin/genomeCoverageBed";
	#my $mosDepBin = getProgPaths("mosdepth");
	my $allowDeleteMap = 1;

	my %params = %{$mapparhr};
	my ($isSorted , $bamFresh,$is2ndMap, $usePairs) =($params{sortedbam},$params{bamIsNew},$params{is2ndMap},$params{usePairs});
	#recreate base pars from map2tar sub
	my $locDoRmDup = $MFopt{MapperRmDup};
	$locDoRmDup = 0 if (!$usePairs);
	my $immediateSubm = $params{immediateSubm} ;
	#print "pairs used:$usePairs\n";

	my $baseN = "$outName$supTag";
	my $bashN = "";	if ($is2ndMap){$bashN = "$outName"; $bashN =~ s/,/./;}
	my $mappingRes = $tmpOut."/$baseN.iniAlignment.bam"; #this is the input, result of previous mapping steps
	my $cramSTO = "$mappDir/$baseN-smd.cram.sto";
	my $nxtBAM = "$tmpOut/$baseN-smd.bam";
	my $sortTMP = $nodeTmp."/$baseN.srt";
	my $sortTMP2 = $nodeTmp."/$baseN.2.srt";
 
	#check if already done
	my $outstat = check_map_done($doCram, $finalD, $baseN, $mappDir);
	my $outstat2 = check_depth_done($doCram, $finalD, $baseN, $mappDir);
	#die "$outstat $mappingRes $nxtBAM\n";
	if ($outstat2){return ("","",$outstat);}

	#die "$mappDir/$baseN-smd.bam.coverage.gz\n$finalD/$baseN-smd.bam.coverage.gz\n" ;#if (-e "$finalD/$baseN-smd.bam.coverage.gz");
		
	my $numCore = ${$dirsHr}{sortCores};#new functionality with sambamba
	#my $biobambamBin = "bamsormadup";
	#depth profile
	#my $cmd = "$smtBin faidx $REF\n $smtBin view -btS $REF.fai $tmpOut/$baseN.sam > $tmpOut/$baseN.bam\n";
	my $cmd = "sleep 1\n";
	#$cmd .= "mkdir -p $nodeTmp\n";
	$cmd .= "mkdir -p $mappDir $tmpOut $nodeTmp\n";
	#sort & remove duplicates
	if ($outstat && !$outstat2 ){
		if ($doCram){
			$cmd .= "#uncramming already stored results..\n" . cram2bsam("$finalD/$baseN-smd.cram",$REF,$mappingRes,1,$numCore) ."\n" ;
		} elsif (-e "$finalD/$baseN-smd.bam") {
			$mappingRes = "$finalD/$baseN-smd.bam";
			$allowDeleteMap = 0;
		} else {
			#might require extra check on tmp dir, if map were not moved yet..
			die "bamDepth:: not sure what to do with input..\n";
		}
	}

	if ($locDoRmDup){
	#bit complicated in samtools now: first sort by name, then fixmate (samtools 1.8)
	#".int($MFopt{mapSortMemGb}/$numCore)."G
		$cmd .= "echo \"Sorting .bam and removing duplicates...\"\n";
		$cmd .= "$smtBin sort -n -m 768M -T $sortTMP -@ $numCore $mappingRes | $smtBin fixmate -m -@ $numCore - - | $smtBin sort -T $sortTMP2 -m 768M -@ $numCore - | $smtBin markdup -s -r -@ $numCore - $nxtBAM\n";
	} else {
		if ($isSorted ){
			$cmd .= "mv $mappingRes $nxtBAM\n";
		} else {
			$cmd .= "echo \"Sorting .bam ...\"\n";
			$cmd .= "$smtBin sort -@ $numCore -m 768M -T $sortTMP $mappingRes  > $nxtBAM\n";
		}
	}
	$cmd .= "echo \"Building .bam index...\"\n";
	$cmd .= "$smtBin index -@ $numCore $nxtBAM\n";
	$cmd .= "[ -s $nxtBAM ] || exit 3\n";#check if output is empty (samtools crashed?)
	$cmd .= "rm -f $mappingRes \n" if (!$isSorted && $allowDeleteMap);

	#$cmd .= "$smtBin view -bS $tmpOut21 > $tmpOut/$baseN.bam\n";
	#my $mdJar = "/g/bork5/hildebra/bin/picard-tools-1.119/MarkDuplicates.jar";
	#$cmd .= "java -Xms1g -Xmx24g -XX:ParallelGCThreads=2 -XX:MaxPermSize=1g -XX:+CMSClassUnloadingEnabled ";		$cmd .= "-jar $mdJar INPUT=$mappDir/$baseN-s.bam OUTPUT=$mappDir/$baseN-smd.bam METRICS_FILE=$mappDir/$baseN-smd.metrics ";		
	#$cmd .= "AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE\n\n";
	
	#$cmd .= "$smtBin index $nxtBAM\n";#$mappDir/$baseN-smd.bam\n";
	#takes too much space, calc on the fly!
	#$cmd .= "$smtBin depth $mappDir/$baseN-s.bam > $mappDir/$baseN.depth.bam.txt\n";#depth file (different estimator)
	my $covCmd = "";
	if (!-e "$nxtBAM.coverage.gz"){
		$covCmd .= "echo \"Creating coverage of reads...\"\n";
		#$covCmd .= "$mosDepBin -t 4 -x $nxtBAM.coverage $nxtBAM\n";
		#$covCmd .= "mv $nxtBAM.coverage.per-base.bed.gz $nxtBAM.coverage.gz\n";
		#$covCmd = "$bedCovBin -ibam $nxtBAM -bg > $nxtBAM.coverage\n";
		#$covCmd .= "rm -f $nxtBAM.coverage.gz\n$pigzBin -f -p $numCore $nxtBAM.coverage\n";
		my $sam2bed = getProgPaths("samcov2bed");
		$covCmd .= "$smtBin depth -@ $numCore $nxtBAM | $sam2bed | $pigzBin -p $numCore -c > $nxtBAM.coverage.gz\n";
		#$covCmd .= "rm -f $nxtBAM.coverage.gz\n$pigzBin -f -p $numCore $nxtBAM.coverage\n";
		if ($map2ndMpde == 3){ #map unmapped to genecat
			$covCmd .= "$readCov_Bin $nxtBAM.coverage.gz - 100\n"; #read length doesn't matter, since no win is given out and no gene is present
		}
	}
	#die "$covCmd\n";
	
	#jgi depth profile - not required, different call to metabat
	$covCmd .= jgi_depth_cmd([$nxtBAM],$nxtBAM,95) if ($MFopt{DoJGIcoverage});
	#$covCmd .= "/g/bork5/hildebra/bin/bedtools2-2.21.0/bin/genomeCoverageBed -ibam $nxtBAM -bg  | ".'awk \'BEGIN {pc=""} {	c=$1;	if (c == pc) {		cov=cov+$2*$5;	} else {		print pc,cov;		cov=$2*$5;	pc=c}';
	#$covCmd .= "} END {print pc,cov}\' $nxtBAM.coverage | tail -n +2 > $nxtBAM.coverage.percontig";
	$covCmd .= "rm -f $nxtBAM.bai;\n";
	my ($CRAMcmd,$CRAMf) = bam2cram($nxtBAM,$REF,1,$doCram,$cramSTO, $numCore);
	$CRAMcmd = "echo \"Building .cram ...\"\n$CRAMcmd" if ($CRAMcmd ne "");
	unless ($mappDir eq $tmpOut){
		$CRAMcmd .= "\nmv $CRAMf* $mappDir 2>/dev/null \n" unless ($CRAMf eq "");
		$CRAMcmd .= "mv $nxtBAM* $mappDir 2>/dev/null \n" unless ($nxtBAM eq "");
	}
	$CRAMcmd .= "echo \"".basename($nxtBAM)."\" > $mappDir/done.sto\n" unless($supportRds);


	my $newJobN = "SBAM$supTag$JNUM$outName";	
	#subsequent jobs not dependent on this one
	my $jobN2 = $jDep; my $retCmds="";
	$covCmd = "" if ( -s "$nxtBAM.coverage");
	
	#my $cleaner = "mv $tmpD $fin
	if (0){
		print "B1 " if ( $doCram && !-e $cramSTO); print "B2 " if (  !-s "$nxtBAM.coverage"); print "B3 " if ( $bamFresh);
		print " ".$cramSTO."\n";
	}
	my $nodeCln = "\nrm -rf $nodeTmp;\necho \"DONE map2\"\n";
	#die "$cmd\n$covCmd\n$CRAMcmd\n";
	if ( ($doCram && !-e $cramSTO) || ( !-s "$nxtBAM.coverage.gz") ){#|| $bamFresh){
		my $baseMem = 1; $baseMem=20 if ($MFopt{largeMapperDB});
		my $preHDDspace=$QSBoptHR->{tmpSpace};		my $baseMapHDD = $HDDspace{mapping} ;  $baseMapHDD =~ s/G$//;
		$QSBoptHR->{tmpSpace} = int($inputFileSizeMB{$curSmpl}*$baseMapHDD/1024)+15  ."G";		if (${$dirsHr}{submit}){
			
		($jobN2,$retCmds) = qsubSystem($qdir.$bashN."map2$supTag.sh",
				$cmd."\n".$covCmd."\n".$CRAMcmd."\n$nodeCln\n"#.$covCmd2
				,$numCore,  int($baseMem +$MFopt{mapSortMemGb}/$numCore)."G",$newJobN,$jDep,"",$immediateSubm,$QSBoptHR->{General_Hosts},$QSBoptHR);
		$QSBoptHR->{tmpSpace}=$preHDDspace;
		} else {
			$cmd =~ s/sleep \d+//;
			$retCmds = $cmd."\n".$covCmd."\n".$CRAMcmd."\n$nodeCln\n";
		}
	} 
	#die();
	return($jobN2,$retCmds,$outstat);
}

sub mergeMP2Table($){
	my ($dir_MP2) = @_;
	return if (!$MFopt{DoMetaPhlan});
	if ($progStats{metaPhl2FailCnts}){
		print "$progStats{metaPhl2FailCnts} / ". ($progStats{metaPhl2ComplCnts}+$progStats{metaPhl2FailCnts}) ." samples with incomplete Metaphlan assignments\n";
		return;
	}
	my $outD = $dir_MP2;
	$outD =~ s/[^\/]+\/?$//;

	my $mergeTblScript = getProgPaths("metPhl2Merge");#"/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
	my $MP2Mstone = "$outD/MP2.cnt.stone";
	my $prevCnts = 0;
	if (-e $MP2Mstone){
		$prevCnts = `cat $MP2Mstone`; chomp $prevCnts; 
	}
	my $redoMP2Tables = 0;
	print "\nAll samples ($progStats{metaPhl2ComplCnts}) have metaphlan assignments.\n";
	if ($prevCnts < $progStats{metaPhl2ComplCnts}){
		print "Redoing metaphlan2 table merge ,due to higher number of samples detected ($progStats{metaPhl2ComplCnts}, prev: $prevCnts)\n";
		$redoMP2Tables = 1;
	} else {
		return;
	}
	my $getHDerCmd = "head -n1 ";
	if ($MFopt{DoMetaPhlan3}){
		$getHDerCmd = "head -n2 ";
	}
	my @slvl = ("k__","p__","c__","o__","f__","g__","s__"); my @llvl = ("kingdom","phylum","class","order","family","genus","species");
	my $mrgCmd = ""; 
	my $tmprawF = " $outD/MePh.Raw.mat";
	my $premrgCmd = "$mergeTblScript $dir_MP2/*MP2.txt > $tmprawF\n" ;
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$outD/MePh.all.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$getHDerCmd  $tmprawF | tail -n1 > $outF\n"; #needed for version 3
			$mrgCmd .= "cat  $tmprawF | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g' >> $outF\n" ;
		}
	}
	
	$mrgCmd .= "echo \"$progStats{metaPhl2ComplCnts}\" > $MP2Mstone";

	#systemW "$mrgCmd";
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
	my ($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/MP2merg.sh",$mrgCmd,1,"80G","MP2mrg","","",1,[],$QSBoptHR) ;
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	return;
	
	
	#this part is no longer used since MATAFILER v 0.23
	#die "$mrgCmd\n";

	$premrgCmd = "$mergeTblScript $dir_MP2/*MP2.noV.noB.txt > $dir_MP2/All.MP2.noV.noB.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$dir_MP2/MP2.noV.noB.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$getHDerCmd  $dir_MP2/tmp.All.MP2.mat > $outF\n";
			$mrgCmd .= "cat $dir_MP2/All.MP2.noV.noB.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g' >> $outF\n"  unless (-e "$outF"&& !$redoMP2Tables);
		}
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.noV.noB.mat\n\n";
	$premrgCmd = "$mergeTblScript $dir_MP2/*MP2.noV.txt > $dir_MP2/All.MP2.noV.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$dir_MP2/MP2.noV.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$$getHDerCmd  $dir_MP2/tmp.All.MP2.mat > $outF\n";
			$mrgCmd .= "cat $dir_MP2/All.MP2.noV.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g'  >> $outF\n"  unless (-e "$outF"&& !$redoMP2Tables);
		}
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.noV.mat\n\n";
	#viruses
	$premrgCmd = "$mergeTblScript $dir_MP2/*MP2.VirusOnly.txt > $dir_MP2/All.MP2.VirusOnly.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$dir_MP2/MP2.VirusOnly.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$$getHDerCmd  $dir_MP2/tmp.All.MP2.mat > $outF\n";
			$mrgCmd .= "cat $dir_MP2/All.MP2.VirusOnly.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g'  >> $outF\n"  unless (-e "$outF"&& !$redoMP2Tables);
		}
	}
	$mrgCmd .= "rm -f $dir_MP2/tmp.All.MP2.mat $dir_MP2/All.MP2.VirusOnly.mat\n\n";

	#die $mrgCmd."\n";
}

sub prepMetaphlan{
	
	return unless ($MFopt{DoMetaPhlan3} || $MFopt{DoMetaPhlan});
	return; #takes too much time..
	my $met3Bin = getProgPaths("metPhl2");
	print "Checking metaphlan version .. ";
	my $vstr = "";
	$vstr = `$met3Bin --version 2>/dev/null`;
	$vstr =~ m/version ([\.\d]+)/;
	print "$1 ";
	if ($MFopt{DoMetaPhlan3}){
		$MFopt{DoMetaPhlan} = 1;
		if (substr($1,0,1) < 3.0){ print "Metaphlan version below 3, reinstall updated metaphlan3\n"; exit(33);}
		print "will use version 3\n";
	}
	#set globally to version 3
	$MFopt{DoMetaPhlan3} = 1 if (substr($1,0,1) >= 3.0);
}


sub metphlanMapping{
	my ($cleanSeqSetHR,$tmpD,$finOutD,$smp,$Ncore,$deps) = @_;
	
	my $inF1a = ${$cleanSeqSetHR}{arp1}; my $inF2a = ${$cleanSeqSetHR}{arp2}; my $inFSa = ${$cleanSeqSetHR}{singAr}; my $mergRdsHr = ${$cleanSeqSetHR}{mrgHshHR};

	
	my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
	
	my $metPhl2Bin = getProgPaths("metPhl2");#"/g/bork3/home/hildebra/bin/metaphlan2/metaphlan2.py";
	#path to metaphlan DB
	my $mpDB = getProgPaths("metPhl2_db",0);#metaphlan2/db_v20/mpa_v20_m200
	my $metPhl2Merge = getProgPaths("metPhl2Merge");#"/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
	
	#split metaphlan DB path for v3
	$mpDB =~ m/^(.*)\/([^\/]+)$/;
	my $mpDB1 = $1."/";my $mpDB2 = $2;
	
	my $stone = $finOutD."$smp.MP2.sto";
	if (!$MFopt{DoMetaPhlan} ||  -e $stone){return;}
	my @car1 = @{$inF1a}; my @car2 = @{$inF2a}; my @sar = @{$inFSa};
	my $inF1 = join(",",@car1); my $inF2 = join(",",@car2); my $inFS = join(",",@sar); 
	system "mkdir -p $finOutD\n" unless (-d $finOutD);
	my $finOut = $finOutD."$smp.MP2.txt";
	my $finOut_noV = $finOutD."$smp.MP2.noV.txt";
	my $finOut_noVB = $finOutD."$smp.MP2.noV.noB.txt";
	my $finOut_Vo = $finOutD."$smp.MP2.VirusOnly.txt";
	my $sam = "$tmpD/metph2.sam";
	my $v3params = "";$v3params = " --sample_id $smp --nproc $Ncore --unknown_estimation --nreads \$readN -o " if ($MFopt{DoMetaPhlan3}); #--unknown_estimation -> requires --nreads
	my $v2params = ""; $v2params = "--ignore_viruses >" if (!$MFopt{DoMetaPhlan3});
	my $taxinfo = "--mpa_pkl $mpDB.pkl ";
	$taxinfo = "--bowtie2db $mpDB1 -x $mpDB2 " if ($v3params ne "");
	my $qsubFile = $logDir."metaPhl2.sh";
	
	
	my $cmd = "mkdir -p $tmpD\n$bwt2Bin --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S $sam -p $Ncore -x $mpDB ";
	$cmd .= "-1 $inF1 -2 $inF2 " if (@car1 > 0 );;
	$cmd .= "-U $inFS " if (@sar > 0);
	$cmd .= "\n\n";
	$cmd .=  "sleep 1\nreadN=\$(grep -v 'Warning:' $qsubFile.etxt |  head -n1 | cut -f1 -d' ')\necho \$readN\n";
	$cmd .= "$metPhl2Bin $sam --input_type sam $taxinfo $v2params $v3params $finOut\n";
	#$cmd .= "$metPhl2Bin $sam --input_type sam $taxinfo $v2params $v3params $finOut_noV\n" if (!$MFopt{DoMetaPhlan3});
	#$cmd .= "$metPhl2Bin $sam --input_type sam --ignore_bacteria --ignore_archaea $taxinfo $v2params $v3params $finOut_noVB\n";
	#$cmd .= "$metPhl2Bin $sam --input_type sam --ignore_bacteria --ignore_eukaryotes --ignore_archaea $taxinfo $v2params $v3params $finOut_Vo\n";
	$cmd .= "rm -f $sam\n";
	my $mergeStr = "$metPhl2Merge *.MP2.txt > $finOutD/comb.MP2.txt";
	$cmd .= "echo \' $mergeStr \' > $stone\n";
	my $jobN = "MP2$JNUM";
	
	my ($jobN2,$tmpCmd) = qsubSystem($qsubFile,
			$cmd,$Ncore,"3G",$jobN,$deps,"",1,[],$QSBoptHR);
	$jobN  = $jobN2;
	return $jobN;
}

sub TaxaTarget{
	my ($cleanSeqSetHR,$tmpD,$finOutD,$smp,$Ncore,$deps) = @_;
		my $inF1a = ${$cleanSeqSetHR}{arp1}; my $inF2a = ${$cleanSeqSetHR}{arp2}; my $inFSa = ${$cleanSeqSetHR}{singAr}; #my $mergRdsHr = ${$cleanSeqSetHR}{mrgHshHR};

	my $stone = $finOutD."$smp.TaxTar.sto";
	if (!$MFopt{DoTaxaTarget} ||  -e $stone){return;}

	my @car1 = @{$inF1a}; my @car2 = @{$inF2a}; my @sar = @{$inFSa};	
	my $inF1 = join(",",@car1); my $inF2 = join(",",@car2); my $inFS = join(",",@sar); 
	my $taxTBin = getProgPaths("TaxaTarget");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
	my $taxTarDir = $taxTBin;
	$taxTarDir =~ s/run_pipeline_scripts\/run_protist_pipeline_fda.py//;
	$taxTarDir =~ s/python //;
	my $tar1 = $car1[0]; my $tar2 = $car2[0];
	my $cmd = "";
	#if ($tar =~ m/\.gz$/){		$tar1 =~ s/\.gz$//;		$cmd .= "gunzip -c $tar > $tar1\n";	}
	#$tmpD 
	$cmd .= "$taxTBin -r $tar1 -r2 $tar2 -e $taxTarDir/run_pipeline_scripts/environment.txt --tmp -t $Ncore -o $taxTarDir\n";
	my $jobN = "TT$JNUM";

	my ($jobN2,$tmpCmd) = qsubSystem($logDir."taxtar.sh",
			$cmd,$Ncore,"3G",$jobN,$deps,"",1,[],$QSBoptHR);
			
	#die $logDir;
	return $jobN2;
}



sub clean_tmp{#routine moves output from temp dirs to final dirs (that are IO limited, thus single process)
	my ($clDar,$cpref,$cpnodel,$jDepe,$tag,$curJname) = @_;
	my @clDa = @{$clDar};
	my @cps = @{$cpref};
	my @cpsND = @{$cpnodel};
	#die;
	#die "@cps ==0 && @clDa ==0 && @cpsND\n";
	if (@cps ==0 && @clDa ==0 && @cpsND==0){return $jDepe;}
	my $cmd = "";
	for (my $i=0;$i<@cps;$i+=2){
		next if (length($cps[$i+1] ) < 5);
		$cmd .= "rm -r -f $cps[$i+1]\nmkdir -p $cps[$i+1]\nrsync -r --remove-source-files $cps[$i] $cps[$i+1]\n" if (@cps > 1);
	}
	for (my $i=0;$i<@cpsND;$i+=2){
		next if (length($cpsND[$i+1] ) < 5);
		$cmd .= "mkdir -p $cpsND[$i+1]\nrsync -r  --remove-source-files $cpsND[$i] $cpsND[$i+1]\n"  if (@cpsND > 1);
	}
	$cmd .= "rm -f -r ".join(" ",@clDa)."\n" if (@clDa > 0  );
	$cmd .= "echo $tag >> $collectFinished\n" unless ($tag eq "");
	$cmd .= "echo \"Done cleaning tmp dir\"\n";
	
	#die "@clDa\n\n$cmd\n";
	my $clnBash = $logDir."clean$curJname.sh";
	if ($curJname eq ""){
		$curJname = "_cln$JNUM";
	}
	#add extra job dependency on d2star
	my $xtraJDep = "";
	#$xtraJDep = $QSBoptHR->{rTag}."_d2met";
	my ($jdep,$tmpCmd) = ("","");
	#die "$jDepe";
	if ($jDepe =~ m/^[\s;]+$/){
		#print "Excecuting copy/rm locally..\n$cmd\n";
		
		systemW $cmd;#." \$" if ($cmd ne "");
	} else {
	#die "X${jDepe}X\n";
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($jdep,$tmpCmd) = qsubSystem($clnBash,$cmd,1,"15G",$curJname,$jDepe.";$xtraJDep","",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}
	#die "jdeps: $jDepe\n$cmd\n";
	return ($jdep );
}



	
sub manageFiles{
	my ($cAssGrp, $cMapGrp, $rmRdsFlag, $doPreAssmFlag, $curOutDir , $jdep, $smplTmpDir,$AssemblyGo) = @_;
	#cleaning of filtered reads & copying of assembly / mapping files
	#$AsGrps{$cAssGrp}{ClSeqsRm} .= $smplTmpDir."seqClean/".";";
	
	my @copiesNoDels=();
	push(@copiesNoDels , @{$AsGrps{$cAssGrp}{PsAssCopies}}) if (exists($AsGrps{$cAssGrp}{PsAssCopies}));
	push(@copiesNoDels , @{$AsGrps{$cAssGrp}{MapSupCopies}} ) if (exists($AsGrps{$cAssGrp}{MapSupCopies}));
	push(@copiesNoDels , @{$AsGrps{$cAssGrp}{MapCopiesNoDel}} ) if (exists($AsGrps{$cAssGrp}{MapCopiesNoDel}));
	my @cleans = ();
	#if ($MappingGo && @bwt2outD==0){ #sync clean of tmp with scnd mapping..
		#push(@cleans , $smplTmpDir);# $smplTmpDir,$metagAssDir."seqClean/filter*",
	#}
	
	#all dependencies before deleting tmp dirs
	#die "$AsGrps{$cAssGrp}{SeqClnDeps}\n";
	my $totJdeps = $jdep . ";" . "$AsGrps{$cAssGrp}{SeqClnDeps}" . ";" . $AsGrps{$cAssGrp}{MapDeps} . ";". $AsGrps{$cAssGrp}{scndMapping}.";".$AsGrps{$cAssGrp}{readDeps}.";".$AsGrps{$cAssGrp}{prodRun};
	
#	for ( ($jdep , $AsGrps{$cAssGrp}{MapDeps} , $AsGrps{$cAssGrp}{scndMapping},$AsGrps{$cAssGrp}{prodRun}) ){
#		push(@sampleDeps, $_ ) if (defined $_ && $_ ne "");
#	}
	#die "$totJdeps\n";
	#print $AsGrps{$cAssGrp}{ClSeqsRm}." FF ".$AsGrps{$cMapGrp}{ClSeqsRm}."\n";
	if ($rmRdsFlag && exists($AsGrps{$cAssGrp}{ClSeqsRm})){push(@cleans,split(/;/,$AsGrps{$cAssGrp}{ClSeqsRm} )); $AsGrps{$cMapGrp}{ClSeqsRm} = "";}
	if ($MFconfig{rmScratchTmp} && $totJdeps =~ m/^[\s;]*$/ && $rmRdsFlag){push(@cleans,$smplTmpDir);}
	#print "test\n\n" if ($totJdeps =~ m/^[\s;]*$/);
	
	#die "\n@cleans\n$totJdeps\n";
	
	my @moves = (@{$AsGrps{$cAssGrp}{MapCopies}},@{$AsGrps{$cAssGrp}{AssCopies}});
	#die "@copies\n";
	#print "@cleans\n@moves\n@copiesNoDels\n";
	my $cln1="";
	#if ($MFopt{DoAssembly} && $boolAssemblyOK || (!$MFopt{DoAssembly}) ){
	if (!$MFconfig{remove_reads_tmpDir} || $doPreAssmFlag == 1){@cleans =();}
	my $finishedClnDir = $curOutDir;  $finishedClnDir = "" if ($doPreAssmFlag);
	#die "XXYZ\n@cleans\n @moves\n @copiesNoDels\n";
	$cln1 = clean_tmp(\@cleans,\@moves, \@copiesNoDels,$totJdeps,$finishedClnDir,"");#$AsGrps{$cAssGrp}{CSfinJobName}); #.";".$contRun
	
	$AsGrps{$cAssGrp}{BinDeps} .= ";$cln1" if ($AsGrps{$cAssGrp}{MapDeps} ne "");
	$QSBoptHR->{LocationCheckStrg}="";
	#clean up assembly groups
	$AsGrps{$cAssGrp}{ClSeqsRm} = ""; @{$AsGrps{$cAssGrp}{MapCopies}} = ();
	@{$AsGrps{$cAssGrp}{MapCopiesNoDel}} = (); @{$AsGrps{$cAssGrp}{MapSupCopies}} = ();
	$AsGrps{$cAssGrp}{MapDeps} = "" if ($AssemblyGo);  
	$AsGrps{$cAssGrp}{readDeps} = ""; $AsGrps{$cAssGrp}{DiamDeps} = "";
	#
	return $cln1;
}



sub remComma($){
	my ($in) = @_;
	$in =~ s/,//g;
	return $in;
}


sub sdmStats($$){
	my ($inF,$inD) = @_;
	my $MaxLengthHistBased=0;
	my $filStats = getFileStr("$inD/LOGandSUB/sdm/filter_lenHist.txt",0);
	if ($filStats ne ""){
		my @tmpSpl = split(/\n/,$filStats);
		if (@tmpSpl != 0){$tmpSpl[$#tmpSpl]=~m/^(\d+)\s/; $MaxLengthHistBased= $1;}
	}
	#tmp deactivate(might still be a bug in some old sdm version):
	#$MaxLengthHistBased=0;
	#die "$MaxLengthHistBased\n";
	$filStats = "";
	$filStats =getFileStr($inF,0,70);
	my ($totRds,$Rejected1,$Rejected2,$Accepted1,$Accepted2,$Singl1,$Singl2,$AvgLen,$MaxLength,$AvgQual,$accErr) = 
		("0","0","0","0","0","0","0","0","0","0","0");
	if ($filStats eq ""){
	}elsif ( $filStats =~ m/Reads processed: ([0-9,]+); ([0-9,]+) \(pa/) ##  paired read mode..
			{ #only do this if newest format
		$totRds =  remComma($1) + remComma($2);
		if ($filStats =~ m/Rejected: ([0-9,]+); ([0-9,]+)\n/){
			$Rejected1 =  remComma($1); $Rejected2 =  remComma($2);
		}
		
		if ($filStats =~ m/Accepted \(High qual\): ([0-9,]+); ([0-9,]+)\s/){
			$Accepted1 =  remComma($1); $Accepted2 =  remComma($2);
		}
		#die "$Accepted1\n";
		#Singletons among these: 269,516; 6,686
		if ($filStats =~ m/Singletons among these: ([0-9,]+); ([0-9,]+)\n/){
			$Singl1 =  remComma($1); $Singl2 =  remComma($2);
		}
		
		if ($filStats =~ m/- [Ss]equence Length :\s*\d+.*\/([^\/]+)\/([^\/]+)\n/){
			$AvgLen = $1;  $MaxLength=$2;if ($MaxLengthHistBased > $MaxLength){$MaxLength = $MaxLengthHistBased;}
			#die "$AvgLen\n$filStats\n\n";
		}
#		$filStats =~ m/- Seq Length :\s*\d+.*\/(\d+.*)\/\d+.*\n/;
#		my $AvgLen = $1;
		if ($filStats =~ m/- Quality :\s*\d+.*\/(\d+.*)\/\d+.*\n/){
			$AvgQual = $1;
		}
		if ($filStats =~ m/- Accum. Error ([\d+\.]+)/){
			my $accErr = $1;
		}
		#$outStr .= "$totRds\t$Rejected1\t$Rejected2\t$Accepted1\t$Accepted2\t$Singl1\t$Singl2\t$AvgLen\t$MaxLength\t$AvgQual\t$accErr\t";
	} elsif ( $filStats =~ m/Reads processed: ([0-9,]+)/){#single end format
		$totRds =  remComma($1);
		if ($filStats =~ m/Rejected: ([0-9,]+)\n/){
			$Rejected1 =  remComma($1); $Rejected2 =  0;
		}
		if ($filStats =~ m/Accepted \(High qual\): ([0-9,]+)\s/ ) {
			$Accepted1 =  remComma($1); $Accepted2 =  0;
			$Singl1 =  $Accepted1; $Singl2 =  0;
		}
		if ($filStats =~ m/- [Ss]equence Length :\s*\d+.*\/([^\/]+)\/([^\/]+)\n/){
			$AvgLen = $1; $MaxLength=$2; if ($MaxLengthHistBased > $MaxLength){$MaxLength = $MaxLengthHistBased;}
			#die "$AvgLen\n$filStats\n\n";
		}
		if ($filStats =~ m/.*- Quality :\s*\d+.*\/(\d+.*)\/\d+.*\n/){
			$AvgQual = $1;
		}
		
		if ($filStats =~ m/.*- Accum\. Error ([\d\.]+)/){
			$accErr = $1;
		}
		#$outStr .= "$totRds\t$Rejected1\t$Rejected2\t$Accepted1\t$Accepted2\t$Singl1\t$Singl2\t$AvgLen\t$MaxLength\t$AvgQual\t$accErr\t";
	}
	if ($MaxLengthHistBased > $MaxLength){$MaxLength = $MaxLengthHistBased;}
	my @ret = ($totRds,$Rejected1,$Rejected2,$Accepted1,$Accepted2,$Singl1,$Singl2,$AvgLen,$MaxLength,$AvgQual,$accErr);
	return @ret;
	
}


sub bwtLogRd($$$){
	my ($splAr,$idx,$rhr) = @_;
	my @spl = @{$splAr};
	my %ret = %{$rhr};
	#DEFAULTS:
	$ret{uniqAlign}=-1;$ret{multAlign}=0;$ret{DisconcAlign}=0;
	if (@spl < 5){return (\%ret);}
	if ($spl[$idx+1] =~ m/\(100.00%\) were unpaired; of these:/){#single end read mapping!
		#die "SE\n";
		if ($spl[$idx+2] =~ m/(\d+) \(.+\) aligned 0 times/){ $ret{notAlign} = $1;} else { $ret{notAlign}=0; print "bwtOut wrg1.1\n";}
		 $ret{uniqAlign}=0;
		 $ret{multAlign}=0;
		$ret{DisconcAlign}=0;
		if ($spl[$idx+3] =~ m/(\d+) \(.+\) aligned exactly 1 time/ ){ $ret{SinglAlign} = $1;} else { $ret{SinglAlign}=0;print "bwtOut wrg1.5 $spl[14]\n";}
		if ($spl[$idx+4] =~ m/(\d+) \(.+\) aligned >1 times/ ){ $ret{SinglAlignMult} = $1;} else { $ret{SinglAlignMult}=0;print "bwtOut wrg1.6 $spl[15]\n";}
		if ( $spl[$idx+5] =~ m/(\d+\.\d+)\% overall alignment rate/ ){ $ret{AlignmRate} = $1;} else {$ret{AlignmRate} = -1; print "bwtOut wrg1.7 $spl[16]\n";}
	} else {
		#die "Pair\n";
		if ($spl[$idx+2] =~ m/(\d+) \(.+\) aligned concordantly 0 times/){ $ret{notAlign} = $1;} else { $ret{notAlign}=0; print "bwtOut wrg1\n";}
		if ($spl[$idx+3] =~ m/(\d+) \(.+\) aligned concordantly exactly 1 time/ ){ $ret{uniqAlign} = $1;} else { $ret{uniqAlign}=0;print "bwtOut wrg2\n";}
		if ($spl[$idx+4] =~ m/(\d+) \(.+\) aligned concordantly >1 times/ ){ $ret{multAlign} = $1;} else { $ret{multAlign}=0;print "bwtOut wrg3  $spl[6]\n";}
		if ($spl[$idx+7] =~ m/(\d+) \(.+\) aligned discordantly 1 time/ ){ $ret{DisconcAlign} = $1;} else { $ret{DisconcAlign}=0;print "bwtOut wrg4 $spl[9]\n";}
		if ($spl[$idx+12] =~ m/(\d+) \(.+\) aligned exactly 1 time/ ){ $ret{SinglAlign} = $1;} else { $ret{SinglAlign}=0;print "bwtOut wrg5 $spl[14]\n";}
		if (@spl>($idx+12) && $spl[$idx+13] =~ m/(\d+) \(.+\) aligned >1 times/ ){ $ret{SinglAlignMult} = $1;} else { $ret{SinglAlignMult}=0;print "bwtOut wrg6 $spl[15]\n";}
		if (@spl>($idx+13) && $spl[$idx+14] =~ m/(\d+\.\d+)\% overall alignment rate/ ){ $ret{AlignmRate} = $1;} else {$ret{AlignmRate} = -1; print "bwtOut wrg7 $spl[16]\n";}
	}
	return (\%ret);
}


sub getMapStats{
	my ($inP) = @_;
	my $inFi = "$inP/map.sh.etxt"; 
	$inFi = "$inP/bwtMap.sh.etxt" if (!-e $inFi && -e "$inP/bwtMap.sh.etxt"); #old MF file names..
	#my $outStrDesc = "";
	my @spl = ();
	my $alignStats = getFileStr($inFi,0);
	if (-s $inFi){
		@spl = split(/\n/,$alignStats);
	}
	#die "$alignStats\n";
	my $idx =-1;my $dobwtStat=0;
	if ($alignStats =~ m/strobealign/){ #m/\[M::worker_pipeline/){
		$dobwtStat=2;
	}elsif($alignStats =~ m/This is strobealign/){
		$dobwtStat=3;
	}elsif (@spl > 12 && $alignStats =~ m/reads; of these:/){
		$dobwtStat=1;
		while($spl[$idx] !~ m/\d+ reads; of these:/){
			$idx++; 
			if ($idx >= @spl){$idx=-1;last;}
		}
		if ($idx>=0 && $spl[$idx+0] =~ m/(\d+) reads; of these:/){
			$locStats{totReadPairs} = $1;
		} else {
			$dobwtStat=0;#die "wrong bwtOut: $inD/LOGandSUB/bwtMap.sh.etxt X $spl[2]";
		}
	}
	#die "$spl[$idx]\n$locStats{totReadPairs}\n";
	my $outStrDesc ="ReadsPaired\tAlignedReads\tOverallAlignment\tUniqueAlgned\tMultAlign\tDisconcAlign\tSingleUniqAlign\tSingleMultiAlign\t";
	my $outStr = "\t" x 8;
	
	my $incoming =0;my $retained =0;my $removed = 0;
	if ($alignStats =~ m/^Inentries: (\d+)/){
		my @matc = $alignStats =~ m/^Inentries: (\d+)/g;		foreach (@matc) {$incoming += int($_);}
		 @matc =$alignStats =~ m/^TotalRetained: (\d+)/g;	foreach (@matc) {$retained += int($_);}
		 @matc = $alignStats =~ m/^TotalRm: (\d+)/g;	foreach (@matc) {$removed += int($_);} 
		 $locStats{totReadPairs} = $incoming;
		 $locStats{uniqAlign} = $retained;
	} else {
		$locStats{totReadPairs} = -1;
		$locStats{uniqAlign} = -1;
	}
	
	if (!$dobwtStat){
		#$outStr .= "\t" x 7;
		#die "X!\n";
	} elsif ($dobwtStat == 1){#$alignStats =~ m/reads; of these:/){
		my ($rhr) = bwtLogRd(\@spl,$idx,\%locStats);
		%locStats = %{$rhr};
		$outStr = "$locStats{totReadPairs}\t$retained\t$locStats{AlignmRate}\t".$locStats{uniqAlign}/$locStats{totReadPairs}*100 ."\t".$locStats{multAlign}/$locStats{totReadPairs}*100 ."\t".$locStats{DisconcAlign}/$locStats{totReadPairs}*100 
		."\t".$locStats{SinglAlign}/$locStats{totReadPairs}*100 ."\t".$locStats{SinglAlignMult}/$locStats{totReadPairs}*100 ."\t";
	} elsif ($dobwtStat == 2){ #minimap2
		my @matches = ($alignStats =~ m/\[M::worker_pipeline::.*\] mapped (\d+) sequences/g);
		#print "@matches\n";
		my $sum=0; $sum += $_ foreach (@matches);
		if ($sum>0){
			my $frac = (1 - ($removed/$sum)) * 100;
			$outStr = "$sum\t$retained\t".  $frac ."\t\t\t\t\t\t";
		}
		#die "minimap!!$sum\n";
	} elsif ($dobwtStat == 3){ #strobealign
		if ($incoming > 0){
			my $frac2 = ($incoming-$removed)/$incoming;
			#$locStats{totReadPairs} = -1 if (!exists($locStats{totReadPairs}));
			my $frac = 0; $frac = $incoming/$locStats{totReadPairs} if( $locStats{totReadPairs}>0);
			$locStats{AlignmRate}=$frac;
			$outStr = "$locStats{totReadPairs}\t$retained\t".  $frac ."\t"  .$locStats{uniqAlign}/$locStats{totReadPairs}*100 . "\t\t\t\t\t";
		}
	}
	return ($outStr,$outStrDesc);
}
sub optiDups{
	my ($inP) = @_;
	my $inFi = "$inP/map2.sh.etxt"; 
	$inFi = "$inP/bwtMap2.sh.etxt" if (!-e $inFi); #old MF file names..
	#my $outStrDesc = "";
	#novocraft sort
	my $doDup = 0;my $alignStats2 = getFileStr($inFi,0);
	if ($alignStats2 ne "" ){if (defined ($alignStats2)){$doDup=1;}}
	if ($doDup){
		
		if ($alignStats2 =~ m/samtools markdup/){#samtools dedup
			my $pairDup =0; my $singlDup=0; my $optiPairDup=0; my $optiSinglDup=0;
			my $Npairs=0; my $Nsingle=0;
			if ($alignStats2 =~ m/PAIRED: (\d+)/){$Npairs=$1;}
			if ($alignStats2 =~ m/SINGLE: (\d+)/){$Nsingle=$1;}
			if ($alignStats2 =~ m/WRITTEN: (\d+)/){$locStats{duplPass}=$1;}
			if ($alignStats2 =~ m/DUPLICATE PAIR: (\d+)/){$pairDup=$1;}
			if ($alignStats2 =~ m/DUPLICATE SINGLE: (\d+)/){$singlDup=$1;}
			if ($alignStats2 =~ m/DUPLICATE PAIR OPTICAL: (\d+)/){$optiPairDup=$1;}
			if ($alignStats2 =~ m/DUPLICATE SINGLE OPTICAL: (\d+)/){$optiSinglDup=$1;}
			if ($alignStats2 =~ m/ESTIMATED_LIBRARY_SIZE: (\d+)/){$locStats{EstLibSize}=$1;}
			$locStats{duplPCR} = ($pairDup+$singlDup);#/($Npairs+$Nsingle);
			$locStats{duplOptic} = ($optiPairDup+$optiSinglDup);#/($Npairs+$Nsingle);
			
		} else {
			#print "$alignStats2 YUYS\n";
			if ($alignStats2 =~ m/Proper Pairs\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/){$locStats{EstLibSize} = $4; $locStats{duplOptic} = $3; $locStats{duplPCR} = $2; $locStats{duplPass} = $1;} else {$doDup=0;}#die "Can't find dupl stats: $alignStats2\n";}
			if ($alignStats2 =~ m/Improper Pairs\s+(\d+)\s+(\d+)\s+(\d+)/){$locStats{duplOptic} += $3; $locStats{duplPCR} += $2; $locStats{duplPass} += $1;} else {$doDup=0;}#{die "Can't find impr. dupl stats: $alignStats2\n";}
		}
	}
	my $outStr = "";
	if ($doDup){
		$outStr = "$locStats{duplOptic}\t$locStats{duplPCR}\t$locStats{duplPass}\t$locStats{EstLibSize}\t";
	} else {
		$outStr = "\t" x 4;
	}
	my $outStrDesc = "OpticalDuplicates\tPCRduplicates\tPassedMD\tEstLibSize\t";
	return ($outStr,$outStrDesc);
}

sub getContamination{
	my ($inFi,$inFi2) = @_;
	my $outStr = "?\t?\t?\t";my $outStrDesc = "";$outStrDesc .= "FilteredContaRdsPerc\tFilteredContaRds\tFilteredNonContaRds\t";
	my $filStats = getFileStr($inFi,0);#`cat $inD/LOGandSUB/KrakHS.sh.etxt`; chomp $filStats;
	#if ($filStats eq "" ){$outStr .= "?\t?\t?\t";	return ($outStr,$outStrDesc);}
	
	my @matches = ($filStats =~ m/\d+ sequences classified \((\d+\.?\d*)%\)/g);
	my @hits = ($filStats =~ m/(\d+) sequences classified \(\d+\.?\d*%\)/g);
	my @nonhits = ($filStats =~ m/(\d+) sequences unclassified \(\d+\.?\d*%\)/g);
	
	if (@nonhits == 0){#check if this was done via hostile..
		$filStats = getFileStr($inFi2,0); #paired reads should be counted as two.. but are counted as one in hostile..
		@hits = ($filStats =~ m/"reads_removed": (\d*),/g); #"reads_removed": 202,
		@nonhits = ($filStats =~ m/"reads_out": (\d*),/g); 
		@matches = ($filStats =~ m/"reads_removed_proportion": (\d*),/g); 
	}
	
	
	if (@nonhits > 0){
		my $totHits=0; my $totNH=0;
		for (my $i=0;$i<@hits;$i++){
			$totHits += $hits[$i];$totNH += $nonhits[$i];
		}
		
		$outStr = join(";",@matches) . "\t$totHits\t$totNH\t";
	} 
	
	
	
	#die "$outStr\n$inFi\n";
	return ($outStr,$outStrDesc);
}

sub getBinnerStats{
	my ($tmpassD,$SmplN) = @_;
	
## binning stats SemiBin
	my $SBbinCM2 = "$tmpassD/Binning/SB/$SmplN.cm2";
	my $addStr = ""; my $addDescr="";
	if (-e $SBbinCM2){
		my $HQbinCnt = 0; my $MQbinCnt = 0; my $totBins=0;
		open I,"<$SBbinCM2" or die $!;
		while (<I>){my @spl = split/\t/; next if ($spl[1] eq "Completeness");
			if ($spl[1] >= 90 && $spl[2] <= 5){$HQbinCnt ++ ;
			} elsif ($spl[1] >= 80 && $spl[2] <= 5){$MQbinCnt ++ ;}
			$totBins ++;
		}
		close I;
		$addStr = "$HQbinCnt\t$MQbinCnt\t";
		#die "$HQbinCnt $MQbinCnt $SBbinCM2\n";
	} else {
		$addStr .= "\t" x 2
	}
	$addDescr .= "HQ_bins_SB\tMQ_bins_SB\t";
	


## binning stats MetaBat2
	$SBbinCM2 = "$tmpassD/Binning/MB2/$SmplN.cm2";
	if (-e $SBbinCM2){
		my $HQbinCnt = 0; my $MQbinCnt = 0;my $totBins=0;
		open I,"<$SBbinCM2" or die $!;
		while (<I>){my @spl = split/\t/; next if ($spl[1] eq "Completeness");
			if ($spl[1] >= 90 && $spl[2] <= 5){$HQbinCnt ++ ;
			} elsif ($spl[1] >= 80 && $spl[2] <= 5){$MQbinCnt ++ ;}
			$totBins ++;
		}
		close I;
		$addStr = "$HQbinCnt\t$MQbinCnt\t";
	} else {
		$addStr .= "\t" x 2
	}
	$addDescr .= "HQ_bins_MB2\tMQ_bins_MB2\t";

## binning stats MetaDecoder
	$SBbinCM2 = "$tmpassD/Binning/MD/$SmplN.cm2";
	if (-e $SBbinCM2){
		my $HQbinCnt = 0; my $MQbinCnt = 0;my $totBins=0;
		open I,"<$SBbinCM2" or die $!;
		while (<I>){my @spl = split/\t/; next if ($spl[1] eq "Completeness");
			if ($spl[1] >= 90 && $spl[2] <= 5){$HQbinCnt ++ ;
			} elsif ($spl[1] >= 80 && $spl[2] <= 5){$MQbinCnt ++ ;}
			$totBins ++;
		}
		close I;
		$addStr = "$HQbinCnt\t$MQbinCnt\t";
	} else {
		$addStr .= "\t" x 2
	}
	$addDescr .= "HQ_bins_MD\tMQ_bins_MD\t";
	return ($addStr,$addDescr);
}

sub getSNPStats{
	my ($inFi) = @_;
	#my $outStrDesc = "";  my $outStr = "";
	my $outStr = "\t" x 4;
	my $outStrDesc = "SNPbpResolved	SNPfastaEntries	SNPconflicts\tSNPconlResolv2ndL\tNumSNPs\t";
	my $geneStats = getFileStr("$inFi",0);
	#Total bp written: 34950480 (0 not resolved) on 26727 entries

	if ($geneStats =~ m/Total bp written: (\d+) \(\d+ not resolved\) on (\d+) entries/){
		my $bps = $1; my $entrs=$2;
		$geneStats =~ m/Conflicting calls: (\d+) Resolved with second line: (\d+)/;
		my $confl=$1; my $resol=$2;
		$geneStats =~ m/Total SNPs detected: (\d+)/; my $snpNum=$1;
		$outStr = "$bps\t$entrs\t$confl\t$resol\t$snpNum\t";
	}
	return ($outStr, $outStrDesc);
}

sub getGeneStats{
	my ($inFi) = @_;
	#my $outStrDesc = "";  my $outStr = "";
	my $outStr = "\t" x 9;
	my $outStrDesc = "GeneNumber	AvgGeneLength	AvgComplGeneLength	BpGenes	BpNotGenes	Gcomplete	G5pComplete	G3pComplete	Gincomplete	";
	my $geneStats = getFileStr("$inFi",0);
	my @spl1 = split("\n", $geneStats);

	#if (!-e "$inFi" || -s "$inFi" == 0){	
	#	return ($outStr, $outStrDesc);
	#}
	#open I,"<$inFi" or die $!; 
	#GeneNumber	AvgGeneLength	AvgComplGeneLength	BpGenes	BpNotGenesGcomplete	G5pComplete	G3pComplete	Gincomplete
	#my $tmp ="";$tmp = <I>;#if (defined($tmp)){chomp($tmp);$outStrDesc .= $tmp."\t";}$tmp = <I>;
	#while (my $tmp = <I>){
	foreach my $tmp (@spl1){
		chomp($tmp);#$outStr = $tmp."\t"; $outStr5 .= $tmp."\t" if ($do500Stat); 
		next if ($tmp =~ m/^GeneNumber/);
		my @spl = split /\t/,$tmp; 
		#print @spl . " @spl\n";
		if (@spl >=9){
			#bit convoluted way of doing this, but better to be conservative with number of tabs used..
			$outStr = "$spl[0]\t$spl[1]\t$spl[2]\t$spl[3]\t$spl[4]\t$spl[5]\t$spl[6]\t$spl[7]\t$spl[8]\t";
		}
	}
	#print $outStr;
	#close I;
	return ($outStr, $outStrDesc);
}
sub getASsemblyStats{
	my ($tmpassD,$assemblStatsFile,$doCirc) = @_;
	my $outStr = "";my $outStrDesc = "";
	my $assStats = getFileStr("$tmpassD/$assemblStatsFile",0);
	if ($assStats ne ""){
		if ($assStats =~ m/N50 contig length\s+(\d+)/){ $locStats{CtgN50} = $1;} else { $locStats{CtgN50} = 0; print "$assemblStatsFile\ncntgn50 wrg1\n";}
		if ($assStats =~ m/Number of scaffolds\s+(\d+)/){ $locStats{NScaff} = $1;} else { die "$assemblStatsFile\nscfcnt wrg1\n";}
		if ($assStats =~ m/Total size of scaffolds\s+(\d+)/){ $locStats{ScaffSize} = $1;} else { die "scfcnt wrg2\n";}
		if ($assStats =~ m/Longest scaffold\s+(\d+)/){ $locStats{ScaffMaxSize} = $1;} else { die "scfcnt wrg3\n";}
		if ($assStats =~ m/N50 scaffold length\s+(\d+)/){ $locStats{ScaffN50} = $1;} else {$locStats{ScaffN50}=-1;}# die "$tmpassD\nscfcnt wrg4\n";}
		if ($assStats =~ m/Number of scaffolds > 1K nt\s+(\d+)/){ $locStats{NScaffG1k} = $1;} else { die "scfcnt wrg5\n";}
		if ($assStats =~ m/Number of scaffolds > 10K nt\s+(\d+)/){ $locStats{NScaffG10k} = $1;} else { die "scfcnt wrg6\n";}
		if ($assStats =~ m/Number of scaffolds > 100K nt\s+(\d+)/){ $locStats{NScaffG100k} = $1;} else { die "scfcnt wrg7\n";}
		if ($assStats =~ m/Number of scaffolds > 1M nt\s+(\d+)/){ $locStats{NScaffG1M} = $1;} else { die "scfcnt wrg8\n";}
		$outStr .= "$locStats{CtgN50}\t$locStats{NScaff}\t$locStats{NScaffG1k}\t$locStats{NScaffG10k}\t$locStats{NScaffG100k}\t$locStats{NScaffG1M}\t$locStats{ScaffN50}\t$locStats{ScaffMaxSize}\t$locStats{ScaffSize}\t";
	} else {
		$outStr .= "\t" x 9;
	}
	$outStrDesc .="ContigN50\tNScaff400\tNScaffG1k\tNScaffG10k\tNScaffG100k\tNScaffG1M\tScaffN50\tScaffMaxSize\tScaffSize\t";
	my $circCtgs=0;my $circg1M=0; 
	if (-e "$tmpassD/scaffolds.fasta.circ"){
		$assStats = getFileStr("$tmpassD/scaffolds.fasta.circ",0);
		$circCtgs = $assStats =~ tr/>//;
		my @matchs = ($assStats =~ m/.*_L=(\d+)=/g);
		foreach(@matchs){$circg1M++ if ($_ > 1000000);}
		#die "$circCtgs  XX $g1M YY @matchs\n".@matchs;
		#$circCtgs = `wc -l $tmpassD/scaffolds.fasta.circ | cut -f1 -d' '`; chomp $circCtgs;
	}
	$outStr .= "$circCtgs\t$circg1M\t";
	$outStrDesc .= "CircCtgs\tCircCtgG1M\t";
	return ($outStr,$outStrDesc);
}

sub smplStats(){
	my ($inD,$assDir,$SmplN) = @_;
	my $do500Stat = 0;
	#print "STATS   \n";
	#my %ret; 
	my $outStr = ""; my $outStrDesc = "";my $outStr5 = ""; my $outStrDesc5 = "";
	

#first report file size and if this is paired or single..
	if (exists($inputFileSizeMB{$curSmpl})){
		$outStr .= int($inputFileSizeMB{$curSmpl}/1024)."G\t";
	} else {$outStr .= "-1\t";}
	$outStrDesc .= "RawInputSize\t";
	$outStr .= $locStats{hasPaired}. "\t" . $locStats{hasSingle}. "\t";
	$outStrDesc .= "InputIsPaired\tInputIsSingle\t";


	my ($t1,$t2) = getContamination("$inD/LOGandSUB/KrakHS.sh.etxt","$inD/LOGandSUB/KrakHS.sh.otxt");
	$locStats{contamination} = $t1;
	$outStr .= $t1;	$outStrDesc .= $t2;
	if ($do500Stat){		$outStr5 .= $t1;	$outStrDesc5 .= $t2;	}

		
	my @sdmStat ;#(0,0,0,0,0,0,0,0,0,0,0);
	if (-s "$inD/LOGandSUB/sdm/filter.log" || -s "$inD/LOGandSUB/sdm/filterS.log"){
		@sdmStat = sdmStats("$inD/LOGandSUB/sdm/filter.log",$inD);# if (-s "$inD/LOGandSUB/sdm/filter.log");
		#@sdmStat = @{$ar};
		#my @sdsm2 =(0,0,0,0,0,0,0,0,0,0,0);
		my @sdsm2 = sdmStats("$inD/LOGandSUB/sdm/filterS.log",$inD);# if (-s "$inD/LOGandSUB/sdm/filterS.log");
		#print "$sdsm2[0]\n$sdsm2[3]\n";
		$sdmStat[0] += $sdsm2[0];$sdmStat[1] += $sdsm2[1];$sdmStat[2] += $sdsm2[2];$sdmStat[3] += $sdsm2[3];
		$sdmStat[4] += $sdsm2[4];$sdmStat[5] += $sdsm2[5];$sdmStat[6] += $sdsm2[6];
		$sdmStat[7] = $sdsm2[7] if ($sdmStat[7] == 0);
		$sdmStat[8] = $sdsm2[8] if ($sdmStat[8] == 0);
		$sdmStat[9] = $sdsm2[9] if ($sdmStat[9] == 0);
		$sdmStat[10] = $sdsm2[10] if ($sdmStat[10] == 0);
		
	} else {#(-s "$inD/LOGandSUB/sdmReadCleaner.sh.etxt") {
		@sdmStat = sdmStats("$inD/LOGandSUB/sdmReadCleaner.sh.etxt",$inD);
	}
	$locStats{totRds} = $sdmStat[0];$locStats{Rejected1} = $sdmStat[1];$locStats{Rejected2} = $sdmStat[2];
	$locStats{Accepted1} = $sdmStat[3];$locStats{Accepted2} = $sdmStat[4];$locStats{Singl1} = $sdmStat[5];
	$locStats{Singl2} = $sdmStat[6];
#	if ($sdmStat[0] == 0){
#		$outStr .= "\t" x 11;
#	} else {
	#print "@sdmStat"." xx".@sdmStat."\n";
	$outStr .= join("\t",@sdmStat)."\t";
#	}
## 11 in total..
	$outStrDesc .= "totRds\tRejected1\tRejected2\tAccepted1\tAccepted2\tSingl1\tSingl2\tAvgSeqLen\tMaxSeqLength\tAvgSeqQual\taccErr\t";
	
	
	#check for flash merged reads
	my $filStats = "";
	$filStats = getFileStr("$inD/LOGandSUB/flashMrg.sh.otxt",0);
	if ($filStats ne ""){
		
		if ($filStats =~ m/\[FLASH\]     Combined pairs:   (\d+)/){$outStr.="$1\t";} else {$outStr.="?\t";}
		if ($filStats =~ m/\[FLASH\]     Uncombined pairs: (\d+)/){$outStr.="$1\t";} else {$outStr.="?\t";}
	} else {
		$outStr .= "\t\t";
	}
	$outStrDesc .= "Merged\tNotMerged\t";
#geno size estimate
	$filStats = getFileStr("$inD/MicroCens/MC.0.result",0);
	if ($filStats ne ""){
		if ($filStats =~ m/average_genome_size:	([\d\.]+)/){$outStr.="$1\t";} else {$outStr.="?\t";}
		if ($filStats =~ m/genome_equivalents:	([\d\.]+)/){$outStr.="$1\t";} else {$outStr.="?\t";}
	} else {
		$outStr .= "\t\t";
	}
	$outStrDesc .= "AvgGenomeSizeEst\tTotalGenomesEst\t";
	
	
	#check if corrected dir still exists..
	system "rm -rf $inD/assemblies/metag/corrected" if (-d "$inD/assemblies/metag/corrected");
	#assembly stats
	my $tmpassD = "";
	if (-e "$inD/assemblies/metag/assembly.txt"){
		$tmpassD = getFileStr("$inD/assemblies/metag/assembly.txt",0); chomp $tmpassD;
	} elsif (-e "$inD/assemblies/metag/AssemblyStats.txt"){
		$tmpassD = "$inD/assemblies/metag/";
	} elsif (-e $assDir."metag/AssemblyStats.txt"){
		$tmpassD = "$assDir/metag/";
	}
#print $tmpassD."\n";
	($t1,$t2) = getASsemblyStats($tmpassD,"AssemblyStats.txt",1);
	$outStr .= $t1;	$outStrDesc .= $t2;
	if ($do500Stat){
		my ($t1,$t2) = getASsemblyStats($tmpassD."AssemblyStats.500.txt",0);
		$outStr5 .= $t1;	$outStrDesc5 .= $t2;
	}
	

	($t1,$t2) = getMapStats("$inD/LOGandSUB/");
	$outStr .= $t1;	$outStrDesc .= $t2;
	if ($do500Stat){$outStr5 .= $t1;	$outStrDesc5 .= $t2;}
	($t1,$t2) = optiDups("$inD/LOGandSUB/");
	$outStr .= $t1;	$outStrDesc .= $t2;
	if ($do500Stat){$outStr5 .= $t1;	$outStrDesc5 .= $t2;}

	#return ($outStrDesc,$outStr,$outStrDesc5, $outStr5); #DEBUG


	# stats on gene number etc
	($t1,$t2) = getGeneStats("$inD/$dir_ContigStats/GeneStats.txt");
	$outStr .= $t1;	$outStrDesc .= $t2;
	if ($do500Stat){$outStr5 .= $t1;	$outStrDesc5 .= $t2;}

	#MB2, SemiBin etc hq & mq MAGs
	($t1,$t2) = getBinnerStats($tmpassD,$SmplN);
	$outStr .= $t1;	$outStrDesc .= $t2;
	if ($do500Stat){$outStr5 .= $t1;	$outStrDesc5 .= $t2;}

	($t1,$t2) = getSNPStats("$inD/LOGandSUB/SNP/ConsAssem.oSNPc.sh.etxt");
	$outStr .= $t1;	$outStrDesc .= $t2;
	if ($do500Stat){$outStr5 .= $t1;	$outStrDesc5 .= $t2;}
	
#clean up strings..	
	chomp $outStrDesc; chomp $outStr;
	if ($do500Stat){
		chomp $outStrDesc5; chomp $outStr5;
	} else {
		$outStrDesc5 = "";  $outStr5 = "";
	}
#	my @gest = split(/\t/,$gsta);
	#my $nGenes = `wc -l $inD/assemblies/metag/genePred/genes.gff | cut -f1 -d ' '`;	chomp $nGenes;
	#die ($outStr."\n$nGenes\n");
#	$outStr .= "$nGenes\t";	$outStrDesc .= "NumGenes\t";

	#print "... END\n";
	
	
	return ($outStrDesc,$outStr,$outStrDesc5, $outStr5);
	#die $outStr."\n";
}

sub spadesHosts{
	#figure out if only certain node subset has enough HDD space
	if (0 && $HDDspace{spades} > 100){
		my $locHosts = `bhosts | grep ok | cut -d" " -f 1 | grep compute | tr "\\n" ","`;
		my $tmpStr = `pdsh -w $locHosts -u 4 "df -l" `;
		my $srchTerm = '/$';
		if (`hostname` =~ m/submaster$/){
			$srchTerm = '/tmp';
		}
		#print "psdh done\n";
		#die $tmpStr;
		#54325072 = 52G
		foreach my $l (split(/\n/,$tmpStr)){
			next if ($l !~ m/compute\S+:/);
			next unless ($l =~ m/$srchTerm/);
			#print $l."\n";
			my @hosts = split(/:/,$l);
			my @spl = split(/\s+/,$hosts[1]);
			#die $spl[1]." XX ".$spl[2]." XX ".$spl[3]." XX ".$spl[4]." \n ";
			#print $spl[4] / 1024 / 1024 . "\n";
			if (($spl[4] / 1024 / 1024) > $HDDspace{spades}){
				push(@{$QSBoptHR->{Spades_Hosts}},$hosts[0]);
			}
			if (($spl[4] / 1024 / 1024) > 40){
				push(@{$QSBoptHR->{General_Hosts}},$hosts[0]);
			}
		}
		print "Found ".scalar @{$QSBoptHR->{Spades_Hosts}}." host machines with > $HDDspace{spades} G space\n";
		print "Found ".scalar @{$QSBoptHR->{General_Hosts}}." host machines with > 40 G space\n";
		if (scalar @{$QSBoptHR->{Spades_Hosts}} ==0 ){die "Not enough hosts for spades temp space found\n";}
		sleep (2);
	}
}


sub readG2M($){
	my ($inf) = @_;
	open K,"<",$inf;
	my %COGs2motus;
	my %motus;
	while (my $li = <K>){
		chomp($li);
		my @spl = split("\t",$li);
		$motus{$spl[0]} = $spl[3];
		$COGs2motus{$spl[0]} = $spl[2];
	}
	close K;
	print "Read gene to map file\n";
}

sub scndMap2Genos{
	# ----------------- 2nd mapping (map to ref genomes supplied by user) ---------------------
	my ($SmplName,$cleanSeqSetHR,$cMapGrp,$cAssGrp,$curOutDir,$nodeSpTmpD,$smplTmpDir,
			$sampleDepsAR,$samplReadLength,$calc2ndMapSNP,$boolScndCoverageOK) = @_;
	my @mapOutXS;my @bamBaseNameS;
	for (my $i=0;$i<@bwt2outD; $i++){
		$bwt2outD[$i] =~ m/\/([^\/]+)\/*$/;
		my $smplXDB = $1;
		push @mapOutXS, $smplTmpDir."xtraMaps/$smplXDB/";
		my $fname = $bwt2ndMapNmds[$i]."_".$SmplName."-0";
		#if (length($fname ) > 20){$bwt2ndMapNmds[$i]."_".$SmplName."-0";
		push @bamBaseNameS, $fname;
	}
#		system "rm -rf ".join(" ",@mapOutXS) if ($redo2ndMapping);

	$make2ndMapDecoy{Lib} = $curOutDir if ($MFopt{DoMapModeDecoy});
	my $cramthebam=0;
	#map to all refs at once		
	my %dirset = 	(nodeTmp=>$nodeSpTmpD,outDir => join(",",@bwt2outD),unalDir=>"",
					sbj => join(",",@DBbtRefX),assGrp => $cAssGrp,
					smplName => join(",",@bamBaseNameS),
					glbTmp => $smplTmpDir."/xtraMaps/", is2ndMap => 1, 
					qsubDir => "$logDir/map2nd/",mapSupport => 0,
					glbMapDir => join(",",@mapOutXS),
					readTec => ${$cleanSeqSetHR}{readTec}, #$map{$curSmpl}{SeqTech}
					submit => 1, submNow => 1, cramAlig => $cramthebam,
					sortCores => $MFopt{bamSortCores}, mapCores => $MFopt{MapperCores});
	#die "XX @DBbtRefX\n";
	# ----------------- 2nd mapping (map to ref genomes supplied by user) ---------------------
	my ($map2CtgsX,$delaySubmCmdX,$mapOptHr) = mapReadsToRef(\%dirset,\%AsGrps,$AsGrps{$cMapGrp}{SeqUnZDeps}.";".$bwt2ndMapDep );#$localAssembly);
	#------------  and calc coverage for each separate
	#different strategy, a bit hacky: deactivate submissions and collect for one call
	$dirset{submit} = 0; my $bigSort = "";my $bigCov = "";
	for (my $i=0;$i<@bwt2outD; $i++){
		 my $mapOutX = $mapOutXS[$i]; 
		$dirset{outDir} = $bwt2outD[$i];$dirset{glbMapDir} =$mapOutXS[$i];
		$dirset{smplName} = $bamBaseNameS[$i];$dirset{sbj} =$DBbtRefX[$i];
		my ($map2CtgsY,$delaySubmCmdY,$mapStat)  = bamDepth(\%dirset,$map2CtgsX,$mapOptHr);
		$bigSort .= "\n\n#----------$i -------------\n$delaySubmCmdY\n" if ($delaySubmCmdY ne "");
		
		#die "$map2CtgsY\n$delaySubmCmdY\n$mapStat\n";
		#call abundance on newly made genes
		if ($MFopt{mapModeCovDo} && !$boolScndCoverageOK){
			if ($mapStat ==2|| $mapStat==0){ #files still on scratch
				($map2CtgsY,$delaySubmCmdY) = calcCoverage("$mapOutX/"."$bamBaseNameS[$i]-smd.bam.coverage.gz",$DBbtRefGFF[$i] ,
					$samplReadLength,$SmplName.".$i",$map2CtgsY.";$bwt2ndMapDep", \%dirset) ;
				push(@{$AsGrps{$cMapGrp}{MapCopiesNoDel}},$mapOutX."/*",$bwt2outD[$i]);
			} else { #files are already copied
				#die "$bwt2outD[$i]/"."$bamBaseNameS[$i]-smd.bam.coverage.gz\nxx\n";
				my $empty;
				($empty,$delaySubmCmdY) = calcCoverage("$bwt2outD[$i]/"."$bamBaseNameS[$i]-smd.bam.coverage.gz",$DBbtRefGFF[$i] ,
					$samplReadLength,$SmplName.".$i",$map2CtgsY.";$bwt2ndMapDep", \%dirset) ;
			}
			#die "ASD\n";
		} elsif ($mapStat == 2 || $mapStat==0){#check if files need to be copied..
			push(@{$AsGrps{$cMapGrp}{MapCopiesNoDel}},$mapOutX."/*",$bwt2outD[$i]);
		}
		#$AsGrps{$cMapGrp}{MapDeps} .= $map2CtgsY.";";  -> not needed if all submissions happen later
		$bigCov .= "\n\n#---------- $i -------------\n$delaySubmCmdY\n" if ($delaySubmCmdY ne "");
		
		
		
		if ($calc2ndMapSNP){ #2nd map SNP calling (consensus)
			my %SNPinfo = (
				assembly => "$bwt2outD[$i]/$bwt2ndMapNmds[$i].fa",#$DBbtRefX[$i],
				MAR => ["$bwt2outD[$i]/$bamBaseNameS[$i]-smd.bam"],
				SNPcaller => $MFopt{SNPcallerFlag},bamcram=>"bam",
				#doesn't work, if contig name and length is not given..
				#depthF => "$bwt2outD[$i]/$bamBaseNameS[$i]-smd.bam.coverage.gz.percontig",
				ofas => "$bwt2outD[$i]/$bamBaseNameS[$i].SNPc.$MFopt{SNPcallerFlag}.fna", #only output needed for this.. unless I later want to add also a gene calling.. (not needed for TEC2 reb)
				#vcfFile => "$bwt2outD[$i]/$bamBaseNameS[$i].$MFopt{SNPcallerFlag}.vcf",
				firstInSample => ($i == 0 ? 1 : 0), 
				SeqTech => $map{$curSmpl}{SeqTech}, SeqTechSuppl => $map{$curSmpl}{seqTechX},
				nodeTmpD => $dirset{nodeTmp}, #$nodeSpTmpD,
				scratch => $dirset{glbTmp}, #"$smplTmpDir/SNP/",
				qsubDir => $dirset{qsubDir}, jdeps => $map2CtgsY.";$bwt2ndMapDep",
				cmdFileTag => $bwt2ndMapNmds[$i], minDepth => $MFopt{consSNPminDepth},
				smpl => $bamBaseNameS[$i], maxCores => $MFopt{maxSNPcores}, memReq => $MFopt{memSNPcall},
				bpSplit => 4e5,	runLocal => 1, split_jobs => $MFopt{SNPconsJobsPsmpl}, overwrite => $MFopt{redoSNPcons},
				minCallQual => $MFopt{SNPminCallQual},
				);
				

			my $consSNPdep = createConsSNP(\%SNPinfo);
			add2SampleDeps($sampleDepsAR, [$consSNPdep]);
			#push(@sampleDeps, $consSNPdep) if (defined $consSNPdep && $consSNPdep ne "");
			#die if ($i==2);
		}
	
	}
	
	#actual job submission of concatenated sort & cov jobs..
	my ($sortJD,$tmpCmd) = qsubSystem($dirset{qsubDir}."SRTB$SmplName.sh",$bigSort,$MFopt{bamSortCores},int(20/$MFopt{bamSortCores})."G",$SmplName."SRT2nd",$map2CtgsX,"",1,[],$QSBoptHR);
	($sortJD,$tmpCmd) = qsubSystem($dirset{qsubDir}."COV$SmplName.sh",$bigCov,1,int(20)."G",$SmplName."COV2nd",$sortJD,"",1,[],$QSBoptHR);

	$AsGrps{$cMapGrp}{MapDeps} .= $sortJD.";";
	#only used for now for the mapping to spec ref
	my $cln1 = clean_tmp([],[], $AsGrps{$cMapGrp}{MapCopiesNoDel},$AsGrps{$cMapGrp}{MapDeps},"",
		"_mcl"."_$JNUM");
		#die"mcl";
	#print "XX $AsGrps{$cMapGrp}{MapDeps}\n";
	$AsGrps{$cMapGrp}{MapCopiesNoDel} = [];#$AsGrps{$cMapGrp}{MapDeps}="";
	$AsGrps{$cAssGrp}{scndMapping} .= $cln1.";";#tmp dir shouldn't be deleted before this is done
	#@cleans = {}; #just deactivate clean up for sec mapping..
}

sub ReadsFromMapping{
	my ($tmpFile,$linkFile) = @_;
	my $newOfile = ""; my $GID = "";
	open TO,">",$linkFile;
	open I,"<",$tmpFile;
	while (my $line = <I>){
		chomp($line);
		my @spl = split(/\t/,$line);
		my $readID = $spl[0];
		my $geneMtch = $spl[2];
		#get read and genomeMatch
		#unless (exists($motus{$geneMtch})){print "could not find gene ID $geneMtch\n"; next;}
		#my $newHd = "@".$COGs2motus{$geneMtch}."___".$motus{$geneMtch};
		$newOfile = $GID.".rdM";
		#die ($newOfile);
		
		$readID=~m/(.*)#/;
		print TO $1."\t"."\t".$newOfile."\n";#$newHd."\n".$spl[9]."\n+\n".$spl[10]."\n";
		#die;
	}
	close TO;
	close I;
}

sub RayAssembly(){
 "mpiexec -n 1 /g/bork5/hildebra/bin/Ray-2.3.1/ray-build/Ray -o test -p test/test_1.fastq test/test_2.fastq -k 31"
 }
 
 
sub buildAssemblyMapIdx{
	my ($finAssLoc,$cAssGrp, $mainRds, $suppRds, $smpl) = @_;
	print "Building mapper index for assembly $finAssLoc\n ";
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
	if ($suppRds){ #if support reads, build a second DB..
		my ($par1,$par2,$parS,$liar,$rear) = getRawSeqsAssmGrp(\%AsGrps,$cAssGrp,$suppRds,$smpl);
		my $MapperProgLoc = decideMapper($MFopt{MapperProg},${$liar}[0]);
		my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx($finAssLoc,$MFopt{MapperCores},$MFopt{largeMapperDB},$MapperProgLoc);#$nCores);
		my ($jname,$tmpCmd) = qsubSystem($logDir."mapperIdxSupp.sh",$cmdDB,(int($MFopt{MapperCores})),(int($MFopt{bwtIdxAssMem}/$MFopt{MapperCores})+1)."G","DBidx$JNUM","","",1,[],$QSBoptHR) ;
		$AsGrps{$cAssGrp}{AssemblJobName} .= ";$jname";
	}
	if ($mainRds){
		my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx($finAssLoc,$MFopt{MapperCores},$MFopt{largeMapperDB},$MFopt{MapperProg});#$nCores);
		my ($jname,$tmpCmd) = qsubSystem($logDir."mapperIdx.sh",$cmdDB,(int($MFopt{MapperCores})),1+(int($MFopt{bwtIdxAssMem}/$MFopt{MapperCores}))."G","bwtIdx$JNUM","","",1,[],$QSBoptHR) ;
		$AsGrps{$cAssGrp}{AssemblJobName} .= ";$jname";
	}
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	
}

 
 sub createPsAssLongReads(){
	my ($cleanSeqSetHR,$jdep, $pseudoAssFile, $Fdir, $smplName) = @_;
	my $arp1 = ${$cleanSeqSetHR}{arp1}; my $arp2 = ${$cleanSeqSetHR}{arp2}; my $singAr = ${$cleanSeqSetHR}{singAr}; 
	
	my @allRds = (@{$arp1},@{$arp2},@{$singAr});
	my $psDir = $pseudoAssFile; $psDir =~ s/[^\/]+$//;
	my $psFinal = $pseudoAssFile; $psFinal =~ s/^.+\//$Fdir\//;
	#die $psFinal;
	my $pseudoAssFileFlag = $pseudoAssFile.".sto";
	
	my $psFile = $pseudoAssFile;#"$finalCommAssDir/longReads.fasta.filt";
	if (-d $psFinal){system "rm -fr $psFinal";}
	if (!-e $pseudoAssFile && -e $psFinal && -e $psFinal.".sto"){
		$pseudoAssFileFlag = $psFinal.".sto"; 
		$psFile = $psFinal;
	}

	#die "@allRds\n";
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	my $cmd = "";
	$cmd .= "mkdir -p $psDir $Fdir\n";
	$cmd .= "$sizFiltScr ".join(",",@allRds)." $MFopt{scaffoldMinSize} -1 $pseudoAssFile\n";
	$cmd .= "$renameCtgScr $psFile $smplName\n";
	$cmd .= "touch $pseudoAssFileFlag\n";
	#die $cmd;
	my ($jname,$tmpCmd) = ("","");
	if (!-e $pseudoAssFileFlag || !-e $psFinal.".sto"){
		$jname = "_PA$JNUM";;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($jname,$tmpCmd) = qsubSystem($logDir."pseudoAssembly.sh",$cmd,1,"5G",$jname,$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;

	} else {$cmd="";}
	my $megDir = $psFile; $megDir =~ s/[^\/]+$//;
	#die "$megDir\n";
	return ($jname,$psFile, $megDir);
 }

#spadesAssembly( \%AsGrps,$cAssGrp,"$nodeSpTmpD/ass",$metagAssDir,$MFopt{spadesBayHam} ,$shortAssembly, $SmplNameX,$hostFilter,$scaffoldFlag) ;
sub spadesAssembly{
	my ($asHr,$cAsGrp,$nodeTmp,$finalOut,$doClean,$helpAssembl,$smplName,$hostFilter,$mateFlag) = @_;

	
	#my $p1ar = $AsGrps{$cAsGrp}{FilterSeq1};my $p2ar = $AsGrps{$cAsGrp}{FilterSeq2};my $singlAr = $AsGrps{$cAsGrp}{FilterSeqS};my $cReadTecAr = $AsGrps{$cAsGrp}{ReadTec};
	my ($p1ar,$p2ar,$singlAr,$cReadTecAr) = getCleanSeqsAssmGrp($asHr, $cAsGrp, 0);
	my $jDepe = $AsGrps{$cAsGrp}{SeqClnDeps};
	
	


	my $spadesBin = getProgPaths("spades");
	my $isCloudSpades = 0;
	foreach my $lRT (@{$cReadTecAr}){if ($lRT =~ m/SLR/){$isCloudSpades = 1;}}
	if ($isCloudSpades){
		$spadesBin = getProgPaths("cloudspades") ;
		print "Using CloudSpades\n";
	}
	

	#print all samples used 
	my $nCores = $MFopt{AssemblyCores};#6
	my $noTmpOnNode = 0; #prevent usage of tmp space on node
	if ($noTmpOnNode){
		$nodeTmp = $finalOut;
	}
 	my $cmd = "rm -rf $nodeTmp\nmkdir -p $nodeTmp\nmkdir -p $finalOut\n\n";
	$cmd .= "\necho '". $AsGrps{$cAsGrp}{AssemblSmplDirs}. "' > $nodeTmp/smpls_used.txt\n\n";
	my $defTotMem = $MFopt{AssemblyMemory};#60;
	if ($defTotMem == -1){ #auto set mem
		$defTotMem = ($inputFileSizeMB{$curSmpl}*4+1e5)/1024;
	}

	my $defMem = ($defTotMem/$nCores);
	$defTotMem = $defMem * $nCores; #total really available mem (in GB)

	$cmd .= $spadesBin;
	my $K = $MFopt{AssemblyKmers} ;
	#insert single reads
	my $errStep = "";
	$errStep = "--only-assembler " if ($doClean == 0);
	my $numInLibs = scalar @{$p1ar};
	my $sprds = inputFmtSpades($p1ar,$p2ar,$singlAr,$logDir,$cReadTecAr);
	#$cmd .= " --meta " ;
	if ($numInLibs <= 1) {$cmd .= " --meta " ;} else {$cmd .= " --sc " ;} #deactivated as never done with other T2 samples..
	$cmd .= " $K $sprds -t $nCores $errStep -m $defTotMem ";#--mismatch-correction "; # --meta  --sc "; #> $log #--meta :buggy in 3.6
	$cmd .= " --mismatch-correction " if ($MFopt{spadesMisMatCor});
	if ($helpAssembl ne ""){
		$cmd .= "--untrusted-contigs $helpAssembl ";
	}
	
	$cmd .= "-o $nodeTmp\n";
	#from here could as well be separate 1 core job
	#cleanup assembly
	$cmd .= "\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector/*\n";
	#dual size filter
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	$cmd .= "$renameCtgScr $nodeTmp/scaffolds.fasta $smplName\n";
	$cmd .= "$sizFiltScr $nodeTmp/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
	
	if ($JNUM > 1 && $helpAssembl ne ""){
		#cluster short reads using CD-HIT, replaces scaffolds.fasta.filt2 file
		my $secAss = $nodeTmp."secondary_shorts/";
		$cmd .= "mkdir -p $secAss\n";
		$cmd .= "cat $nodeTmp/scaffolds.fasta >> $nodeTmp/scaffolds.fasta2\n";
		$cmd .= $spadesBin." -s $nodeTmp/scaffolds.fasta2 --only-assembler -m 1000 -t $nCores $K -o $secAss\n";
		#cleanup
		$cmd .= "cp  $secAss/contigs.fasta $nodeTmp/scaffolds.fasta2\nrm -f -r $secAss\n";
	}
	#die "$cmd\n\n";
	my $assStatScr = getProgPaths("assStat_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/assemblathon_stats.pl";
	$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.500.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.ini.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta.filt > $nodeTmp/AssemblyStats.txt\n";

	my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$nodeTmp/scaffolds.fasta.filt",$nCores,$MFopt{largeMapperDB},$MFopt{MapperProg});#$nCores);
	$cmd .= $cmdDB unless($mateFlag || !$MFopt{map2Assembly}); #doesn't need bowtie index
	
	#clean up
	$cmd .= "\n $pigzBin -f -p $nCores -r $nodeTmp/scaffolds.fasta $nodeTmp/misc/\n $pigzBin -f -p $nCores -r $nodeTmp/contigs.paths  $nodeTmp/*contigs.fa* \n";
	$cmd .=  "rm -rf $nodeTmp/assembly_graph*.gfa $nodeTmp/corrected $nodeTmp/*.fastg $nodeTmp/before_rr*\n";
	unless ($noTmpOnNode){
		$cmd .=  "mkdir -p $finalOut\ncp -r $nodeTmp/* $finalOut\n" ;
		$cmd .= "rm -rf $nodeTmp\n";
	}
	my $jname = "";
	
	if (-e $logDir."spaderun.sh.otxt"){	#check for out of mem
		open I,"<$logDir/spaderun.sh.otxt" or die "Can't open old assembly logfile $logDir\n"; my $str = join("", <I>); close I;
		if ($str =~ / Error in malloc(): out of memory/ ||$str =~ m/TERM_MEMLIMIT: job killed after reaching LSF memory usage limit/){ #memory error for real
			my $replMem  = "";
			if ($str =~ /\n    Max Memory :     (\d+) MB\n/){	$replMem = int($1*1000/$nCores*1.7);
			} elsif ($str =~ /\nMAX MEM (\d+)G\n/){	$replMem = int($1/$nCores*1.7);}
			unless ($replMem eq ""){
				if (($replMem *$nCores)< 50){$replMem = 12;} 
				$defMem = $replMem;
				print $defMem."G: new MEM\n"; #die $defMem."\n";
			}
		}
	}
	$cmd .= "echo \"MAX MEM ".$defTotMem."G\"\n";
	$cmd .= "echo \"SPADES\" > $finalOut/$STOassmbleDone\n";
	#print "in Assembly\n$jDepe\n";
	#print "$finalOut/scaffolds.fasta.filt\n";
	my $locDiskSpace = $HDDspace{assembler};
	if ($locDiskSpace eq "-1"){
		$locDiskSpace = $HDDspace{spades};
	}
	#die "$locDiskSpace\n";
	unless (-e "$finalOut/scaffolds.fasta.filt" && !-z "$finalOut/scaffolds.fasta.filt" && !-z "$finalOut/AssemblyStats.txt"){
		#my $size_in_mb = (-s $fh) / (1024 * 1024);
		my $tmpCmd="";
		$jname = "_A$JNUM";#$givenJName;
		#$QSBoptHR->{useLongQueue} = 1;
		if ($hostFilter || $MFopt{SpadesAlwaysHDDnode}){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};
			$QSBoptHR->{useLongQueue} = $MFopt{SpadesLongtime};
			$QSBoptHR->{tmpSpace} = $locDiskSpace;
			#$QSBoptHR->{tmpSpace} = $HDDspace{spades}; #set option how much tmp space is required, and reset afterwards
			($jname,$tmpCmd) = qsubSystem($logDir."spaderun.sh",$cmd,(int($nCores/2)+1),int($defMem*2)."G",$jname,$jDepe,"",1,$QSBoptHR->{Spades_Hosts},$QSBoptHR) ;
			$QSBoptHR->{tmpSpace} = $tmpSHDD;
			$QSBoptHR->{useLongQueue} = 0;
		} else {
			($jname,$tmpCmd) = qsubSystem($logDir."spaderun.sh",$cmd,(int($nCores/2)+1),int($defMem*2)."G",$jname,$jDepe,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		}
		#$QSBoptHR->{useLongQueue} = 0;
	} else {
		print "Spades: Assembly still on tmp dir\n";
	}
	#die("SPADE\n");
	return ($jname);
}


sub longRdAssembly{
	my ($asHr,$cAsGrp,$nodeTmp,$finalOut,$helpAssembl,$smplName, $useSupportRds) = @_;
	
	
	my ($p1ar,$p2ar,$singlAr,$cReadTecAr) = getCleanSeqsAssmGrp($asHr, $cAsGrp, $useSupportRds);
	if (@{$p1ar} > 0){print "Paired reads defined (@{$p1ar}), but long read assemblies rely on singleton reads!\nAborting\n";die;}
	my $numInLibs = scalar @{$singlAr};
	#print "$numInLibs libs\n";
	
	if ($MFopt{DoAssembly} == 5){#hybrid mode.. check that all required files are present or stop here
		if (${$cReadTecAr}[0] ne "PB"){print "Hybrid Assembly.. expected \"PB\" reads for support reads! (found \
		${$cReadTecAr}[0]\"\nAborting..\n";die;}
		return "" unless (-e $helpAssembl || $helpAssembl eq "hybridmMDBG");
	}
	#my ($p1arX,$p2arX,$singlArX,$cReadTecArX) = getCleanSeqsAssmGrp($asHr, $cAsGrp, 1);

	#my $singlAr = $AsGrps{$cAsGrp}{FilterSeqS};#my $cReadTec = $AsGrps{$cAsGrp}{ReadTec};
	my $jDepe = $AsGrps{$cAsGrp}{SeqClnDeps};
	
	my $nameProg= "flye"; $nameProg="mMDBG"if($MFopt{DoAssembly}==4 || $MFopt{DoAssembly}==5);
	if (${$cReadTecAr}[0] =~ m/SLR/i){die "Can't use synthetic long reads (SLR) with $nameProg\n";}
	my $nCores = $MFopt{AssemblyCores};#6
	
	my $nodeTmp2 = "$nodeTmp/tmpRawRds/";
 	my $cmd = "rm -rf $nodeTmp\nmkdir -p $nodeTmp $finalOut $nodeTmp2\n\n"; #\n  mkdir -p $nodeTmp/tmp\n
	#input reads for assembly .. expected unpaired, long reads
	my @inRds = @{$singlAr};
	
	#metaMDBG "hack" to impute illumina assemblies:
	my $cmdPre = "";
	if ($helpAssembl eq "hybridmMDBG" && $MFopt{DoAssembly} == 5){
		if (${$cReadTecAr}[0] ne "PB"){print"Expected PacBio (\"PB\") readTech for metaMDBG, found \"${$cReadTecAr}[0]\"\n";die;}
		my $spl4m = getProgPaths("split_fasta4metaMDBG_scr");
		#my $illPathS = ;
		my @illDirs = @{$AsGrps{$cAsGrp}{preAsmblDir}}; #split /,/,$illPathS;
		$cmdPre .= "#presplitting helper assembly:\n";
		die "preLib num (" .@illDirs . ") != read libs (" . @inRds . ")!" if (@illDirs != @inRds);
		for (my $i=0;$i<@illDirs;$i++){
			my $illD = $illDirs[$i];
			die "longRdAssembly:: $illD is not a dir!" unless (-d $illD);
			my $contigCov = "$illD/Coverage.median.percontig.gz"; my $dupiAssmbl = "$nodeTmp2/assmbl.$i.pre.fastq.gz";
			my $preAssmbl = "$illD/scaffolds.fasta.filt";
			$cmdPre .= "$spl4m $preAssmbl $contigCov $dupiAssmbl;\n";
			#merge this split with single reads in tmp dir..
			$cmdPre .= "cat $inRds[$i] >> $dupiAssmbl;\n";
			$inRds[$i] = $dupiAssmbl;
		}
		#transfer commands to main #cmd stream..
		$cmd .= $cmdPre."#presplitting done\n\n";
	}
	
	#die "@inRds\n";
	#die $cmd;
	
	
	my $contigRecovery = "";
	if ($MFopt{DoAssembly}==3){#FLYE
		if (${$cReadTecAr}[0] ne "ONT"){print"Expected Oxford Nanopore (\"ONT\") readTech for flye, found \"${$cReadTecAr}[0]\"\n";die;}
		my $flyeBin = getProgPaths("flye");
		$cmd .= $flyeBin;
		$cmd .= " --nano-raw $inRds[0] -t $nCores --meta -g 3g ";
		if ($helpAssembl ne ""){
			$cmd .= "--subassemblies $helpAssembl ";
		}
		$cmd .= "-o $nodeTmp\n";  #--tmp-dir $nodeTmp/tmp/
		#$outAssemblyF = "$nodeTmp/assembly.fasta";
		$contigRecovery .= "\nrm -fr $nodeTmp/00-assembly/ $nodeTmp/10-consensus/ $nodeTmp/20-repeat/ $nodeTmp/30-contigger/ $nodeTmp/40-polishing/ $nodeTmp/assembly_graph.gv\n";
		$contigRecovery .= "\nmv $nodeTmp/assembly.fasta $nodeTmp/scaffolds.fasta\n\n";
	} elsif($MFopt{DoAssembly}==4 || $MFopt{DoAssembly}==5){ #metaMDBG
		my $inFileFlag = "--in-hifi";
		$inFileFlag = "--in-ont" if (${$cReadTecAr}[0] eq "ONT");
		my $mMDBG = getProgPaths("metaMDBG");
		$cmd .= "$mMDBG asm --threads $nCores --out-dir $nodeTmp $inFileFlag " . join(" ",@inRds) . "\n";
		$cmd .= "rm -rf $nodeTmp/tmp/;\n";
		$contigRecovery .= "zcat $nodeTmp/contigs.fasta.gz > $nodeTmp/scaffolds.fasta; rm $nodeTmp/contigs.fasta.gz\n\n";
		
		#contigs.fasta.gz
	} else {
		die "Wrong assembler selected\n";
	}
	$cmd .= "\nrm -rf $nodeTmp2\n";
	
	#$cmd .= getProgPaths("activateBase")."\n"; #not needed any longer.. all assemblers are in base env
	
	
	#from here could as well be separate 1 core job
	#cleanup assembly
	#$cmd .= "\necho '". chomp($AsGrps{$cAsGrp}{AssemblSmplDirs}). "' > $nodeTmp/smpls_used.txt\n\n";
	$cmd .= "\necho '". $AsGrps{$cAsGrp}{AssemblSmplDirs}. "' > $nodeTmp/smpls_used.txt\n\n";

	$cmd .= $contigRecovery;#
	#dual size filter
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	$cmd .= "$renameCtgScr $nodeTmp/scaffolds.fasta $smplName\n";
	$cmd .= "$sizFiltScr $nodeTmp/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
	
	my $assStatScr = getProgPaths("assStat_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/assemblathon_stats.pl";
	$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.500.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.ini.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta.filt > $nodeTmp/AssemblyStats.txt\n";
	
	my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$nodeTmp/scaffolds.fasta.filt",$nCores,$MFopt{largeMapperDB},$MFopt{MapperProg});#$nCores);
	$cmd .= $cmdDB unless( !$MFopt{map2Assembly}); #doesn't need bowtie index
	
	#clean up
	$cmd .= "\n $pigzBin -f -p $nCores $nodeTmp/scaffolds.fasta\n";
	$cmd .=  "mkdir -p $finalOut\ncp -r $nodeTmp/* $finalOut\n" ;
	$cmd .= "rm -rf $nodeTmp\n";
	my $jname = "";

	my $defTotMem = $MFopt{AssemblyMemory};#60;
	if ($defTotMem == -1){ #auto set mem
		$defTotMem = ($inputFileSizeMB{$curSmpl}*8+1e4)/1024;
	}

	my $defMem = ($defTotMem/$nCores);

	$cmd .= "echo \"MAX MEM ".$defTotMem."G\"\n";
	$cmd .= "echo \"$nameProg\" > $finalOut/$STOassmbleDone\n";
	#die "$cmd\n\n";
	if (!-e "$finalOut/scaffolds.fasta.filt"  || !-e "$finalOut/AssemblyStats.txt"){
		#my $size_in_mb = (-s $fh) / (1024 * 1024);
		my $tmpCmd="";
		$jname = "$nameProg$JNUM";#$givenJName;
		$QSBoptHR->{useLongQueue} = 0;#super fast, doesn't need long queue
		if ( $MFopt{SpadesAlwaysHDDnode}){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};
			$QSBoptHR->{tmpSpace} = $HDDspace{riboFind};
			($jname,$tmpCmd) = qsubSystem($logDir."$nameProg.sh",$cmd,(int($nCores)),int($defMem)."G",$jname,$jDepe,"",1,$QSBoptHR->{Spades_Hosts},$QSBoptHR) ;
			$QSBoptHR->{tmpSpace} = $tmpSHDD;
		} else {
			($jname,$tmpCmd) = qsubSystem($logDir."$nameProg.sh",$cmd,(int($nCores)),int($defMem)."G",$jname,$jDepe,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		}
		$QSBoptHR->{useLongQueue} = 0;
	} else {
		print "longReadAssm:Assembly still on tmp dir\n";
	}
	return ($jname);

	
	
}
#( \%AsGrps,$cAssGrp,"$nodeSpTmpD/ass",$metagAssDir ,$shortAssembly, $SmplNameX,$hostFilter,$scaffoldFlag)
sub megahitAssembly{
	my ($asHr,$cAsGrp,$nodeTmp,$finalOut,$helpAssembl,$smplName,$hostFilter,$mateFlag) = @_;
	my $megahitBin = getProgPaths("megahit");
	my $stoneAssmbl = "$STOassmbleDone";
	
#	my $p1ar = $AsGrps{$cAsGrp}{FilterSeq1};my $p2ar = $AsGrps{$cAsGrp}{FilterSeq2};my $singlAr = $AsGrps{$cAsGrp}{FilterSeqS};
	my $jDepe = $AsGrps{$cAsGrp}{SeqClnDeps};
#	my $cReadTec = $AsGrps{$cAsGrp}{ReadTec};
	my ($p1ar,$p2ar,$singlAr,$cReadTec) = getCleanSeqsAssmGrp($asHr, $cAsGrp, 0);

	if (${$cReadTec}[0] =~ m/SLR/i){die "Can't use synthetic long reads (SLR) with megahit\n";}
	#print all samples used 
	my $nCores = $MFopt{AssemblyCores};#6
	my $noTmpOnNode = 0; #prevent usage of tmp space on node
	if ($noTmpOnNode){
		$nodeTmp = $finalOut;
	}
	my $nodePreD = $nodeTmp;$nodePreD=~s/\/[^\/]+\/*$/\//;
	#die "$nodePreD\n$nodeTmp\n";
	
 	my $cmd = "rm -rf $nodeTmp\nmkdir -p $nodePreD $finalOut\n\n"; #\n  mkdir -p $nodeTmp/tmp\n

	my $inputSizeloc = spaceInAssGrp($curSmpl);
	#die "DS$inputSizeloc\n";
	my $defTotMem = $MFopt{AssemblyMemory};#60;
	if ($defTotMem == -1){ #auto set mem
		#die "ts:$inputSloc\n";
		$defTotMem = ($inputSizeloc*1.85 + 5e4)/1024;
	}
	my $defMem = int($defTotMem/$nCores);
	$defTotMem = $defMem * $nCores; #total really available mem (in GB)
	
	my $locDiskSpace = $HDDspace{assembler};
	if ($locDiskSpace eq "-1G" || $locDiskSpace eq "-1"){
		$locDiskSpace = int(25+(4*$inputSizeloc/1024))."G";
		#$locDiskSpace = $HDDspace{megaHit};
	}
	#print "diskspace: $locDiskSpace\n";die;


	my $K = $MFopt{AssemblyKmers} ;
	$K =~ s/-k //;
	#check that k's are closer than 28
	my @spl  = sort {$a <=> $b}(split /,/,$K);
	for (my $i=1;$i<@spl;$i++){
		if ($spl[$i] - $spl[$i-1] > 28){die "kmers for megahit: $K: steps must be <28\n";}
	}
	$K = join(",",@spl);
	#insert single reads
	my $numInLibs = scalar @{$p1ar};
	my $sprds = inputFmtMegahit($p1ar,$p2ar,$singlAr,$logDir);
	$cmd .= $megahitBin;
	$cmd .= " --k-list $K $sprds -t $nCores -m ". int($defTotMem*1024*1024*1024*0.8) ." --out-prefix megaAss ";
	if ($helpAssembl ne ""){
		if ($helpAssembl eq "preAssmbl"){ #check for keywords
			$stoneAssmbl = $STOpreAssmblDone;
			$helpAssembl = "";
		} elsif (-f $helpAssembl) {
			$cmd .= "--untrusted-contigs $helpAssembl ";
		} else { die "Can't decipher helpAssmbl input to megahitAssembly: $helpAssembl\n";}
	}
	$cmd .= "-o $nodeTmp  \n";  #--tmp-dir $nodeTmp/tmp/
	#from here could as well be separate 1 core job
	#cleanup assembly
	$cmd .= "\necho '". $AsGrps{$cAsGrp}{AssemblSmplDirs}. "' > $nodeTmp/smpls_used.txt\n\n";
	$cmd .= "\nrm -fr $nodeTmp/tmp/ $nodeTmp/intermediate_contigs/ \nmv $nodeTmp/megaAss.contigs.fa $nodeTmp/scaffolds.fasta\n\n";
	#dual size filter
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	$cmd .= "$renameCtgScr $nodeTmp/scaffolds.fasta $smplName\n";
	$cmd .= "$sizFiltScr $nodeTmp/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
	
	if ($JNUM > 1 && $helpAssembl ne ""){
		die("helper assemblies not supported with megahit: $helpAssembl\n");
		#cluster short reads using CD-HIT, replaces scaffolds.fasta.filt2 file
		my $secAss = $nodeTmp."secondary_shorts/";
		#$cmd .= "mkdir -p $secAss\n";
		#$cmd .= "cat $nodeTmp/scaffolds.fasta >> $nodeTmp/scaffolds.fasta2\n";
		#$cmd .= $spadesBin." -s $nodeTmp/scaffolds.fasta2 --only-assembler -m 1000 -t $nCores $K -o $secAss\n";
		#cleanup
		#$cmd .= "cp  $secAss/contigs.fasta $nodeTmp/scaffolds.fasta2\nrm -f -r $secAss\n";
	}
	
	
	my $assStatScr = getProgPaths("assStat_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/assemblathon_stats.pl";
	$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.500.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.ini.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta.filt > $nodeTmp/AssemblyStats.txt\n";
	
	my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$nodeTmp/scaffolds.fasta.filt",$nCores,$MFopt{largeMapperDB},$MFopt{MapperProg});#$nCores);
	$cmd .= $cmdDB unless($mateFlag || !$MFopt{map2Assembly}); #doesn't need bowtie index
	#print "DB:: $cmdDB\nMATE::$mateFlag\n";
	#clean up
	$cmd .= "\n $pigzBin -f -p $nCores $nodeTmp/scaffolds.fasta.filt2 $nodeTmp/scaffolds.fasta.lnk\n";
	$cmd .= "rm -f $nodeTmp/scaffolds.fasta\n";
	unless ($noTmpOnNode){
		$cmd .=  "mkdir -p $finalOut\ncp -r $nodeTmp/* $finalOut\n" ;
		$cmd .= "rm -rf $nodeTmp\n";
	}
	my $jname = "";
	
#	if (-e $logDir."megahitrun.sh.otxt"){	#check for out of mem
#	}

	$cmd .= "echo \"MAX MEM ".$defTotMem."G\"\n";
	$cmd .= "echo \"MEGAHIT\" > $finalOut/$stoneAssmbl\n";
	#die "$cmd\n\n";
	if (!-e "$finalOut/scaffolds.fasta.filt"  || !-e "$finalOut/AssemblyStats.txt"){
		#my $size_in_mb = (-s $fh) / (1024 * 1024);
		my $tmpCmd="";
		$jname = "mA$JNUM";#$givenJName;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = $locDiskSpace;#$HDDspace{megaHit};
		($jname,$tmpCmd) = qsubSystem($logDir."megahitrun.sh",$cmd,$nCores,$defMem."G",$jname,$jDepe,"",1,$QSBoptHR->{Spades_Hosts},$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
	} else {
		print "MegaHit: Assembly still on tmp dir\n";
	}
	#die("SPADE\n");
	return ($jname);
}


sub metagAssemblyRun{
	my ( $cAssGrp,$cleanSeqSetHR,$nodeTmp,$metagAssDir ,$shortAssembly, $SmplNameX,$scaffoldFlag,$metaGscaffDir,
				$assemblyFlag,$AssemblyGo,$ePreAssmbly, $doPreAssmFlag, $postAssmblGo) = @_;
				#metagAssembly( $cAssGrp,"$nodeSpTmpD/ass",$metagAssDir ,$shortAssembly, $SmplNameX,$scaffoldFlag,
				#	$assemblyFlag,$AssemblyGo,$ePreAssmbly, $doPreAssmFlag, $postPreAssmblGo);
	print "Assembly step";
	my $hostFilter = 0;$hostFilter = 1 if ($AsGrps{$cAssGrp}{CntAimAss} > 3);#reset required HDD space
	my $tmpN ="";
	my $LasseP = $MFopt{DoAssembly};
	
	if ($LasseP == 5 && $map{$curSmpl}{"SupportReads"} eq ""){#decide on single tech
		$LasseP = 2; #go for megahit by default..
		my ($p1ar,$p2ar,$singlAr,$cReadTecAr) = getCleanSeqsAssmGrp(\%AsGrps, $cAssGrp, 0);
		$LasseP  = 4 if (${$cReadTecAr}[0] eq "PB");
		$LasseP  = 3 if (${$cReadTecAr}[0] eq "ONT");

	}
	
	if ($LasseP == 5){ #hybrid assembly: first megahit(2), then metaMDBG(4)
		if ($doPreAssmFlag == 1){
			print "Preassembly step: ";
			$tmpN = megahitAssembly( \%AsGrps,$cAssGrp,$nodeTmp,$metagAssDir ,
				"preAssmbl", $SmplNameX,$hostFilter,$scaffoldFlag) if (!$ePreAssmbly);
		} elsif ($doPreAssmFlag == 0 && $postAssmblGo) {
			print "Final combining long assembly step: ";
			system "rm -fr $metagAssDir\n"; #just clean up assembly dir..
			#die;
			$tmpN = longRdAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$metagAssDir,
				"hybridmMDBG",$SmplNameX,1) ; #$metaGpreAssmblDir, 
		} elsif ($doPreAssmFlag == 0 ){
			print "Preassembly: waiting for all samples\n";
		} else {
			die "Unknown preassembly state: $doPreAssmFlag\n";
		}
#					die;
	}elsif($LasseP == 1){
		$tmpN = spadesAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$metagAssDir,$MFopt{spadesBayHam} ,
			$shortAssembly, $SmplNameX,$hostFilter,$scaffoldFlag) ;
	}elsif($LasseP == 2){
		$tmpN = megahitAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$metagAssDir ,
			$shortAssembly, $SmplNameX,$hostFilter,$scaffoldFlag) ;
	} elsif( ($LasseP == 3 ||  $LasseP == 4) && ${$cleanSeqSetHR}{is3rdGen} ){
		$tmpN = longRdAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$metagAssDir,
			$shortAssembly, $SmplNameX,0) ;
	}
	#die ;
	#if mates available, do them here
	$AsGrps{$cAssGrp}{AssemblJobName} .= ";".$tmpN; #always add in dep on read extraction
	
	
	# 2nd assembly step: scaffolding; maybe move later further down?
	if ($scaffoldFlag){
	#$finalCommScaffDir "$finalCommScaffDir/scaffDone.sto" $finAssLoc $metaGassembly
		#my $curAssLoc = $metaGassembly;
		#$curAssLoc = $finAssLoc if ($efinAssLoc);
		my ($newScaff,$sdep) = scaffoldCtgs(\%AsGrps,$cAssGrp, #$AsGrps{$cMapGrp}{RawSeq1},$AsGrps{$cMapGrp}{RawSeq2},$AsGrps{$cMapGrp}{Libs},
				[],[],
				$metagAssDir,$nodeTmp."/scaff/",$metaGscaffDir,$AsGrps{$cAssGrp}{AssemblJobName},$MFopt{MapperCores}, $SmplNameX,1,"");
		unless ($sdep eq ""){
			$AsGrps{$cAssGrp}{AssemblJobName} .= ";$sdep";
		}
	}
	
	#external contigs to be scaffolded (e.g. TEC2 extracts)
	if ($scaffTarExternal ne ""){
	#die "inscaff\n";
	#die "@scaffTarExternalOLib1\n";
		my $metaGscaffDirExt = "$metaGscaffDir/$scaffTarExternalName/";
		my ($newScaff,$sdep) = scaffoldCtgs(\%AsGrps,$cAssGrp, #$AsGrps{$cMapGrp}{RawSeq1},$AsGrps{$cMapGrp}{RawSeq2},$AsGrps{$cMapGrp}{Libs},
				#\@scaffTarExternalOLib1,\@scaffTarExternalOLib2,
				[],[],
				$scaffTarExternal,$nodeTmp."/SCFEX$scaffTarExternalName/",$metaGscaffDirExt,$AsGrps{$cAssGrp}{AssemblJobName},$MFopt{MapperCores}, 
				$SmplNameX,0,$scaffTarExternalName);
	#my($ar1,$ar2,$scaffolds,$GFdir_a) = @_; #.= "_GFI1";
		if (@scaffTarExternalOLib1 > 0 ){
			#die "in gapfill\n$newScaff\n";
			GapFillCtgs(\@scaffTarExternalOLib1,\@scaffTarExternalOLib2,$newScaff,$metaGscaffDirExt."GapFill/",$sdep,$scaffTarExternalName);
		}

		#last;
	}

	
	
	return;
} 


sub genePredictions($ $ $ $ $) {
	my ($inputScaff, $outDir, $jobDepend,$finDir,$specJname,$scrathD,$locSubm) = @_;
	my $tmpGene = $outDir."/tmpCalls/";
	my $numThr = 8; #defaul 4 cores..
	#system("mkdir -p $outDir");
	my $scrDir = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/euk_gene_caller/bin";
	my $output_format_prodigal = "gff"; 
	my $expectedD = "$finDir/genePred";
	my $augCmd  = ""; my $cmpCmd = "";
	my $tmpCmd=""; #placeholder for qsub return
	my $splitDep = $jobDepend;
	my $inputBac = $inputScaff; my $inputEuk = "$scrathD/euk.kraken.fasta";
	my $bacmark="";
	if ($MFopt{DoEukGenePred} ){ $bacmark = ".bac";	}
#	die;
	
	my $pprodigalBin = getProgPaths("pprodigal");
	#my $prodigalBin = getProgPaths("prodigal");
	
	#print "$outDir/proteins$bacmark.shrtHD.faa\n";
	
	if ( (-s "$expectedD/proteins$bacmark.shrtHD.faa" && -s "$expectedD/genes$bacmark.gff") ||
			(-s "$outDir/proteins$bacmark.shrtHD.faa" && -s "$outDir/genes$bacmark.gff") ){
				#die "-s $expectedD/proteins$bacmark.shrtHD.faa || -s $expectedD/genes$bacmark.gff";
		return "";
	}
	#use kraken to classify contigs..
	if ($MFopt{DoEukGenePred} ){ #stupid way of doing things..
		#first split euk vs bac contigs :: 22.6.23:: deprecated: use whokaryote that requires prodigal annotations first
		die;
		my $splCores = 10;
		my ($refDB,$shrtDB,$clnCmd) = prepDiamondDB("NOG",$globaldDiaDBdir,$splCores,1);
		my $scmd = "mkdir -p $outDir\n";
		my $splitKgdContig = getProgPaths("contigKgdSplit_scr");
		$scmd .= "$splitKgdContig $inputScaff $scrathD $krakenDBDirGlobal $globaldDiaDBdir $splCores\n";
#		my $splitKgdContig2 = getProgPaths("contigKgdSplit_scr2");
#		$scmd .= "$splitKgdContig2 $inputScaff $scrathD $gff $splCores\n";
		$scmd .= "touch $scrathD/bac.euk.split.sto\n";
		#die $scmd."\n";
		$inputBac = "$scrathD/bact.kraken.fasta";
		#die "$inputEuk || !-e $inputBac\n";
		if (!-e $inputEuk || !-e $inputBac || !-e "$scrathD/bac.euk.split.sto"){
			($splitDep,$tmpCmd) = qsubSystem($logDir."splitAssembly.sh",$scmd,$splCores,int(50/$splCores)."G","_SC$JNUM","$jobDepend;$krakDeps","",1,$QSBoptHR->{General_Hosts},$QSBoptHR); 
		}
	}

#die "toofar";
#run prodigal on bacterial contigs..
	my $prodigal_cmd = "";
	$prodigal_cmd .= "rm -rf $outDir;\nmkdir -p $tmpGene $outDir;\n";
	
	$prodigal_cmd.="if [ -e $inputBac ] && [ ! -s $inputBac ] ; then\n";
	$prodigal_cmd.="touch $tmpGene/proteins$bacmark.shrtHD.faa $tmpGene/genes$bacmark.per.ctg $tmpGene/genes$bacmark.shrtHD.fna $tmpGene/genes$bacmark.gff\n";
	$prodigal_cmd.="else\n";
	
	$prodigal_cmd .= ("$pprodigalBin -i $inputBac -o $tmpGene/genes$bacmark.gff -a $tmpGene/proteins.faa -d $tmpGene/genes.fna -f $output_format_prodigal -p meta -T $numThr \n");
	$prodigal_cmd .= "sleep 1;cut -f1 -d \" \" $tmpGene/proteins.faa > $tmpGene/proteins$bacmark.shrtHD.faa\n" ;
	$prodigal_cmd .= "cut -f1 -d \" \" $tmpGene/genes.fna > $tmpGene/genes$bacmark.shrtHD.fna\n" ;
	$prodigal_cmd .= "cut -f1 $tmpGene/genes$bacmark.gff | sort | uniq -c | grep -v '#' |  awk -v OFS='\\t' {'print \$2, \$1'} > $tmpGene/genes$bacmark.per.ctg\n";
	$prodigal_cmd .= "rm -f $tmpGene/proteins.faa $tmpGene/genes.fna\n";	
	#copy everything to the right place (for now, might have to change this later to take care of Eukarya)
	$prodigal_cmd .= "mkdir -p $expectedD\n";
		
	$prodigal_cmd.="fi\n";
	
	
	
	if ($MFopt{DoEukGenePred} ){ #stupid way of doing things..
		#1st set mark that prodigal genes are on bacteria only
		my $augustusBin = getProgPaths("augustus");#"/g/bork3/home/hildebra/bin/augustus-3.2.1/bin/augustus";
		$prodigal_cmd .= "touch $tmpGene/genes.prodigal.bact.only\n";
		$augCmd .= "mkdir -p $tmpGene\n";
		my $eukmark = ".euk";

		$augCmd.="if [ ! -s $inputEuk ];then\n";
		$augCmd.="touch $tmpGene/genes$eukmark.gff $tmpGene/genes$eukmark.codingseq $tmpGene/genes$eukmark.fna\n";
		$augCmd.="else\n";
		$augCmd .= "$augustusBin --species=ustilago_maydis $inputEuk --protein=off --codingseq=on > $tmpGene/augustus.gff;\n";
		$augCmd .= "perl $scrDir/getAnnoFast.pl --seqfile=$inputEuk $tmpGene/augustus.gff;";
		$augCmd.="mv $tmpGene/augustus.gff $tmpGene/genes$eukmark.gff\n mv $tmpGene/augustus.codingseq $tmpGene/genes$eukmark.codingseq\nmv $tmpGene/augustus.cdsexons $tmpGene/genes$eukmark.fna\n";
		$augCmd .= "cut -f1 -d \" \" $tmpGene/genes$eukmark.fna > $tmpGene/genes$eukmark.shrtHD.fna\n" ;
		$augCmd .= "rm -f $tmpGene/genes$eukmark.fna\n";
		$augCmd.="fi\n";

		
		my $axl = "$tmpGene/augustus.list"; my $pxl = "$tmpGene/genes.list";
		my $cmpCmd = "grep -v \"^#\" $tmpGene/genes.gff > $pxl;\n";
		$cmpCmd .= "grep -w \"transcript\" $tmpGene/augustus.gff > $axl;\n";
		$cmpCmd .= "grep \">\" $inputEuk > $tmpGene.scaff.list;\n";
		$cmpCmd .= "awk \'!/^>/ { printf \"%s\", $0; n = \"\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }\' $tmpGene/augustus.codingseq > $axl.nolines.fa;\n";
		$cmpCmd .= "awk \'!/^>/ { printf \"%s\", $0; n = \"\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }\' $tmpGene/genes.shrtHD.fna > $pxl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/process_nodes.pl $pxl $axl;\n";
		$cmpCmd .= "perl $scrDir/match_nodes_euk.pl $axl.output.txt $axl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/match_nodes_bac.pl $pxl.output.txt $pxl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/overlap_nodes.pl $tmpGene.scaff.list $axl.output.txt $pxl.list.output.txt ;\n"; #$f.blast.matched.txt
		$cmpCmd .= "perl $scrDir/finalgeneparsing.pl $tmpGene.scaff.list.ALLmatched.txt\n";
		#now, let's pull the appropriated gene. for each contig, need to retrieve ALL genes for that contig.
		$cmpCmd .= "cut -f 4 $tmpGene.scaff.list.ALLmatched.txt.decision.txt | sort |uniq> $tmpGene/formatching.txt\n";
		#$cmpCmd .= "perl ./bin/getgenes.pl $tmpGene/formatching.txt $f.aug.nolines.fa.renamed.fa \n";
		#$cmpCmd .= "perl ./bin/getgenes.pl $tmpGene/formatching.txt $f.mgm.nolines.fa.renamed.fa\n";
		#$cmpCmd .= "cat $f.aug.nolines.fa.renamed.fa.output.txt $f.mgm.nolines.fa.renamed.fa.output.txt > $f.genecalled.fa\n";
		#die $augCmd;
	}
	
	#finish up file transfer etc
	$prodigal_cmd .= "mv $tmpGene/* $outDir\nrm -fr $tmpGene\n";

	#if (!-e "$inputScaff") {die "Input scaffold $inputScaff does not exists.\n";}
	#die $prodigal_cmd."\n";
	#if (system $prodigal_cmd){die "Finished prodigal with errors";}
	my $jname = "";
	#print "$expectedD/proteins.shrtHD.faa\n$expectedD/genes.gff\n";
	if ( (!-s "$expectedD/proteins$bacmark.shrtHD.faa" || !-s "$expectedD/genes$bacmark.gff") &&
			(!-s "$outDir/proteins$bacmark.shrtHD.faa" || !-s "$outDir/genes$bacmark.gff") ){
		$jname = "_GP$JNUM"; 
		$jname = $specJname.$JNUM if ($specJname ne "");
		#print "genesubm\n";
		if ($locSubm){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			($jname,$tmpCmd) = qsubSystem($logDir."prodigalrun.sh",$augCmd.$prodigal_cmd,$numThr,"1G",$jname,$splitDep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR); 
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
		} else { #just return run command, can be exectuted from outside
			$jname = $augCmd.$prodigal_cmd."\n";
		}

	}
	return $jname;
}




sub setupHPC{
	my $QSBoptHR1 = emptyQsubOpt($doSubmit,"",$MFconfig{submSytem});
	my $currentJobs = numUserJobs($QSBoptHR1,1);
	print "Found $currentJobs jobs registered to user.";
	#could be only the submitting job is active? is this within a submission?
	if ($currentJobs == 0 && $MFconfig{rmSmplLocks} ==0){print " Auto-removing sample locks (to override set \"-rmSmplLocks -1\").\n";$MFconfig{rmSmplLocks}=1;}
	elsif ($MFconfig{rmSmplLocks} == -1 ){$MFconfig{rmSmplLocks}=0;} #clearly a debug option.. not documented
	else {print "\n";}
	#die;
	#my %QSBopt = %{$QSBoptHR1};
	$QSBoptHR1->{tmpSpace} = $MFconfig{nodeHDDspace};
	$QSBoptHR1->{wcKeysForJob} = $MFconfig{wcKeysForJob};
	$QSBoptHR1->{excludeNodes} = $MFconfig{excludeNodes};
	$QSBoptHR1->{Spades_Hosts} = []; $QSBoptHR1->{General_Hosts} = [];
	spadesHosts();
	#my $LocationCheckStrg=""; #command that is put in front of every qsub, to check if drives are connected, sub checkDrives
	#queing capability
	return $QSBoptHR1
}



sub setDefaultMFconfig{
	
	print "Setting default parameters..  ";
	
	#progStats
	#progStats: object to track progress of programs/submissions
	$progStats{riboFindFailCnts}=0; $progStats{riboFindComplCnts} = 0; 
	$progStats{taxTarFailCnts}=0;  $progStats{taxTarComplCnts}=0; 
	$progStats{mOTU2FailCnts}=0;  $progStats{mOTU2ComplCnts}=0; 
	$progStats{metaPhl2FailCnts}=0; $progStats{metaPhl2ComplCnts}=0;
	$progStats{KrakTaxFailCnts} =0;


	#HDDspace: object to handle HDD usage: Always format as "XXG" XX = space requirements in Gb. Excecption: "-1"
	$HDDspace{kraken} = "120G"; 
	$HDDspace{assembler} = "-1";
	$HDDspace{spades} = "250G"; 
	$HDDspace{megaHit} = "140G";
	$HDDspace{riboFind} = "100G";
	$HDDspace{SNPcall} = "60G";
	$HDDspace{metabat2} = "120G";
	$HDDspace{mapping} = "4G"; #scales linear per input 1 GB filesize
	$HDDspace{diamond} = "80G";
	$HDDspace{prepPub} = "80G";
	$HDDspace{Ribos} = "80G";


	#MFcontstants: object to store essential paths/file endings
	$MFcontstants{bwt2IdxFileSuffix} = ".bw2";
	$MFcontstants{mini2IdxFileSuffix} = ".mmi";
	$MFcontstants{kmaIdxFileSuffix} = ".kma";
	#locking samples during processing (so no other jobs are started in these)
	$MFcontstants{DefaultSampleLock} = "MGTK.locked";



	#MFopt: global object with options for MATAFILER. Added in MF v0.5, slowly rebuild MF around this system

	#non-asembly based tax + functional assignments #DoMetaPhlan3 merges into $DoMetaPhlan , after version check
	#my $DoMetaPhlan = 0;   my $DoMetaPhlan3 = 0; my $DoMOTU2 = 0; my $DoTaxaTarget = 0; 
	$MFopt{DoMetaPhlan}=0;$MFopt{DoMetaPhlan3}=0;
	$MFopt{DoMOTU2}=0;$MFopt{DoTaxaTarget}=0;
	$MFopt{PABtaxChk} =0;


	#my $DoRibofind = 0; my $doRiboAssembl = 0; my $RedoRiboFind = 0; #ITS/SSU/LSU detection
	$MFopt{DoRibofind}=0; $MFopt{doRiboAssembl}=0; $MFopt{RedoRiboFind}=0;#ITS/SSU/LSU detection
	#my $riboLCAmaxRds = 250000; #set to 50k by default, larger gets slower too fast; MF0.45: 250k due to usage of lambda3
	$MFopt{riboLCAmaxRds} = 250000;
	$MFopt{riboStoreRds}=0;$MFopt{RedoRiboAssign}=0;$MFopt{RedoRiboThatFailed}=0;$MFopt{checkRiboNonEmpty}=0;
	$MFopt{globalRiboDependence} = {DBcp => ""};

	#Binning related options
	$MFopt{useBinnerScratch} = 0;$MFopt{BinnerMem} = 0;$MFopt{useCheckM2} = 1;
	$MFopt{DoBinning} = 0;$MFopt{useCheckM1} = 0;$MFopt{BinnerCores} = 9;
	$MFopt{DoMetaBat2} = 0; $MFopt{BinnerRedoEmpty} = 0;

	#read preprocessing
	$MFopt{unzipCores} = 3; 
	$MFopt{trimAdapters} =1;
	$MFopt{sdmOpt} = "";
	$MFopt{usePorechop} = 0;
	$MFopt{useSDM} = 2;
	$MFopt{SDMlogQualvsLen} = 0; #sdm log of qual per read vs length (eg for PacBio qual checks..)
	$MFopt{sdmCores} = 4; #sdm specific cores  #currently set to 1, sdm multi thread instability
	$MFopt{sdm_opt} = {}; #empty object that can be used to modify default sdm parameters
	$MFopt{tmpSdmminSL} =0; $MFopt{tmpSdmmaxSL}=0;
	$MFopt{gzipSDMOut} = 1;#zip sdm filtered files
	$MFopt{sdmProbabilisticFilter} =1;


	#Assembly related options
	$MFopt{doReadMerge} = 0;
	$MFopt{DoAssembly} = 1;  #1=Spades, 2=MegaHIT, 3= flye, 4=metaMDBG, 5=hybrid ill-PB (megahit, metaMDBG)
	$MFopt{SpadesAlwaysHDDnode} = 1;$MFopt{spadesBayHam} = 0; 
	$MFopt{spadesMisMatCor} = 0; $MFopt{redoAssembly} =0 ; $MFopt{SpadesLongtime} = 0;
	$MFopt{pseudoAssembly} = 0; #in case no assembly is possible (soil single reads), just filter for reads X long 
	$MFopt{AssemblyCores} = 8; $MFopt{AssemblyMemory} = -1; #in GB
	$MFopt{AssemblyKmers} = "27,43,67,87,101,127" ; #"27,33,55,71";
	$MFopt{kmerPerGene} = 0; #calculate kmer frequencies for each gene instead of per scaffold
	$MFopt{kmerAssembly} = 0; #calculate kmer frequencies for each scaffold
	$MFopt{scaffoldMinSize} = 500; #all scaffolds/contigs below this will be dropped

	#mapping related options
	$MFopt{MapperProg} = -1;#1=bowtie2, 2=bwa, 3=minimap2, 4=kma, 5=strobealign, -1=auto (bowtie2 short, minimap2 long reads)
	$MFopt{map2Assembly} = 1; $MFopt{mapSortMemGb} = 20; #in Gb
	$MFopt{SaveUnalignedReads} =0; $MFopt{useUnmapped} = 0;
	$MFopt{mapSupport2Assembly} = 0;
	$MFopt{bwtIdxAssMem} = 40; #total mem in GB, not core adjusted, for building index from assembly
	$MFopt{doBam2Cram} = 1; $MFopt{redoAssMapping} =0;
	$MFopt{DoJGIcoverage} = 0; #only required for metabat binning, not required any longer..
	$MFopt{bamfilterIll} = "0.05 0.75 20"; $MFopt{bamfilterPB} = "0.15 0.75 10";
	$MFopt{mapSaveCram} = 1; #by default, keep the back-mapping bams/crams 
	$MFopt{MapperCores} = 8;  $MFopt{MapperRmDup} = 1; #mapping cores; ??? ; remove Dups (can be costly if many ref seqs present)
	$MFopt{bamSortCores} = -1;
	$MFopt{MapperMemory} = -1; #total mem for mapping job in Gb, Default -1
	$MFopt{redoMapping} = 0; #rewrite mapping?
	$MFopt{mapModeActive}=0; $MFopt{mapModeCovDo}=1;#get the coverage per gene etc
	$MFopt{mapModeTogether} = -1; #-1: comb map & report, 0: separate bwt runs, 1:competitive, 2: non-competitive mapping, sep reports per file
	$MFopt{largeMapperDB} = 0; #use flags in mapper index built for large ref DBs?

	$MFopt{DoMapModeDecoy} =1;  $MFopt{MapRewrite2nd} = 0;
	$MFopt{Do2ndMapSNP} = 0;
	$MFopt{refDBall} =""; $MFopt{bwt2NameAll}="" ; #user options for refDBs

	#extra programs to be used?
	$MFopt{DoNonPareil}=0; #non-pareil
	$MFopt{DoCalcD2s} = 0;  #D2s distance .. bit outdated, no longer used..
	$MFopt{calcOrthoPlacement} =0; #Jaime's 6 frame translation and hmmersearch (AB production)
	$MFopt{DoGenoSizeEst} =0; #genome size estimate via microcensus

	#Kraken related config
	$MFopt{DoKraken} = 0; $MFopt{RedoKraken} = 0; 
	$MFopt{krakenCores} = 9;
	$MFopt{completeContaStats} = 1;
	$MFopt{humanFilter} = 0; ##0: no, 1:kraken2, 2: kraken1, 3:hostile
	$MFopt{filterHostDB1} = "";  #customize host org to filter (e.g. human, chicken ..)
	$MFopt{krakHostConf} = 0.01; #confidence needed to assign read to host ref
	$MFopt{filterHostKr2QuickMode} = "";# "--quick "; deactivated for now..
	$MFopt{hostileIndex} = "human-t2t-hla";
	$MFopt{globalKraTaxkDB} = "";
	$MFopt{globalDiamondDependence} = {CZy=>"",MOH2 => "", MOH=>"",NOG=>"",ABR=>"",ABRc=>"",KGB=>"",KGE=>"",ACL=>"",KGM=>"", PTV=>"", PAB => ""};
	



	#SNPs
	$MFopt{DoConsSNP}=0; $MFopt{DoSuppConsSNP}=0; $MFopt{redoSNPcons} = 0; $MFopt{redoSNPgene} =0; $MFopt{SNPconsJobsPsmpl} = 1; 
	$MFopt{SNPminCallQual} = 20;
	$MFopt{saveVCF} = 0; $MFopt{maxSNPcores} = 5; $MFopt{memSNPcall} = 23; $MFopt{consSNPminDepth} = 0;
	$MFopt{SNPcallerFlag} = "MPI"; #"MPI" mpileup or ".FB" for freebayes

	#Func annotation
	$MFopt{DoDiamond} = 0; $MFopt{rewriteDiamond} =0; $MFopt{redoDiamondParse} = 0; #redoes matching of reads; redoes interpretation
	$MFopt{rewriteAllIfAnyDiamond}=0;
	$MFopt{maxReqDiaDB} = 6; #max number of databases supported by METAFILER
	$MFopt{reqDiaDB} = "";#,NOG,MOH,ABR,ABRc,ACL,KGM,PTV,PAB";#,ACL,KGM,ABRc,CZy";#"NOG,CZy"; #"NOG,MOH,CZy,ABR,ABRc,ACL,KGM"   #old KGE,KGB
	$MFopt{diaEVal} = "1e-7"; $MFopt{diaCores} = 12; ; $MFopt{DiaRmRawHits} = 0; $MFopt{diaRunSensitive} = 0;
	$MFopt{DiaMinAlignLen} = 20; $MFopt{DiaMinFracQueryCov} = 0.1; $MFopt{DiaPercID} =40;

	#gene prediction related
	$MFopt{DoEukGenePred} = 0;
	
	
	#MFconfig configuration with defaults
	$MFconfig{mapFile} = "";
	$MFconfig{configFile} = "";
	$MFconfig{nodeHDDspace} = 30; #30 Gb; default value
	$MFconfig{maxUnzpJobs} = 20; #how many unzip jobs to run in parallel (not to overload HPC IO)

	$MFconfig{rawFileSrchStr1} = '.*1\.f[^\.]*q\.gz$';
	$MFconfig{rawFileSrchStr2} = '.*2\.f[^\.]*q\.gz$';
	$MFconfig{rawFileBamSrchSing} = ""; #inputBAMregex 
	$MFconfig{rawFileSrchStrSingl} = "";
	$MFconfig{rawFileSrchStrXtra1} = '.*1_sequence\.f[^\.]*q\.gz$';
	$MFconfig{rawFileSrchStrXtra2} = '.*2_sequence\.f[^\.]*q\.gz$';
	$MFconfig{submSytem} = ""; #user supplied submission system flag.. default is empty and autodetect


	$MFconfig{mateInsertLength} = 20000; #controls expected mate insert size , import for bowtie2 mappings


	#more specific control: unfiniRew=rewrite unfinished sample dir; $redoCS = redo ContigStats completely; 
	#removeInputAgain=remove unzipped files from scratch, after sdm; remove_reads_tmpDir = leave cleaned reads on scratch after everything finishes
	$MFconfig{unfiniRew}=0; $MFconfig{redoCS}=0; $MFconfig{removeInputAgain}=1; $MFconfig{remove_reads_tmpDir}=1; 
	$MFconfig{OKtoRWassGrps} = 0; $MFconfig{rmBinFailAssmbly} = 0;
	$MFconfig{skipWrongPairedSmpls} =1; #check in sdm output if wrong read pairs present
	$MFconfig{defaultContigSubs} = "as"; #default subprts for contigstats..

	$MFconfig{silent} = 0;
	$MFconfig{redoFails} = 0;$MFconfig{XfirstReads} = -1;
	$MFconfig{killDepNever} = 0;  $MFconfig{checkMaxNumJobs} = 0;  #slurm related.. $killDepNever=1 kills jobs in state "DependencyNeverFinished" (happens a lot), while $checkMaxNumJobs=X halts the pipeline if more than X jobs are already queued up
	$MFconfig{excludeNodes} = ""; #excluding certain nodes..
	$MFconfig{readsRpairs} =-1; #are reads given in pairs? default: -1 = no clue
	#my $useTrimomatic=0; #trimmomatic step now replaced by sdm solution -> $MFopt{trimAdapters}

	$MFconfig{abortOnEmptyInput} = 0; 
	#my $relaxedSmplNames = 0;#don't abort when SMPLID in map contains strange chars  -> this is now handled by adding to map "#RelaxSMPLID	TRUE"
	$MFconfig{ignoreSmpl} = "";
	$MFconfig{splitFastaInput} = 0; #assembly as input..
	$MFconfig{importMocat} = 0; $MFconfig{mocatFiltPath} = "reads.screened.screened.adapter.on.hg19.solexaqa/"; 
	$MFconfig{alwaysDoStats} = 1; 
	$MFconfig{rmScratchTmp}=1;#Default; extremely important option as this adds a lot of overhead and scratch usage space, but reduces later overhead a lot and makes IO more stable
	$MFconfig{unpackZip} = 0; #only goes through with read filtering, needed to get files for luis
	$MFconfig{filterFromSource}=1; #powerful option that skips the unzip step.. use careful
	$MFconfig{doDateFileCheck} = 0; #very specific option for Moh's reads that were of different dates..
	$MFconfig{DoFreeGlbTmp} = 0; 
	$MFconfig{defaultReadLength} = 150; $MFconfig{defaultReadLengthX} = 5000;
	$MFconfig{oldStylFolders} =0; #0=smpl name as out folder; 1=inputdir as out foler (legacy)
	$MFconfig{mocatLinkDir} = "";
	$MFconfig{wcKeysForJob} = ""; #EI specific system to register jobs under certain flag
	$MFconfig{prefSinglFQgreps} = 0; #if grep of files (rawSrchString) has multi assignments, which grep to trust more?
	$MFconfig{rmSmplLocks} = 0;
	$MFconfig{uploadRawRds} = ""; #prepare raw input fastq's for upload 2 EBI? Clean reads will be stored in this dir

	
	print "Done. ";

}



sub help {
	print "Help for MG-TK version $MATFILER_ver\n";
	print "TODO :( ref to README.md for now..\n\n";
	exit(0);
}



sub getCmdLineOptions{
	
	print "Reading command line options.. ";
	
	GetOptions(

	#base options
		"help|?|h" => \&help,
		"checkInstall" => sub { checkMFFInstall("",1) },
		"map=s"      => \$MFconfig{mapFile},
		"config=s" => \$MFconfig{configFile},

	#flow related
		"redoFails=i" =>\$MFconfig{redoFails}, #if any step of requested analysis failed, just redo everything (extraction etc)
		"redoContigStats=i" => \$MFconfig{redoCS}, #runContigStats (coverage per gene, kmers, GC content) will be deleted & started again
		"submSystem=s" => \$MFconfig{submSytem},  #qsub,SGE,bsub,LSF.. by default will try to autodetect
		"submit=i" => \$doSubmit,  #submit any jobs at all? (0= no submission, just for trying if everything is correctly set up)
		"from=i" => \$FROM1,  #start at which samples from map file?
		"to=i" => \$TO1,   #stop at which samples from map file?
		"loopTillComplete=s" => \$loop2completion, #dangerous flag, script will loop over the assigned samples until all jobs are finished.
		#use synatx "X:Y" where X is num loops, Y is the window size, eg "6:250" would run 6 loops of max 250 samples, then move on to next 250 samples
		"excludeNodes=s" => \$MFconfig{excludeNodes}, #exclude certain nodes?
		"maxConcurrentJobs=i" => \$MFconfig{checkMaxNumJobs}, #max jobs in queue, useful for large samples sets, currently only works on slurm 
		"killDepNever=i" => \$MFconfig{killDepNever}, #kill jobs in "Dependency never finished" state? 
		"requireInput=i" => \$MFconfig{abortOnEmptyInput},  #in case input reads are no longer present, 0 will continue pipeline, 1 will abort
		"ignoreSmpls=s" => \$MFconfig{ignoreSmpl},  #skip a certain sample (sample id)
		"rmSmplLocks=i" => \$MFconfig{rmSmplLocks}, #remove existing sample locks (useful if jobs have crashed, leaving abondened sample locks)
		#"rmRawRds=i" => \$MFconfig{DoFreeGlbTmp}, #rm raw sequences, once all jobs have finished <-- redundant with reduceScratchUse
		"silent" => \$MFconfig{silent},
		"maxUnzpJobs=i" => \$MFconfig{maxUnzpJobs}, #how many unzip jobs to run in parallel (not to overload HPC IO). Default:20
		
	#input FQ related

	#file strucuture
		#"relaxedSmplNames=i" => \$relaxedSmplNames, #add instead to map: #RelaxSMPLID	TRUE
		"rm_tmpdir_reads=i" => \$MFconfig{remove_reads_tmpDir}, #Default 1, remove tmpdir with reads
		"rm_tmpInput=i" => \$MFconfig{removeInputAgain},#remove raw, human / adaptor filtered reads, if sdm clean created? (and not needed any longer)
		"reduceScratchUse=i" => \$MFconfig{rmScratchTmp}, #should always be 1, unless debugging..
		"globalTmpDir=s" => \$sharedTmpDirP,#absolute path to global shared tmp dir (like a scratch dir)
		"nodeTmpDir=s" => \$nodeTmpDirBase,#absolute path to tmp dir on local HDD of each executing node
		"nodeHDDspace=s" => \$MFconfig{nodeHDDspace},#HDD tmp space to be requested for each node (in Gb). Some systems don't support this
		"legacyFolders=i" => \$MFconfig{oldStylFolders}, #legacy option, controls if output folders will use the read dir as name (1) or the name in the mapping file (0). Default=0

	#preprocessing (cleaning reads etc)
	#	"useTrimomatic=i" => \$useTrimomatic, #rm adapter seq from input reads
		"usePorechop=i" => \$MFopt{usePorechop}, #adapter rm for Nanopore.. should actually be automatically with newer sdm (not implemented)
		"inputFQregex1=s" => \$MFconfig{rawFileSrchStr1}, #regex for detecting read pair 1 in input fastq files
		"inputFQregex2=s" => \$MFconfig{rawFileSrchStr2}, #regex for detecting read pair 2 in input fastq files
		"inputFQregexSingle=s" => \$MFconfig{rawFileSrchStrSingl}, #regex for detecting single end reads in input fastq files
		"inputFQregexTrustSingle=i" => \$MFconfig{prefSinglFQgreps} , #if grep of files (rawSrchString) has multi assignments, which grep to trust more?
		"inputBAMregex=s" => \$MFconfig{rawFileBamSrchSing}, #bams that will be converted to fastq
		"splitFastaInput=i" => \$MFconfig{splitFastaInput},
		"mergeReads=i" => \$MFopt{doReadMerge},  #merge read pair 1+2 before assembly etc? (usually doesn't help assembly, but useful for mapping to ref database in some rare instances)
		"ProbRdFilter=i" => \$MFopt{sdmProbabilisticFilter},
		"pairedReadInput=i" => \$MFconfig{readsRpairs}, #determines if read pairs are expected in each in dir
		"inputReadLength=i" => \$MFconfig{defaultReadLength},
		"inputReadLengthSuppl=i" => \$MFconfig{defaultReadLengthX},
		"filterHostRds|filterHumanRds=i" => \$MFopt{humanFilter}, #0: no, 1:kraken2, 2: kraken1, 3:hostile
		"filterHostKrak2DB=s" => \ $MFopt{filterHostDB1}, #customize host org to filter (e.g. human, chicken ..)
		"filterHostKr2Conf=s" => \$MFopt{krakHostConf},
		"filterHostKr2Quick=s" => \$MFopt{filterHostKr2QuickMode},
		"hostileIndex=s" => \$MFopt{hostileIndex},
		"onlyFilterZip=i" => \$MFconfig{unpackZip},
		"mocatFiltered=i" => \$MFconfig{importMocat},
		"logQualvsLen=i" => \$MFopt{SDMlogQualvsLen}, #sdm log file.. can be quite large; logs qual of read vs read length

	#sdm related
		"gzipSDMout=i" => \$MFopt{gzipSDMOut},
		"XfirstReads=i" => \$MFconfig{XfirstReads},
		"minReadLength=i" => \$MFopt{tmpSdmminSL},
		"maxReadLength=i" => \$MFopt{tmpSdmmaxSL},
		"filterAdapters=i" => \$MFopt{trimAdapters},
		"customSDMopt=s"  => \$MFopt{sdmOpt},

	#assembly related
		"spadesCores|assemblCores=i" => \$MFopt{AssemblyCores},
		"spadesMemory|assemblMemory=i" => \$MFopt{AssemblyMemory}, #in GB
		"spadesKmers|assemblyKmers=s" => \$MFopt{AssemblyKmers}, #comma delimited list
		"reAssembleMG=i" => \$MFopt{redoAssembly},
		"asssemblyHddSpace=i" => \$HDDspace{assembler},
		"assembleMG=i" => \$MFopt{DoAssembly}, #1=Spades, 2=MegaHIT, 3= flye, 4=metaMDBG, 5=hybrid ill-PB (megahit, metaMDBG)
		"assemblyLongTime=i" => \$MFopt{SpadesLongtime},
		"assemblyScaffMinSize=i" => \$MFopt{scaffoldMinSize},
	#binning
		"Binner|MetaBat2|binSpeciesMG=i" => \$MFopt{DoMetaBat2}, #0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder
		"BinnerCores=i" => \$MFopt{BinnerCores}, #cores used for Binning process (and checkM)
		"BinnerMem=i" => \$MFopt{BinnerMem}, # define binning memory, Gb
		"checkM2=i" => \$MFopt{useCheckM2},
		"checkM1=i" => \$MFopt{useCheckM1},
		"BinnerScratchTmp=i" => \$MFopt{useBinnerScratch}, #very specific (undocumented) use of scratch instead of nodetmp dir
		#"binSpeciesMG=i" => \$MFopt{DoBinning},#deactivated, replaced by MetaBat2
		"redoEmptyBins=i" => \$MFopt{BinnerRedoEmpty}, #debug option; redo bins that are empty (no bin detected). Note: this can sometimes happen for metagenomes
		"redoBinning=i" => \$MFopt{BinnerRedoAll}, 
	#gene prediction on assembly
		"predictEukGenes=i" => \$MFopt{DoEukGenePred},#severely limits total predicted gene amount (~25% of total genes)
		"kmerPerGene=i" => \$MFopt{kmerPerGene}, #calculate kmer frequencies for each gene instead of per scaffold
	#mapping
		"mapper=i" => \$MFopt{MapperProg}, ##1=bowtie2, 2=bwa, 3=minimap2, 4=kma, 5=strobealign -1=auto (bowtie2 short, minimap2 long reads), -2=auto(strobealign short, minimap2 long) 
		"mapUnmapped=i" => \$MFopt{useUnmapped},
		"mappingCoverage=i" => \$MFopt{mapModeCovDo},
		"mappingMem=i" => \$MFopt{MapperMemory}, #total mem for mini2/kma/bwa/bwt2 in GB
		"mapSortMem=i" => \$MFopt{mapSortMemGb}, #total mem for samtools sort in GB
		"rmDuplicates=i" => \$MFopt{MapperRmDup},
		"mappingCores=i" => \$MFopt{MapperCores},
		"mapperFilterIll=s" => \$MFopt{bamfilterIll}, #defaults to "0.05 0.75 20", meaning: <=5% ANI, >=75% of read aligned, >=20 mapping quality
		"mapperFilterPB=s" => \$MFopt{bamfilterPB},
		"mapSaveCRAM=i" => \$MFopt{mapSaveCram},
		#"redoMapping=i" =>\$MFopt{redoMapping},
	#mapping related (2) (assembly)
		"remap2assembly|redoMap2assembly|redoMapping=i" => \$MFopt{redoAssMapping},
		"JGIdepths=i" => \$MFopt{DoJGIcoverage},
		"mapReadsOntoAssembly=i" => \$MFopt{map2Assembly} ,  #map original reads back on assembly, to estimate abundance etc
		"mapSupportReadsOntoAssembly=i" => \$MFopt{mapSupport2Assembly}, # (1) map "SupportReads" onto assembly. Default: 0
		"saveReadsNotMap2Assembly=i" => \$MFopt{SaveUnalignedReads},
	#map2tar / map2DB / map2GC
		"decoyMapping=i" => \$MFopt{DoMapModeDecoy},	#1: "Decoy mapping": map against reference genome AND against assembly of metagenome (drawing obvious better hits to metagenome, the "decoy")
		"competitive2ndmap=i" => \$MFopt{mapModeTogether}, #1: Competitive, 2: combined but report separately per input genome, -1: combined and report all together 
		"ref=s" => \$MFopt{refDBall},
		"mapperLargeRef=i" => \$MFopt{largeMapperDB}, #use flags in mapper index built for large ref DBs?
		"mapnms=s" => \$MFopt{bwt2NameAll}, #name for this final files
		"redo2ndmap=i" => \$MFopt{MapRewrite2nd},
	#SNPs 
		"get2ndMappingConsSNP=i" => \$MFopt{Do2ndMapSNP},#SNPs (onto mapping)
		"getAssemblConsSNP=i" => \$MFopt{DoConsSNP},  #SNPs (onto self assembly) #calculates consensus SNP of assembly (useful for checking assembly gets consensus and Assmbl_grps)
		"getAssemblConsSNPsuppRds=i" => \$MFopt{DoSuppConsSNP}, #same as getAssemblConsSNP, but SNP calling for support reads
		"redoAssmblConsSNP=i" => \$MFopt{redoSNPcons},
		"redoGeneExtrSNP=i" => \$MFopt{redoSNPgene},
		"SNPjobSsplit=i" => \$MFopt{SNPconsJobsPsmpl}, #how many parallel jobs are run on each 
		"SNPminCallQual=i" => \$MFopt{SNPminCallQual},
		"SNPsaveVCF=i" => \$MFopt{saveVCF},
		"SNPcaller=s" => \$MFopt{SNPcallerFlag},
		"SNPcores=i" => \$MFopt{maxSNPcores},
		"SNPmem=i" => \$MFopt{memSNPcall}, #memory for consensus SNP job in Gb
		"SNPconsMinDepth=i" => \$MFopt{consSNPminDepth}, #how many reads coverage to include position for consensus call?
		
	#functional profiling (diamond)
		"profileFunct=i"=> \$MFopt{DoDiamond},
		"reParseFunct=i" => \$MFopt{redoDiamondParse},
		"reProfileFunct=i" => \$MFopt{rewriteDiamond},
		"reProfileFuncTogether=i" => \$MFopt{rewriteAllIfAnyDiamond}, #if any func database needs to be redone, than redo all indicated databases (useful if number of reads used changes..)
		"diamondCores=i" => \$MFopt{diaCores},
		"DiaParseEvals=s" => \$MFopt{diaEVal}, #evalues at which to accept hits to func database
		"DiaSensitiveMode=i" => \$MFopt{diaRunSensitive},
		"rmRawDiamondHits=i" => \$MFopt{DiaRmRawHits},
		"DiaMinAlignLen=i" => \$MFopt{DiaMinAlignLen},
		"DiaMinFracQueryCov=f" =>  \$MFopt{DiaMinFracQueryCov},
		"DiaPercID=i" => \$MFopt{DiaPercID},
		"diamondDBs=s" => \$MFopt{reqDiaDB},#NOG,MOH,ABR,ABRc,ACL,KGM,CZy,PTV,PAB,MOH2
	#functional profiling (Jaime tree)
		"orthoExtract=i" => \$MFopt{calcOrthoPlacement},
	#ribo profiling (miTag)
		"profileRibosome=i" => \$MFopt{DoRibofind},
		"riobsomalAssembly=i"  => \$MFopt{doRiboAssembl},
		"reProfileRibosome=i" => \$MFopt{RedoRiboFind} ,  
		"reRibosomeLCA=i"=> \$MFopt{RedoRiboAssign},
		"riboMaxRds=i" => \$MFopt{riboLCAmaxRds},
		"saveRiboRds=i" => \$MFopt{riboStoreRds},
		"thoroughCheckRiboFinish=i" => \$MFopt{checkRiboNonEmpty},
	#other tax profilers..
		"profileMetaphlan2=i"=> \$MFopt{DoMetaPhlan},
		"profileMetaphlan3=i"=> \$MFopt{DoMetaPhlan3},
		"profileMOTU2=i" => \$MFopt{DoMOTU2},
		"profileKraken=i"=> \$MFopt{DoKraken},
		"profileTaxaTarget=i" => \$MFopt{DoTaxaTarget},
		"estGenoSize=i" => \$MFopt{DoGenoSizeEst}, #estimate average size of genomes in data
		"krakenDB=s"=> \$MFopt{globalKraTaxkDB}, #"virusDB";#= "minikraken_2015/";
	#D2s distance
		"calcInterMGdistance=i" => \$MFopt{DoCalcD2s},
	#IO for specific uses
		"newFileStructure=s" => \$MFconfig{mocatLinkDir},#just relink raw files for use in mocat
		"upload2EBI=s" => \$MFconfig{uploadRawRds}, #copy human read removed raw files to this dir, named after sample
	#institute specific: EI
		"wcKeyJobs=s" => \$MFconfig{wcKeysForJob},
	#DEBUG
		"OKtoRWassGrps=i" => \$MFconfig{OKtoRWassGrps}, # can delete assemblies, if suspects error in them
	);
	
	
	# ------------------------------------------ options post processing ------------------------------------------
	setConfigFile($MFconfig{configFile});

	die "No mapping file provided (-map)\n" if ($MFconfig{mapFile} eq "");
	if (!$MFopt{DoAssembly}){
		$MFopt{mapSupport2Assembly}=0;$MFopt{map2Assembly}=0;
	}

	die "\"-mappingMem\" argument contains characters: $MFopt{MapperMemory}" if ($MFopt{MapperMemory} !~ m/[\d-]+/);
	die "\"-mapSortMem\" argument contains characters: $MFopt{mapSortMemGb}" if ($MFopt{mapSortMemGb} !~ m/[\d-]+/);
	die "\"-assemblMemory\" argument contains characters: $MFopt{AssemblyMemory}" if ($MFopt{AssemblyMemory} !~ m/[\d-]+/);
	die "\"-BinnerMem\" argument contains characters: $MFopt{BinnerMem}" if ($MFopt{BinnerMem}  !~ m/[\d-]+/);
	die "\"-SNPmem\" argument contains characters: $MFopt{memSNPcall}" if ($MFopt{memSNPcall} !~ m/[\d-]+/);
	if ($MFopt{MapperMemory} == -1 ){
		if($MFopt{MapperProg} >2){$MFopt{MapperMemory} = 35 ;
		} else {$MFopt{MapperMemory} = 20 ;}	
	}

	if ($MFopt{DoDiamond} && $MFopt{reqDiaDB} eq ""){die "Functional profiling was requested (-profileFunct 1), but no DB to map against was defined (-diamondDBs)\n";}
	$MFopt{AssemblyKmers} = "-k $MFopt{AssemblyKmers}" unless ($MFopt{AssemblyKmers} =~ m/^-k/);
	$MFopt{sdm_opt}->{minSeqLength}=$MFopt{tmpSdmminSL} if ($MFopt{tmpSdmminSL} > 0);
	$MFopt{sdm_opt}->{maxSeqLength}=$MFopt{tmpSdmmaxSL} if ($MFopt{tmpSdmmaxSL} > 0);
	$MFconfig{remove_reads_tmpDir} = 1 if ($MFconfig{DoFreeGlbTmp} || $MFconfig{rmScratchTmp});
	if ($MFopt{bamSortCores} == -1){$MFopt{bamSortCores} = $MFopt{MapperCores};}
	$MFconfig{filterFromSource}=1 if ($MFconfig{unpackZip} );
	@filterHostDB = split /,/,$MFopt{filterHostDB1};
	die "SNPcaller argument invalid, has to be \"MPI\" or \"FB\"\n" if ($MFopt{SNPcallerFlag} ne "MPI" && $MFopt{SNPcallerFlag} ne "FB");
	foreach my $k (keys (%HDDspace)){
		$HDDspace{$k} .= "G" unless ($HDDspace{$k} =~ m/G$/);
	}
	
	if ($MFopt{DoAssembly} == 5 && $MFopt{mapSaveCram}){
		print "deactivating \"-mapSaveCram\", not supported for hybrid assemblies\n";
		$MFopt{mapSaveCram} = 0;
	}
	
	
	print "Done. ";

}


