#!/usr/bin/perl
#This will build a phylo tree for each MGS (between samples)
#./strain_within.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//MB2.clusters.ext.can.Rhcl.mgs 10 /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//between_phylo/prunned.nwk 1 1
#./strain_within.pl /g/bork3/home/hildebra/data/SNP/GCs/alienGC2 /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Binning/MetaBat/MB2.clusters.obs 12 
use warnings;
use strict;

use Getopt::Long qw( GetOptions );
use List::Util qw/shuffle/;


use Mods::GenoMetaAss qw(fileGZe fileGZs readClstrRev systemW median mean readMapS readFasta getAssemblPath getAssemblGFF getAssemblContigs checkSeqTech);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystem2 qsubSystemJobAlive qsubSystemWaitMaxJobs);
use Mods::IO_Tamoc_progs qw(getProgPaths truePath);
use Mods::TamocFunc qw(readTabbed getFileStr checkMF);
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::math qw(quantileArray);

sub extractFNAFAA2genes;
sub histoMGS;
sub readGenesSample_Singl;
sub reportingsMGS;
sub prepRun;
sub filterMultiCopyGenes;
sub createAGlist; sub preComputeConsSNP;
sub timeNice;


die "Not enough args!\n" unless (@ARGV > 1);

#v.14: reworked massively how many genes get included
#v.15: included lessons learned from MGS.pl v0.21 
#v.16 added familyVar and groupStabilityVars arguments for stability calculations
#v.17: considerations to improve speed of intial fna/faa extractions..
#v.18: 16.11.24: handling genes occurring >1 in a single sample/assembly
#v.19: 17.11.24: stricter filtering of genes, removing entire MGS if too many "bad genes" in them; added abundance based filtering of genes/sample
#v0.20: 22.11.24: fixed bug with v0.19 no longer accepting assmblGrps. code refactor that makes it a lot easier to understand
#v0.21: 2.1.25: v0.20 fix, to only select single gene instead of COG; further changed how genes are selected, to reomve potentially conspecific MGS per sample (instead of removing entire gene)
#v0.22: added per sample (not assmblGrp) MGS filtering based on multigenes
#v.23: removed MGS conspecific filter: was too harsh and didn't make sense to have a global filter: MGS are conspecific in a single sample, not all samples..
#.24: 31.10.25: on-the-fly creation of SNP consensus fastas, if correct vcf present
#.25: 22.12.25: precompute for vcf2fna added
#.26: 28.12.25: code refactor to later enable parallelization of main gene-collecting routine
#.27: 12.2.26: made code faster and more stable. changed default MSA aligner
my $version = 0.27;

#input args..
my $GCd = "";#$ARGV[0];
my $MGSfile = "";#$ARGV[1];
my $numCores = 4;#$ARGV[2];
my $subJob=0;#if 0, is main submitting job..
my $maxSubJob = 0;#into how many subjobs to split??
my $outDpre = "";
my $locTmpDir = ""; my $locTmpDir1 = "";
my $maxCores = -1;
my $onlySubmit =0;#extract genes anew?
my $reSubmit=0;
my $treeFile = "";
my $doSubmit=0;
my $subMode="";
my $multiGeneSmplMax = 0.15; #no higher than this rate in single samples conspec genes..
my $conspGeneSmplMax = 0.05; #no higher than this conspecific genes/MGS/sample

my $maxNGenes = 400;
my $MSAprog = 2; ##(0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5
my $presortGenes = 1200;
my $checkMaxNumJobs = 600;
my $useGTDBmg = "GTDB";
my $selfMemGb = 10;
my $redoSubmissionData = 0;
my $PhylosExist = 1;
my $deepRepair = 0;
my $rmMSA = 1; #argument passed to buildTree5.pl 
my $contTests = ""; my $discTests = ""; #stat tests to be given to strain_within_2.2.pl
my $familyVar = ""; my $groupStabilityVars = "";

#SNP calling
my $minSNPDepth = 2; #changed to two: seems to give better results
my $minSNPCallQual = 20;
my $forceVCF2FNA = 0; #force the recalc of cons fasta from vcf..
my $SNPconsLOGs = ""; #logs for recalculating cons SNPs
my $preCompCons=0; #if >0, precompute in these blocks

my $takeAll = 0;
my $conspecificSpThr = 0.1; #higher fraction of genes being two copies in the same sample (abundance >0), and the whole MGS is removed from that sample
my $multiCpyThr = 0.2; #kick out specific genes if too many copies / genome: should be relatively high (0.2+) as there are less intrusive mechanisms for removing such genes..
my $multCpyMGSthr = 0.7; #kick out entire MGS, if it has too many multicopy genes at above threshhold
my $MGStoolowGsThr = 10; #less genes than this in a single sample -> rm MGS from sample for strains
my $mode = "MGS";
my $appendWriteTrigger = 120; #every Xth samples, genes are written (to manage memory)
my $startSubFromMGS = ""; #debug option: only start resubmitting tree building from this MGS (e.g. "MGS.1382" )
#define local files..
my $lSNPdir="SNP"; my $lMAPdir = "mapping";
my $lConsFNA = "genes.shrtHD.SNPc.MPI.fna.gz";
my $lConsCTG = "contig.SNPc.MPI.fna.gz";
my $lConsFAA = "proteins.shrtHD.SNPc.MPI.faa.gz";
my $SNPcaller = "MPI";
my $lConsVCF = "allSNP.${SNPcaller}.vcf.gz";
my $lConsVCFsup = "allSNP.${SNPcaller}-sup.vcf.gz";


#set up some base paths specific to pipeline..
my $FNAstdof = "allFNAs.fna"; my $FAAstdof = "allFAAs.faa";
my $LINKstdof = "link2GC.txt"; my $CATstdof = "all.cat";
my $abundF="/assemblies/metag/ContigStats/Coverage.pergene.gz";
my $bamDepthFsuffix = "-smd.bam.coverage.gz";
my $bamDepthFsuffixSup = ".sup-smd.bam.coverage.gz";


checkMF();
#$treeFile = $ARGV[3] if (@ARGV > 3);$onlySubmit = $ARGV[4] if (@ARGV > 4);
#$doSubmit = $ARGV[5] if (@ARGV > 5);$subMode = $ARGV[6] if (@ARGV > 6);


GetOptions(
	"GCd=s"          => \$GCd,
	"outD=s"         => \$outDpre,
	"MGS=s"          => \$MGSfile,
	"nodeTmp|tmpD=s" => \$locTmpDir1, 
	"submit=i"       => \$doSubmit,
	"selfMemGb=i"    => \$selfMemGb,
	"onlySubmit=i"   => \$onlySubmit, #submit only jobs, or also recreate input fna/faa files? (can take days)
	"reSubmit=i"     => \$reSubmit, #for all MGS: resubmit tree phylo building
	#workflow HPC usage
	"subjob=i"       => \$subJob,
	"maxSubJob=i"    => \$maxSubJob,
	#"cores=i"        => \$numCores, #not used any longer..
	"maxCores=i"     => \$maxCores, #superseedes -cores, will dynamically allocate num cores based on input file size, if defined
	"presortGenes=i" => \$presortGenes, #how many potential genes to include, of the original MGS (receovered will vary strongly  between samples)
	"maxGenes=i"     => \$maxNGenes, #how many genes to try to include? -> will be decided on each samples
	"forceSNPcalls=i"  => \$forceVCF2FNA,
	"preCompConsSNP=i"   => \$preCompCons,
	
	#used genes fine tuning..
	"MGSminGenesPSmpl=i" => \$MGStoolowGsThr, #less genes than this in a single sample -> rm MGS from sample for strains. default 10
	"multiGeneSmplMax=f" => \$multiGeneSmplMax, #default 0.15
	"conspGeneSmplMax=f" => \$conspGeneSmplMax, #default 0.05

	"MGSphylo=s"     => \$treeFile,
	"submissionMode=s"      => \$subMode,
	"MSAprog=i"      => \$MSAprog,
	"MGset=s"        => \$useGTDBmg,
	"redoSubmissionData=i" => \$redoSubmissionData,  #for all MGS: will resubmit phylo and rebuild the fna/faa files..
	"deepRepair=i"   => \$deepRepair, #for missing MGS phylos: will resubmit phylo and rebuild fna/faa 
	"rmMSA=i"        => \$rmMSA, #remove MSA, to save diskspace
	"ContTests=s"      => \$contTests, #continous stat tests to be handed to next step (just a passthrough)
	"DiscTests=s"      => \$discTests, #discrete stat tests to be handed to next step (just a passthrough)
	"familyVar=s"      => \$familyVar, #column name in metadata containing family id
	"groupStabilityVars=s"      => \$groupStabilityVars, #column names of categories used for calculation of resilience and persistence
	"treeSubFromMGS=s" => \$startSubFromMGS, #debug option..
	
	#SNP calling
	"minSNPDepth=i"  => \$minSNPDepth,
	"minSNPCallQual=i"  => \$minSNPCallQual,

);



#die timeNice(20) ." ".timeNice(12252)."\n"; #TEST

#define global vars
my $QSBoptHR = emptyQsubOpt($doSubmit,"",$subMode);
my $MGSfileOri = $MGSfile; #save for later..

my $bindir;my $outD;my $scratchD;my $preConDir;my $LOGDIR;my $mapF;
my %map; my %AsGrps;my @samples;#map and assembly groups
my %ConspecificMGS; #list of conspecific MGS

my $gene2taxF; #where to find info what genes (gene cat)

prepRun();


my %AGlist; #list of assembly groups that need to be processed together;
createAGlist();
#foreach (sort keys %AGlist) {   print "$_ : @{$AGlist{$_}}\n";}die;

my %preCompSNPs;
preComputeConsSNP();


my %replN; #my %genesWrite; #keep stats/track

#my %allFNA; my %allFAA; #big hash with all genes in @allGenes
#my %gene2genes; #no longer needed
my %cl2gene2; #contains link from GCgene to fasta header assembly, cleaned up for multi copy already..
#my %SIcat;


my $SIgenes; my $Gene2COG; my $Gene2MGS; my $COGprios;
my %SIdirs; #unified storage of dirs per SI (SI==MGS)
my %MGSsmplConsp; #saves single samples within a MGS that seem to have too high rates of conspecific genes (controlled by $multiGeneSmplMax )





#key step to determine with set of genes (representing MGS) is to be MSA'd for strain phylos
#these might be very limited number of genes here..
($SIgenes,$Gene2COG,$Gene2MGS,$COGprios) = readGene2tax($gene2taxF,$presortGenes);#$maxNGenes);
#%SIgenes=%{$hr1};%Gene2COG=%{$hr2}; %Gene2MGS = %{$hr3}; %COGprios = %{$hr4};
my @specis = sort(keys(%{$SIgenes}));
#sort specis by numbers, so start with MGS1, MGS2 etc
my %sis; foreach (@specis){m/(\d+)$/; $sis{$_}=$1;}
@specis = sort {$sis{$a} <=> $sis{$b} } keys %sis;


#die "specis::\n@specis\n";
my $cnt=0; my $SaSe = "|"; my $dirsArePrepped = 1; my $allCatFileE = 1;

foreach my $SI (@specis){ #loop creates per specI file structure to run buildTreeScript on..
	#PART I: create fasta files required by tree
	my $outD2 = "$outD/$SI/";
	$SIdirs{$SI} = $outD2;
	#print "$outD2\n";
	$dirsArePrepped =0 unless (-d $outD2 &&  (-e "$outD2/data.log" || -e "$outD2/all.cat.tmp") ); # first phase only has "all.cat.tmp" file..
	$allCatFileE = 0 if (!-e  "$SIdirs{$SI}/$CATstdof");
	if (-d $outD2 && $onlySubmit == 0){#don't delete folders if we want to submit a job later..
		system "rm -rf $outD2/*";
	}
	system "mkdir -p $outD2" unless (-d $outD2);
	
}
$PhylosExist = 0 if (!$allCatFileE || !$dirsArePrepped);

my %smplsPerMGS; #stats: MGS is represented in how many different samples?

#hashes of strings that keep results to be written for each species..
my %OCstrH ; my %OFstrH ; my %OAstrH ; my %OLstrH ;


if ($dirsArePrepped == 0 || $onlySubmit == 0 
			|| $subJob){
	$PhylosExist=0;
	
	print "\n\n----------------------------------------------------\nPart I:: extracting relevant core MGS genes (SNP consensus called) from original assemblies\n----------------------------------------------------\n\n";
	
	filterMultiCopyGenes();
	$Gene2COG = {}; #delete, no longer needed..
	
	reportingsMGS();
	
	my @jobsMain;

	if ($maxSubJob && !$subJob){
		#here needs to submit itself maxSubJob times
		my $strain1scr = getProgPaths("MGS_strain1_scr"); #self reference
		#my $MGSfile1 = $MGSfile; $MGSfile1 =~ s/\.srt$//; #check that this isnt' the sorted MGS file being used in a subjob..
		my $selfCmd = "$strain1scr -GCd $GCd -outD $outD -MGS $MGSfileOri -submit $doSubmit -onlySubmit 0 -reSubmit 0  -maxSubJob $maxSubJob -MGSminGenesPSmpl $MGStoolowGsThr -multiGeneSmplMax $multiGeneSmplMax -conspGeneSmplMax $conspGeneSmplMax -MGSphylo $treeFile -presortGenes $presortGenes -maxGenes $maxNGenes -MGset $useGTDBmg -redoSubmissionData 0 -deepRepair 0 -rmMSA 0 -minSNPDepth $minSNPDepth -minSNPCallQual $minSNPCallQual -forceSNPcalls 0 -preCompConsSNP $preCompCons";
		$selfCmd .= "-tmpD $locTmpDir1" if ($locTmpDir1 ne "");
		
		my $tmpHDD=$QSBoptHR->{tmpSpace} ; $QSBoptHR->{tmpSpace} =15; #request some basic amount
		
		#submission of self-subjobs..
		for (my $sj = 1; $sj < $maxSubJob; $sj ++){
			my $cmdX = "$selfCmd -subjob $sj;\n";
			my $checkF = "$LOGDIR/mainExtr.${sj}.stone";
			$cmdX .= "touch $checkF\n";
			#die "$cmdX\n\n";
			print $LOGDIR."Strain1_B${sj}.sh\n";
			my ($dep,$qcmd) = qsubSystem($LOGDIR."Strain1_B${sj}.sh",$cmdX,1,"${selfMemGb}G","Str1.$sj","","",1,[],$QSBoptHR);
			push(@jobsMain,$dep);
		}
		$QSBoptHR->{tmpSpace} = $tmpHDD;
	}
	
	#and extract the corresponding fna/ faa from every other dir.. main single core work
	#this will also determine how many genes per MGS are now extracted..
	extractFNAFAA2genes();#@allGenes);
	%cl2gene2 = (); #no longer needed, delete
	#write logs to found genes etc.
	writeLogsStep1();
	
	if ($subJob){
		print "Finished subJob ${subJob}/$maxSubJob. Exiting..\n";
		exit(0);
	}

	if (@jobsMain && $maxSubJob && !$subJob){ # second part for main worker: check that everything else is finished..
		qsubSystemJobAlive( \@jobsMain,$QSBoptHR );
		#check if all required files present
		for (my $sj = 1; $sj < $maxSubJob; $sj ++){ #job 0 doesn't have stone (is this job..)
			my $checkF = "$LOGDIR/mainExtr.${sj}.stone";
			die "Can't find checkfile $checkF .. abort\n" unless (-e $checkF);
			unlink $checkF; #and delete..
		}
		
	}
	
	print "\nGene extraction & redistribution finished, ready to proceed to phylogeny jobs\n";

} else {
	print "Skipping Part I, outdir already prepared.\n";
}

#die;


#load some log files..
#if (scalar(keys(%genesWrite)) == 0) { #load genes found..
#	#read logs of found genes etc.
#	foreach my $SI (@specis){
#		my $outD2 = $SIdirs{$SI}; my $llogF="$outD2/geneFnd.log";
#		next unless (-e $llogF);
#		my $Lstr = `cat $llogF`; $Lstr =~ m/Total genes write (\S+): (\d+)/; 
#		$genesWrite{$1} = $2;
#		die "$llogF incorrect: $1 != $SI\n" if ($1 ne $SI);
#		$PhylosExist =0 if (!-d "$outD2/pjylo/");
#	}
#}
if (scalar(keys(%ConspecificMGS)) == 0){
	my $conlog = "$bindir/LOGandSUB/ConspecificMGS.log";
	open I,"<$conlog" or die "Can't open conspecific $conlog\n";
	while (my $l = <I>){my @spl = split /\t/,$l;$ConspecificMGS{$spl[0]} = [split(/,/,$spl[1])];}
	close I;
}




my %FNAref; my %FAAref;
my %SIgenes_OG; #later reads in SIgenes again, but no restriction to length 

my $geneCatLoaded=0;
#read in genecat to create outgroup fasta sequences..
if (1 && (!$allCatFileE || $deepRepair || !$dirsArePrepped || $onlySubmit == 0 || $redoSubmissionData == 1)){
	#also read reference gene seqs (for outgroup)
	my $refFNA = ""; my $refFAA = ""; my $refNameL = "unknw";
	if ($mode eq "MGS" || $mode eq "MGSall"){
		print "Reading reference genecat genes, to create outgroup sequences\n";
		$refFNA = "$GCd/compl.incompl.95.fna"; $refFAA = "$GCd/compl.incompl.95.prot.faa";
		$geneCatLoaded=1;$refNameL = "geneCat";
	} elsif ($mode eq "FMG"){
		print "reading FMG ref genes..";
		$refFNA = "$GCd/FMG/COG*.fna"; $refFAA = "$GCd/FMG/COG*.faa";
		$refNameL = "FMG ref";
	}
	
	my ($hr1,$Gene2COG_OG,$hr3,$hr4) = readGene2tax($gene2taxF);
	%SIgenes_OG = %{$hr1};
	#%SIgenes_OG=%{$hr1}; my %Gene2COG_OG=%{$hr2}; 
	$hr1 = readFasta($refFAA,1,"\\s",$Gene2COG_OG); %FAAref = %{$hr1};
	$hr1 = readFasta($refFNA,1,"\\s",$Gene2COG_OG); %FNAref = %{$hr1};
	print "read ". scalar(keys %FNAref)." genes from $refNameL\n";
	print "done\n";
}



print "\n\n----------------------------------------------------\nPart II:: resort .cat files, submit intraStrain phylogenies for " . scalar(@specis) . " MGS\n----------------------------------------------------\n\n";


die "Tree for outgroup specified, but file not found:$treeFile\nAborting..\n" if  ($treeFile ne "" && !-e $treeFile);

#die;
#go through every SpecI;
$cnt=0; my $lcnt=0; my @jobs;
foreach my $SI (@specis){ #loop creates per specI file structure to run buildTreeScript on..
	last if (!$reSubmit && !$redoSubmissionData && $PhylosExist);
	if (exists($ConspecificMGS{$SI}) && $ConspecificMGS{$SI}->[0] =~ m/multicopy/){print "Skipping $SI due to inclusion in conspecific MGS list.\n";next;}
	$lcnt++;
	if ($startSubFromMGS ne "" ){
		if ($SI ne $startSubFromMGS){next;
		} else { $startSubFromMGS = "";} #deactivate now
	}
	qsubSystemWaitMaxJobs($checkMaxNumJobs);
	#print "$SI  XX "; die;
	#print "$SI\n";
	#next unless ($lcnt>10);
	#PART I: create fasta files required by tree
	my $outD2 = $SIdirs{$SI};
	my $treeStone = "$outD2/treeDone.sto";
	#next if (-e "$outD2/phylo/IQtree_allsites.treefile");
	#print "$outD2\n";
	my $multiSmpl=0;
	system "mkdir -p $outD2" unless (-d $outD2);
	my $FNAtf = "$outD2/$FNAstdof"; my $FAAtf = "$outD2/$FAAstdof";
	my $CATtf = "$outD2/$CATstdof"; #my $Linkf = "$outD2/$LINKstdof";
	my $IQtreef= "$outD2/phylo/IQtree_allsites.treefile";
	
	#job done..
	print "${SI}::"; 
	#if (!exists($genesWrite{$SI}) ) { print "$SI does not exist in genesWrite object\n";
	#}elsif ($genesWrite{$SI} <10){print "WARNING: $SI has too few genes that could be found! Skipping..\n";next;}
	if (-e $treeStone && -e $IQtreef ){print "Skipping (tree exists?)..\n";next;}
	
	my $outgS = "";my $OG = "";
	if (-e "$outD2/data.log"){$OG = `cat $outD2/data.log`; chomp $OG; $OG=~s/^OG://;}
	my $contPhylo = 1; $contPhylo = 0 if ($reSubmit || $redoSubmissionData);
	
	#main command to build within species strain tree.. missing outgroup so far ($outgS)
	
	my $inputFNAsize = fileGZs($FNAtf) / (1024 * 1024); #size in MB
	if ($inputFNAsize ==0){print "empty input $FNAtf .. next.\n";next;} #empty input
	if (  ($MSAprog==4 && $inputFNAsize>700) ){ #only if FNA is > X mb   #$inputFNAsize  > 1200 ||
		$QSBoptHR->{useLongQueue} = 1 ;
		#deactivate for now, not needed really..
	}
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
	my $totMem = int($inputFNAsize /250*40)+10;$totMem = 20 if ($totMem < 20);
	my $numCoreL = $numCores;	
	if ($maxCores >0){ #scale cores according to used memory size
		$numCoreL = int($maxCores * (sqrt($totMem) / sqrt(180)));
		$numCoreL = 4 if ($numCoreL < 4);		$numCoreL = $maxCores if ($numCoreL > $maxCores);
	}
	
	my $subsSmpl = -1;my $useSuperTree = 0;
	my $bts = getProgPaths("buildTree_scr");
	my $Tcmd= "$bts -fna $FNAtf -aa $FAAtf -smplSep '\\$SaSe' -cats $CATtf -outD $outD2  -runIQtree 1 -runFastTree 0 -cores $numCoreL  "; 
	$Tcmd .= "-AAtree 0 -bootstrap 0 -NTfiltCount 400 -NTfilt 0.07 -NTfiltPerGene 0.5 -GenesPerSpecies 0.1 -runRaxMLng 0 -minOverlapMSA 2 ";
	$Tcmd .= "-subsetSmpls $subsSmpl -fracMaxGenes90pct 0.7 "; #concentrate on almost complete gene groups.. can yield more samples overall and speeds up calc..
	$Tcmd .= "-rmMSA $rmMSA -gzInput 1 "; #save more diskspace..
	$Tcmd .= "-SynTree 0 -NonSynTree 0 -MSAprogram $MSAprog -continue $contPhylo -AutoModel 0 -iqFast 1 -superTree $useSuperTree ";
	$Tcmd .= "-runDNDS 0 -runTheta 0 -tmpD $scratchD/$SI/ -map $mapF ";
	my $postCmd = "\n\ntouch $treeStone\n";
		#die "$cmd\n" if ($cnt ==10);
	
	
	#early submission.. no further work needed here!
	if ( (!$deepRepair || $OG ne "") && $redoSubmissionData == 0 && $onlySubmit==1 && !-e "$CATtf.tmp" && $dirsArePrepped  && fileGZe($FNAtf) && fileGZe($FAAtf) && fileGZe($CATtf)){
		unless (fileGZe($FNAtf) && fileGZe($FAAtf) && fileGZe($CATtf)){
			print "Can't find required input files:\n$FNAtf\n$FAAtf\n$CATtf\n" ;
			die if ($cnt <= 1);
			next;
		}
		$cnt ++;
		$outgS = " -outgroup $OG " ;
		#die "$totMem ;; $inputFNAsize\n\n";
		my ($dep,$qcmd) = qsubSystem($outD2."treeCmd.sh",$Tcmd.$outgS.$postCmd,$numCoreL,int($totMem/$numCoreL)+1 ."G","FT$cnt","","",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
		push (@jobs,$dep);
#		qsubSystem2($outD2."treeCmd.sh",$QSBoptHR,{cores => $numCores*2 .",".$numCores});
		#die "$outD2";
		next;
	}
	die "Gene cat wasn't loaded, check program logic.\n!$deepRepair && $redoSubmissionData == 0 && $onlySubmit==1 && $dirsArePrepped && !-e $CATtf.tmp \n" if (!$geneCatLoaded);
	my %SIcatLoc;
	if (-e "$CATtf.tmp"){
		open ICT,"<$CATtf.tmp" or die "Can't open cat file $CATtf.tmp\n";
		while (<ICT>){
			chomp; my @spl = split /\t/;
			#translation SIcat:
			##print OC "$SI\t$cog\t$sd3\t$ng\n";
			##$SIcat{$SI}{$cog}{$sd3} = $ng;
			die "$_\n$spl[0] eq $SI\n" unless ($spl[0] eq $SI);
			$SIcatLoc {$spl[1]} {$spl[2]} = $spl[3];
		}
		close ICT;
	} elsif (-e $CATtf) {
		#print OC $SIcatLoc{$cog}{$smpl};				print OC "\t".$SIcatLoc{$cog}{$smpl};		my $ng = "$OG$SaSe$cog";
		#print "Reconstructing tmp cat file.";
		open ICT,"<$CATtf" or die "Can't open (precompiled) cat file $CATtf\nConsider deleting strain dir and rerunning strainMGS script\n";
		my $catLines=0;my $cntItems=0;
		while (<ICT>){
			chomp; my @spl = split /\t/;
			foreach my $tags (@spl){
				#my @spl2 = split (/\\$SaSe/,$tags);
				#$SIcatLoc {$spl2[1]}{$spl2[0]}  = $tags;print "$spl2[1] : $spl2[0]  = $tags\n";
				$tags =~ m/^(.*)\|(.*)$/;		
				$SIcatLoc {$2}{$1}  = $tags;				
				$cntItems++;
				#print "$2 $1  = $tags\n";
			}
			$catLines++;
		}
		close ICT;
		print "${SI}:: $catLines cat lines, $cntItems items: $CATtf\n";
	} else {
		print "WARNING:: ${SI}:: possible error: neither .cat nor .cat.tmp exists in $outD2\n";
		next;

	}
	
	
	#my @curCogs = sort keys %{$SIcat{$SI}};
	my @curCogs = sort keys %SIcatLoc;
	if (scalar(@curCogs) < 10){
		die "$SI error: \@curCogs is empty or < 10 members\n$CATtf\n";
	}
	#print "COGs: $curCogs[0] $curCogs[1]\n";
	
	#include outgroup?
	if ($treeFile ne ""){
		my $neiTree = getProgPaths("neighborTree");
		my $call = "$neiTree $treeFile $SI";
		#print "$call\n";
		my $OG1 = `$call`; chomp $OG1;
		die "Can't find outgroup from call $call\n\n" if (!defined $OG1);
		my @sspl = split /\s/,$OG1; $OG = "";
		next if (@sspl == 0);
		$OG = $sspl[0] if (@sspl); my $cntX=-1;
		my $cntShrCogs=0;
		while ( defined($OG) && $cntShrCogs < 10 && $cntX < @sspl){
			$cntX++; 
			$OG=$sspl[$cntX];
			if (!exists($SIgenes_OG{$OG})){
				next;
			}
			#$cntShrCogs=0;
			if (exists($SIgenes_OG{$OG})){
				foreach my $cog (@curCogs){
					next unless (exists($SIgenes_OG{$OG}{$cog}));
					next unless (exists($FNAref{$SIgenes_OG{$OG}{$cog}}));
					next if ($cog =~ m/uniq\d+/);
					$cntShrCogs ++;
				}
			}
			#last if ($cntX >= @sspl || $cntShrCogs >= 10);
		}
		if ($cntShrCogs < 1 && $cntX>0){
			print "Could not find outgroup for $SI!!\n@sspl\n@curCogs[1..10]\n";
		}
		unless (exists($SIgenes->{$OG})){
			print "can't find speci $OG\n$OG1\n";
			$OG="";
		} else {
			$outgS = "-outgroup $OG " ;
		}
		print "outgroup $OG ";
		#next;
	}
	
	#and fasta/faa/cat files..
	#open OL,">$Linkf" or die "Can't open link file $Linkf\n";
	#append to FNA/FAA (for outgroups)
	
	my %uniqSmpls;my $OGgenesUsed=0;
	my $tmpFAAog = ""; my $tmpFNAog = "";
	foreach my $cog (@curCogs){
		if ($OG ne "" && exists($SIgenes_OG{$OG}{$cog})){#deal with outgroup..
			next unless (exists($SIgenes_OG{$OG}{$cog}) && exists($FNAref{$SIgenes_OG{$OG}{$cog}}));
			next if ($cog =~ m/uniq\d+/);
			my $ng = "$OG$SaSe$cog";
			$tmpFNAog .= ">$ng\n$FNAref{$SIgenes_OG{$OG}{$cog}}\n";
			$tmpFAAog .= ">$ng\n$FAAref{$SIgenes_OG{$OG}{$cog}}\n";
			#$SIcat{$SI}{$cog}{$OG} = $ng;
			$SIcatLoc{$cog}{$OG} = $ng;
			$OGgenesUsed++;
			#if ($new){ print OC "$ng";$new=0;
			#} else { print OC "\t$ng"; }
		}
	}
	open OF,">>$FNAtf" or die "Can't open NT file $FNAtf\n";
	open OA,">>$FAAtf" or die "Can't open AA file $FAAtf\n";
	print OF "$tmpFNAog"; print OA "$tmpFAAog";
	close OA; close OF;
	
	#print "used $OGgenesUsed genes  ";
		
	open OC,">$CATtf" or die "Can't open cat file $CATtf\n";
	foreach my $cog (@curCogs){
		my $cntL=0;
		foreach my $smpl (sort keys %{$SIcatLoc{$cog}}){
			if ($cntL==0){
				#print OC $SIcat{$SI}{$cog}{$smpl};
				print OC $SIcatLoc{$cog}{$smpl};
			} else {
				print OC "\t".$SIcatLoc{$cog}{$smpl};
			}
			$cntL++;
			$uniqSmpls{$smpl} = 1;
		}
		print OC "\n";
	}
	close OC;
	print "Generated CAT file ";
	if ($OGgenesUsed ==0 && $OG ne ""){
		die "Couldn't include any outgroup genes! $OG\n$FNAtf\n";
	}
	
	#note done somewhere how many genes these actually are..
	system "echo \"OG:$OG\" > $outD2/data.log";
	$multiSmpl = scalar(keys %uniqSmpls);
	if ($multiSmpl>2){
		print "$SI: multiSmpls:\t$multiSmpl\tpotential genes: ". scalar(@curCogs) ."\tcores:$numCoreL mem:$totMem\n";
	} else {print "\n$SI: too few samples ($multiSmpl) for tree stats\n";next;}
	
	system "rm -f $CATtf.tmp\n";
	#PART II: qsub tree build command
	
	#die "$cmd\n" if ($cnt ==10);
	my ($dep,$qcmd) = qsubSystem($outD2."treeCmd.sh",$Tcmd.$outgS.$postCmd,$numCoreL,int($totMem/$numCoreL)+1 ."G","FT$cnt","","",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	$QSBoptHR->{useLongQueue} = 0;
	$cnt ++;
	push (@jobs,$dep);
	#die $outD2."treeCmd.sh\n";

}
#too many jobs to use as job dependency..
qsubSystemJobAlive( \@jobs,$QSBoptHR );
print "\nAll done for $cnt Bins\nRun strain_within_2.pl for summary stats:\n";

my $outDX =  $MGSfile;#"$GCd/$mode/intra_phylo/";
$outDX =~ s/[^\/]+$//;
my $MGSabundance = "$GCd/Anno/Tax/GTDBmg_MGS/specI.mat";
$MGSabundance = "$bindir/Annotation/Abundance/MGS.matL7.txt";

my $strain2Scr = getProgPaths("MGS_strain2_scr");

my $nxtCmd = "$strain2Scr -GCd $GCd -FMGdir $outD -MGSmatrix $MGSabundance -map $mapF -cores 4 -Hcores $maxCores -reSubmit 0 -DiscTests \"$discTests\" -ContTests \"$contTests\" -familyVar \"$familyVar\" -groupStabilityVars \"$groupStabilityVars\" \n"; #$GCd/MB2.clusters.ext.can.Rhcl.matL0.txt
	my ($dep,$qcmd) = qsubSystem($outD."/strainAnalysis2.sh",$nxtCmd,1,"60G","2StrainSub","","",1,[],$QSBoptHR);
print "\n". $nxtCmd."\n";


#cleanup
system "rm -rf $locTmpDir";
system "rm -rf $preConDir" if ($preCompCons);


exit(0);

 

#########################################################################################
#########################################################################################



sub writeLogsStep1{


	#print log file
	my $conlog = "$bindir/LOGandSUB/ConspecificMGS.log";
	open LO,">>$conlog" or die "Can't open conspecific log file: $conlog\n";
	foreach my $MGS (keys %ConspecificMGS){
		print LO $MGS . "\t" . join(",",@{$ConspecificMGS{$MGS}}) . "\n";
	}
	close LO;
}


sub filterMultiCopyGenes{
	
	print "Preparing base strain alignments, per MGS\nThis might take a good while..\n";
	my ($hr1,$cl2gene) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",0,$Gene2COG);
	$hr1 = {}; #my %cl2gene = %{$hr2}; 
	#stores alt names (wihtout M4__ at end)
	#read binning based on SpecI's 

	#my @allGenes;
	my $goodGene =0; my $badGene =0; my $geneCntd=0;
	my %multiGnInMGS; #stores genes kicked out in MGS due to multi copy filtering done in this routine
	my %singleGnInMGS;
	my %multiGeneMGSsmpl;
	my %MGSgeneCnt;
	my $NgenesTotal = scalar(keys %{$Gene2COG});
	my $NMGSTotal = scalar(keys %{$COGprios});
	foreach my $MGS(keys %{$COGprios}){
		$multiGnInMGS{$MGS} =0;$singleGnInMGS{$MGS} =0;
	}
	#1: in this process we can also check for multi genes (genes represent by >1 gene in an assembly)
	foreach my $gene (keys %{$Gene2COG}){
		my $geneStr = $cl2gene->{$gene}; $geneStr =~ s/>//g;
		my %tmpGen; 
		#my @genegenes = split /,/,$geneStr;die "$gene - @genegenes\n$geneStr\n";#2828988 - SMPL6M31__C1607771_L=26964=_19 SMPL2M32__C332244_L=10947=_18
		#finding the sample might be slow in the start.. but will pay off long term
		foreach my $sg (split /,/,$geneStr){
			my @spls = split /__/,$sg;#die "{$gene}{$spls[0]}\n";
			push(@{$tmpGen{$spls[0]}}, $sg); #gene in several samples..
		}
		my $MGSn = $Gene2MGS->{$gene};
		$MGSgeneCnt{$MGSn}++;
		foreach my $sm(keys %tmpGen){ #go over samples found..
			if (scalar(@{$tmpGen{$sm}}) > 1){
				$multiGeneMGSsmpl{$MGSn}{$sm}++;
				#print "@{$tmpGen{$sm}}\n";
			}
		}
	}
	
	
	#2: catalogue if an MGS should not be used in a certain samples..
	my $MGSwithConspSmpl=0; my $MGSwithoutConspSmpl=0; 
	my $totalConsSmpls=0; my $totalNONConsSmpls=0;
	foreach my $MGS(keys %{$COGprios}){
		my $totGenes = $MGSgeneCnt{$MGS} ;
		#print "$MGS - $totGenes\n";
		if ($MGSgeneCnt{$MGS} == 0 || !exists($multiGeneMGSsmpl{$MGS})){$MGSwithoutConspSmpl++;next;}
		my $locConspec =0; my $locNonCons=0;
		foreach my $sm (keys %{$multiGeneMGSsmpl{$MGS}}){
			my $lfrac = ($multiGeneMGSsmpl{$MGS}{$sm} / $totGenes);
			if ($lfrac > $multiGeneSmplMax){
				$locConspec++;
				$MGSsmplConsp{$MGS}{$sm}=1;
				
				#print "MGSsmplConsp{$MGS}{$sm};\n";#MGSsmplConsp{MGS.66}{IL140M24};
			} else {$locNonCons++;}
		}
		if ($locConspec>0){$MGSwithConspSmpl++ ;
		} else {$MGSwithoutConspSmpl++;
		}
		$totalNONConsSmpls += $locNonCons ;$totalConsSmpls+=$locConspec;
	}
print $MGSwithConspSmpl ."/" . ($MGSwithConspSmpl + $MGSwithoutConspSmpl) . " MGS with >1 multi-copy gene sample. " . $totalConsSmpls . "/" . ($totalConsSmpls+$totalNONConsSmpls) . " samples with MGS >$multiGeneSmplMax multi-copy gene threshold across all-sample-MGS combinations.\nTotal MGS: ${NMGSTotal}, total genes: $NgenesTotal\n";

		
	#3: check if an MGS is constantly represented by multi copies, taking into account that difficult samples were removed in (2)
	foreach my $gene (keys %{$Gene2COG}){
		my $geneStr = $cl2gene->{$gene}; $geneStr =~ s/>//g;
		my %tmpGen; 
		foreach my $sg (split /,/,$geneStr){
			my @spls = split /__/,$sg;push(@{$tmpGen{$spls[0]}}, $sg); 
		}
		my $MGS = $Gene2MGS->{$gene};
		my $multiGeneCnt =0;
		foreach my $sm(keys %tmpGen){ #go over samples found..
			if (scalar(@{$tmpGen{$sm}}) > 1 && !exists($MGSsmplConsp{$MGS}{$sm}) ){
				$multiGeneCnt++ ;
			}
		}
		

		#fraction of samples where gene was double. Genes were this is obviously too high need to be excluded (multi copy)
		#my $frac = $multiGeneCnt/$totSmpl;
		my $frac = $multiGeneCnt/keys(%tmpGen); #print "$frac = $multiGeneCnt/$totSmpl \n"; #fraction of samples with >1 gene   #@genegenes;
		if ($frac > $multiCpyThr){
			$badGene++;$multiGnInMGS{$MGS}++;
			next;
		}
		$goodGene ++;$singleGnInMGS{$MGS}++; #in principle good gene (single copy mostly).. but could still be multicopy in single samples..
		foreach my $sm(keys %tmpGen){
			#next if (scalar(@{$tmpGen{$sm}}) > 1); #don't copy this gene in this sample .. no decission which is the right one.. -> too harsh?? check for each sample within assembly group instead..
			#here gene is saved for later usage..
			@{$cl2gene2{$sm}{$gene}} = @{$tmpGen{$sm}};
			$geneCntd++;
			$smplsPerMGS{$MGS}{$sm}++;
		}
	}
	$cl2gene = {}; #lessen mem
	print "Found $goodGene / ".($badGene + $goodGene)." genes useable (<$multiCpyThr, \"single copy\"). $geneCntd gene-sample combis added.\n";
	
	#find out if any MGS had a too high rate of multicopy genes
	my %multiCpyRateMGS;
	foreach my $MGS(keys %{$COGprios}){
		my $totGenes = $multiGnInMGS{$MGS} + $singleGnInMGS{$MGS};
		if ($totGenes>0){
			$multiCpyRateMGS{$MGS} = $multiGnInMGS{$MGS} / ($totGenes );
		} else {
			$multiCpyRateMGS{$MGS} = 0;
		}
	}
	my @MGSsrt = sort {$multiCpyRateMGS{$b} <=> $multiCpyRateMGS{$a} } keys %multiCpyRateMGS;
	print "Too high multicopy genes in MGS (>$multCpyMGSthr): "; my $lcnt=0; 
	foreach my $MGS (@MGSsrt){
		next if ($multiCpyRateMGS{$MGS} < $multCpyMGSthr ); 
		print "$MGS $multiCpyRateMGS{$MGS}; "; $lcnt++; 
		#push(@{$ConspecificMGS{$MGS}},"multicopy:$multiCpyRateMGS{$MGS}"); -> way too strict!
	} 
	print "(none)" if (!@MGSsrt);
	print "\n";
}

sub prepRun{

	$mode = "FMG" if ($MGSfile eq "");
	if ($mode eq "FMG"){$takeAll = 0;}
	if ($takeAll){$maxNGenes = -1;$mode="MGSall"; $doSubmit = 0;}


	$bindir = $MGSfile;$bindir =~ s/[^\/]+$//; 
	$outD =  $bindir."/intra_phylo/";#"$GCd/$mode/intra_phylo/";
	if ($outDpre ne ""){
		$outD = $outDpre ; 
		$outD .= "/" unless ($outD =~ m/\/$/);
		}
	$LOGDIR = "$outD/LOGandSUB/";
	$SNPconsLOGs = "$outD/SNPconsCalls.log" if ($SNPconsLOGs eq "");

	$GCd =~ m/.*\/([^\/]+)\/*/; my $GCname = $1;
	$outD =~ m/.*\/([^\/]+)\/*/; my $outDname = $1;
	$scratchD = getProgPaths("globalTmpDir",0);
	$scratchD .= "/strainsScr1/$GCname.$outDname/";
	#die "$scratchD  :$GCname :$GCd\n";
	if ($locTmpDir1 eq ""){
		my $locTmpN = getProgPaths("nodeTmpDir",0) ;
		my $suffix = ""; $suffix = "/SJ.${subJob}/" if ($subJob);
		if ($locTmpN eq ""){
			$locTmpDir =  "$outD/strainsScr1/$GCname.$outDname/$suffix" ; 
		} else {
			#my $tmp = `echo \$SLURM_LOCAL_SCRATCH`;
			#print STDERR "echo $locTmpN\n$tmp\n";
			#$locTmpN =~ s/\$/\\\$/;$locTmpN = `echo $locTmpN;`; #eval in sys
			$locTmpN=truePath($locTmpN,1);
			#die $locTmpN."\n";
			$locTmpDir =  "$locTmpN/strainsScr1/$GCname.$outDname/$suffix" ; 
		}
	} else { $locTmpDir = $locTmpDir1;}
	
	$preConDir = "$scratchD/preComp/";
	system "mkdir -p $locTmpDir $scratchD";


	print "\n!! WARNING !!: RESUBMISSION mode selected (will resubmit MSA + phylos even for already completed MGS) !!\n" if ($reSubmit);
	print "\n!! WARNING !!: REDOSUBMISSIONDATA mode selected (will redo and resubmit MSA + phylos even for already completed MGS) !!\n" if ($redoSubmissionData);

	$mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;chomp $mapF;
	
	#read info gene <-> taxonomy from this file, depends on config..
	$gene2taxF = "$GCd/FMG/gene2specI.txt";
	$gene2taxF = "$GCd/GTDBmg/gene2specI.txt" if ($useGTDBmg eq "GTDB");
	if ($mode eq "MGS" || $mode eq "MGSall"){
		$gene2taxF = createGene2MGS($MGSfile,$GCd);
		print "Using MGS from $MGSfile, adding eggNOG in: $gene2taxF\n";
	} 
	#$mapF = $GCd."LOGandSUB/inmap.txt" if ($mapF eq "");
	my ($hr1,$hr2) = readMapS($mapF,-1);
	%map = %{$hr1}; %AsGrps = %{$hr2};
	#get all samples in assembly group, but only last in mapgroup
	@samples = @{$map{opt}{smpl_order}};


	


	#die;

	#---------------
	#everything after is only for main submission job..
	if ($subJob){
		print "=============\n=============\nStrain_within v$version, subjob ${subJob}/$maxSubJob\n=============\n=============\n";
		return;
	} else {
			print "============= Strain_within v$version =============\n";
		print "Creating within species strains for ${mode}s in $GCd\n";
		print "GC dir: $GCd\nIn Cluster: $MGSfile\nCores: $numCores (max: ${maxCores})\n";
		print "MAP: $mapF\n";
		#print "Ref tree: $treeFile\n";
		print "Using tree $treeFile to create automatically outgroups\n" if ($treeFile ne "");
		print "Outdir: $outD\nTmpDir: $locTmpDir\n";
		print "MGs: $useGTDBmg\nGene2Tax: $gene2taxF\n";
		print "Using $presortGenes genes from each MGS for location\n";
		print "Deep repariing remaining submission files\n" if ($deepRepair);
		print "Pre-creating ConsSNPs in $preConDir in $preCompCons runs\n" if ($preCompCons);
		print "minSNPDepth=$minSNPDepth, minSNPCallQual=$minSNPCallQual\n";
		print "DiscTests=$discTests\n" unless ($discTests eq "");
		print "ContTests=$contTests\n" unless ($contTests eq "");
		print "familyVar=$familyVar\n" unless ($familyVar eq "");
		print "groupStabilityVars=$groupStabilityVars\n" unless ($groupStabilityVars eq "");
		print "MSAaligner: $MSAprog\n";
		
		
		if ($takeAll){print "**************** Take all genes MGS mode\n";}
		else {print "Using first $maxNGenes genes found per sample\n";}
		print "==============================================\n";
		if ($onlySubmit){print "Only submission mode\n";
		} elsif (!$subJob) {
			print "Creation of strain genes, old data might be deleted!\nDo you want to continue? (10s wait, use Ctrl-c to abort)\n"; sleep 10;
		}
	}
	#prep sorted MGS gene file
	if (!-e $MGSfile.".srt" && !$subJob){
		die "In rewriting loop while in a subjob.. aborting\n" if ($subJob);
		print "base files missing.. preparing complete resubmission and recalc of data\n";
		system "rm -fr $outD";
		system "rm -f $MGSfile.srt*";
		my $sortMGSgenes = getProgPaths("sortMGSGeneImport_scr");
		my $cmd = "$sortMGSgenes $GCd $MGSfile $useGTDBmg\n";
		print "$cmd\n";
		systemW $cmd;
	} else {
		print "Continuing on preped .srt files\n";
	}
	$MGSfile .= ".srt";
	print "\nnew MGS file: $MGSfile\n\n";

	system "mkdir -p $outD" unless (-d $outD);
	system "mkdir -p $LOGDIR" unless (-d $LOGDIR);
	
	#DEBUG
	system "rm -rf $preConDir;mkdir -p $preConDir" if ($preCompCons && !$subJob);

	#STONES
	system "mkdir -p $outD/stones/" unless (-d "$outD/stones/");
	
	if (system "mkdir -p $locTmpDir" ){ die "mkdir tmp failed"; }
	system("touch $locTmpDir/test.txt");
	if ( ! -e "$locTmpDir/test.txt"){die "Couldn't create test file in local dir $locTmpDir\n";}
#die "passed $locTmpDir\n";


	return;
}

sub preComputeConsSNP{
	my $inputChkd = 0;
	my $inputChk = "$outD/stones/0.fileChk.sto";
	$inputChkd =1 if (-e "$inputChk");
	my $fileAbsent = 0; 
	my $submPreComp = 1;#DEBUG
	$submPreComp = 0 if ($subJob);

	
	my @accumVCFcmds; my $BatchCnt=0;my @jobsPre;
	foreach my $smpl (@samples){ # just check that files are there..
		#check if SNP file is present
		last if ($onlySubmit && $inputChkd && !$preCompCons);
		my $cD = $map{$smpl}{wrdir}."/";
		#my $tarF = $cD."/SNP/genes.shrtHD.SNPc.MPI.fna.gz";
		my $tarF = $cD."/$lSNPdir/$lConsFNA";
		my $tarF2 = $cD."/$lSNPdir/$lConsFAA";
		my $tarVCF = $cD."/$lSNPdir/$lConsVCF";
		if (! fileGZe($tarVCF) || (-e $tarF && -e $tarF2 ) ){
			print "Can't find SNP files: $cD\n" unless (-e "$cD/SMPL.empty");
			$fileAbsent = 1;
			#die "$tarVCF\n";
		}
		if ($forceVCF2FNA && ! fileGZe($tarVCF)){
			die "Option -forceSNPcalls 1 used, but missing vcf $tarVCF\n";
		}
		
		if ($preCompCons && ( $forceVCF2FNA  || (! fileGZe( $tarF ) && fileGZe($tarVCF)))){
			#store these in scratch, uncompressed (much faster)
			my $fastaf = "$preConDir/$smpl.cons.genes.fna.gz";
			my $fastafAA = "$preConDir/$smpl.cons.prots.faa.gz";
			my $vcf2fnaCmd = createConsFastas($cD, $smpl, $fastaf, $fastafAA, 0, 1);
			$preCompSNPs{$smpl}{NT}=$fastaf;$preCompSNPs{$smpl}{AA}=$fastafAA;

			push(@accumVCFcmds,$vcf2fnaCmd);
			
			if (@accumVCFcmds > $preCompCons){
				print "submitting precomp batch $BatchCnt\n" if ($submPreComp);
				my $cmdX = "\necho \"BATCH $BatchCnt\"\nmkdir -p $preConDir/;\n\n" . join("\n",@accumVCFcmds);
				my $tmpSHDD=$QSBoptHR->{tmpSpace} ; $QSBoptHR->{tmpSpace} =0;
				my ($dep,$qcmd) = qsubSystem($LOGDIR."PreCompConsSNP_B${BatchCnt}.sh",$cmdX,1,"10G","ConsSNP$BatchCnt","","",$submPreComp,[],$QSBoptHR);
				$QSBoptHR->{tmpSpace} =$tmpSHDD;

				push(@jobsPre,$dep);
				#reset counters etc
				$BatchCnt++; @accumVCFcmds=();

				#die;
			}
		}
	}
	#last batch of jobs..
	if (@accumVCFcmds){
		
		my $cmdX = "\necho \"BATCH $BatchCnt\"\nmkdir -p $preConDir/;\n\n" . join("\n",@accumVCFcmds);
		my ($dep,$qcmd) = qsubSystem($LOGDIR."PreCompConsSNP_B${BatchCnt}.sh",$cmdX,1,"10G","ConsSNP$BatchCnt","","",$submPreComp,[],$QSBoptHR);
		push(@jobsPre,$dep);

	}
	if (@jobsPre){
		qsubSystemJobAlive( \@jobsPre,$QSBoptHR );
	}
	unless (-e "$inputChk"){
		if ($fileAbsent){
			print "Not all required input present\n" ;
		} else {
			print "All samples have SNP calls\n";
		}
		system "touch $inputChk" ;
	}
}


sub createAGlist{
	foreach my $smpl (@samples){ #fill up AGlist
	#and fill %AGlist .. so always let run..
		next if ($map{$smpl}{AssGroup} eq "-1");
		my $cAssGrp = $map{$smpl}{AssGroup};
		my $cMapGrp = $map{$smpl}{MapGroup};
		#print "$smpl $cAssGrp $cMapGrp\n";
		$AsGrps{$cMapGrp}{CntMap} ++;
		if (!exists($AsGrps{$cMapGrp}{CntMap})){ die "Can;t find CntMap for $cMapGrp";}
		next if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} );
		push(@{$AGlist{$cAssGrp}} , $smpl);
		if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} ){
			next;
		}
	}
}

sub histoMGS{#specifically for MGS..
	my ($aref,$msg) = @_;
	my @cnts = @{$aref};
	my @binSiz = (10,20,30,50,70,100,200,300,500,700,1000,2000,5000,10000,1e6);
	my %binC; #my $prevC=0;
	foreach (@binSiz){$binC{$_} = 0;}
	foreach my $c(@cnts){
		#print "$c ";
		my $bs=0;
		while ($c > $binSiz[$bs] && $bs < @binSiz){
			$bs++;
		}
		#print " X$c:${bs}X ";
		$binC{$binSiz[$bs]} ++;
	}
	#display bin counts..
	print $msg.": ";#"Bin size distribution: ";
	foreach (@binSiz){print " <$_:$binC{$_} " if ($binC{$_}>0);}
	print "\n";
	#DEBUG
	#print @cnts." : @cnts\n";
}



sub appendWriteMGSgenes{
	#write genes to respective MGS intra phyla..
	my ($writeLink) = @_;
	my $wrMGS=0;
	my @smplSet = keys(%OFstrH);
	my @smplSetS = shuffle(@smplSet); #shuffle to further reduce chance of multiple jobs writing consistently to the same files..

	foreach my $SI (@smplSetS){
		next if ($OFstrH{$SI} eq "");
		my $hasSlept=0;
		#handle file paths..
		my $outD2 = $SIdirs{$SI};
		my $FNAtf = "$outD2/$FNAstdof"; my $FAAtf = "$outD2/$FAAstdof";my $Linkf = "$outD2/$LINKstdof";
		my $CATtf = "$outD2/$CATstdof.tmp";
		my $blockF = "$outD2/block.tmp";
		
		#block sample for other writes..
		while (-e $blockF){sleep(5);$hasSlept=1;}
		while ($hasSlept && -e $blockF){sleep(8);}#second security layer..
		system "touch $blockF";
		sleep(2) if ($hasSlept);#security that other process has finished writes completely

		#writing strings out..
		open OF,">>$FNAtf" or die "Can't append NT file $FNAtf\n";print OF $OFstrH{$SI}; close OF;
		open OA,">>$FAAtf" or die "Can't append AA file $FAAtf\n";print OA $OAstrH{$SI}; close OA;
		if ($writeLink){open OL,">>$Linkf" or die "Can't append link file $Linkf\n" ; print OL $OLstrH{$SI}; close OL;}
		#this is only a temp file, that needs to be rewritten later..
		open OC,">>$CATtf" or die "Can't append to CAT file $CATtf\n";print OC $OCstrH{$SI} ; close OC;
		
		
		$OCstrH{$SI} = ""; $OFstrH{$SI} = ""; $OAstrH{$SI} = ""; $OLstrH{$SI} = "";
		$wrMGS++;
		
		system "rm -f $blockF";
	}
	print "\nwrote for $wrMGS MGS data..\n";
}


sub reportingsMGS{
	#eval #sample/MGS
	my %smplPmgs;
	foreach my $MGS (keys %smplsPerMGS){
		foreach my $sm (keys %{$smplsPerMGS{$MGS}}){
			$smplPmgs{$MGS}++ if ($smplsPerMGS{$MGS}{$sm} > 10);
		}
	}
	my $qt50=quantileArray(0.5,values(%smplPmgs));my $qt90=quantileArray(0.90,values(%smplPmgs));
	my @smplNs = values(%smplPmgs);@smplNs = sort { $b <=> $a}  @smplNs;
	print "Samples/MGS: QTL 50,90: $qt50 $qt90 . Top 5: $smplNs[0] $smplNs[1] $smplNs[2] $smplNs[3] $smplNs[5]\n";
	#die;
	
}

sub timeNice($){
	my ($tIN) = @_;
	$tIN = int($tIN);
	if ($tIN > (3600)){
		my $remMin = ($tIN%3600);
		return int($tIN/3600)."h".int($remMin/60)."m" . ($remMin%60) . "s";
	}
	if ($tIN > 60){
		return int($tIN/60)."m" . ($tIN%60) . "s";
	}
	return $tIN . "s";
}




#this routine hast to get genes out of each sample, that are needed
#and save them to be later written per specI
sub extractFNAFAA2genes{
	my %perMGScnts;
	my $gnCnt=0;
	#my %totGnes;
	#create gene to genes list
	foreach my $sm (keys %cl2gene2){
		#my @locGenes;
		#print "$sm ";
		my $MGSgeneCnt=0;
		foreach my $gn (keys %{$cl2gene2{$sm}}){
			#$totGnes{$gn} = 1;
			$gnCnt++;
			if (exists($Gene2MGS->{$gn})){
				$perMGScnts{$Gene2MGS->{$gn}}{$gn}=1;
				#print "1";
				$MGSgeneCnt++;
			}
		}
		#print "$sm  $gnCnt $MGSgeneCnt \n";
	}
	my @histoMGScnts ;#= values %perMGScnts;
	foreach my $MGS (keys %perMGScnts){
		my $perMGSgenes = scalar( keys( %{$perMGScnts{$MGS}} ) );
		push(@histoMGScnts,  $perMGSgenes);
		if ($perMGSgenes < 10){print "WARNING: only $perMGSgenes genes/COGs for MGS $MGS - MGS genes might be multi copy?\n";}
	}
	#DBUG
	print "Genes per MGS (prefiltering, N= ". $gnCnt  ." genes, " .scalar(keys(%perMGScnts)) . " MGS, avg " . int(0.5+ $gnCnt /scalar(keys(%perMGScnts))) . " genes/MGS):\n";	
	histoMGS(\@histoMGScnts,"Theorectical best Bin sizes: ");
	#some stats on genes/MGS
	my @srtdSmpls = sort (keys %cl2gene2);
	
	#subjob? split up what samples to process..
	if ($maxSubJob){
		my $Ndirs = scalar(@srtdSmpls);
		my $Nsmpls=0;
		foreach my $sd(keys %AGlist){
			$Nsmpls += scalar (@{$AGlist{$sd}}); #@{$AGlist{$cAssGrp}}
		}
		print "total samples: $Nsmpls , ASgrps: $Ndirs\n";
		#my $start = int( (($subJob-1)/$maxSubJob) * $Ndirs ) ;
		#my $end = int( (($subJob)/$maxSubJob) * $Ndirs ) ;
#		for (my $i=$start;$i < $end; $i+=$maxSubJob){
		my @smp2;my %smp2H;
		for (my $i=($subJob);$i < $Ndirs; $i+=$maxSubJob){
			push(@smp2,$srtdSmpls[$i]); 
			$smp2H{$srtdSmpls[$i]} =1;
		}
		@srtdSmpls = @smp2;
		print "\nSUBJOP: choosen samples: @smp2\n\n";
		#clean up hashes..
		foreach my $sm (keys %cl2gene2){
			delete($cl2gene2{$sm}) unless (exists ($smp2H{$sm}));
		}

	}
	
	
	
	print "Extracting GC genes from " . scalar(@srtdSmpls). " dirs\n";

	
	#different way to go over genes..
	 my $smCnt=1;
	 #storage hash for raw fasta/faa/link files, needs to be written separately
	#goes over every assembly group to extract SNP corrected genes that fall into each MGS
	my $writeLink = 1; my $appCnt=0;
		#DEBUG	@srtdSmpls = ("PDB3.F");
	
	
	my $sttime = time;	
	foreach my $sm (@srtdSmpls){

		print "\nAT SMPL:: $smCnt/" . scalar(@srtdSmpls) ." $sm - ". "Elapsed time : ", timeNice(time - $sttime) . "\n";
		#readGenesSample_Singl($sm, $OFstrHR, $OAstrHR, $OCstrHR, $OLstrHR, $writeLink,$sttime);
		readGenesSample_Singl($sm, $writeLink,$sttime);
		$smCnt++; $appCnt++;
		
		if ($appCnt >= $appendWriteTrigger){
			appendWriteMGSgenes( $writeLink);
			$appCnt=0;
		}
	}
	
	
	appendWriteMGSgenes($writeLink);
	print "Done writing all genes to subdirs, elapsed time: " . timeNice(time - $sttime)  . "\n";
	$appCnt=0;
	#done at the point with gene extractions
	return;
}

sub createConsFastas{
	my ($cD,$sm, $oFNA, $oFAA,$append2LOG,$returnCmd) = @_;
	my $vcf2fnaBin = getProgPaths("vcf2fna");
	my $vcf2fnaOpt = "";
	#my $seqPlatf = "hiSeq"; #-> get this from .map ..
	my $refFA = getAssemblContigs($cD); my $refGFF = getAssemblGFF($cD);
	my $depthFile = "$cD$lMAPdir/$sm$bamDepthFsuffix";
	my $ofasCons = "$cD/$lSNPdir/$lConsCTG";
	my $vcfFile = "$cD/$lSNPdir/$lConsVCF";
	
	#DEBUG
	
	my $secSeqTechS = "";#secondary reads..
	if (  $map{$sm}{"SupportReads"} =~ m/PB:/){$secSeqTechS = "PB" ;
	} elsif (  $map{$sm}{"SupportReads"} =~ m/ONT:/) {$secSeqTechS = "ONT" ;}
	my $seqPlatf =$map{$sm}{SeqTech}; #primary reads

	my $cmd ="";
	if ($secSeqTechS eq ""){
		#in case of only illumina:
		if ($seqPlatf eq ""){$seqPlatf = "hiSeq";} #if empty, assume hiSeq
		#checkSeqTech($seqPlatf);
		$vcf2fnaOpt = "-seqPlatform $seqPlatf -t 1 -minCallDepth $minSNPDepth -minCallQual $minSNPCallQual ";
		$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref $refFA -inVCF $vcfFile -depthF $depthFile  ";#-oCtg /dev/null " ;
	} else {
		#die;
		#in case of both PacBio and illumina:
		#$vcf2fnaOpt = "-seqPlatform $SNPIHR->{SeqTech},$SNPIHR->{SeqTechSuppl} -t 1 -minCallDepth $minDepth,$minDepth -minCallQual $minCallQual ";
		#$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref $refFA -inVCF $vcfFile,$vcfFileS -depthF $depthFile,$depthFileS ";# -oCtg $ofasCons.gz " ;
		my $vcfFileS = "$cD/$lSNPdir/$lConsVCFsup";
		my $depthFileS = "$cD$lMAPdir/$sm$bamDepthFsuffixSup";
		$vcf2fnaOpt = "-seqPlatform $seqPlatf,$secSeqTechS -t 1 -minCallDepth $minSNPDepth -minCallQual $minSNPCallQual ";
		$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref $refFA -inVCF $vcfFile,$vcfFileS -depthF $depthFile,$depthFileS  -oCtg /dev/null " ;
	}

	$cmd .= "-gff $refGFF -oGeneNT $oFNA -oGeneAA $oFAA";
	if ($append2LOG){$cmd.=" >> $SNPconsLOGs\n";
	} else {$cmd .= "\n";}
	if ($returnCmd){ #don't excecute
		return $cmd;
	}
	
	#local excecution.. probably takes forever..
	#print "$cmd\n";
	#system "echo \$SLURM_LOCAL_SCRATCH";
	system $cmd;
}

sub readGenesSample_Singl{
	#go into curSpl dir and extract all marked gene reps.. 
	#write to correct format so they can be used in phylo later
	my ($sm, $writeLink,$sttime) = @_;
	#my %subG = %{$subGHR};#$_[0]};
	
	my %subG; my %locMGScnt;
	my %locCl2G2 = %{$cl2gene2{$sm}};
			
	foreach my $gn (keys %locCl2G2){
		#put genes into hash to avoid duplicates..
		$subG{$_} = 1 foreach(@{$locCl2G2{$gn}}) ;
		#stats collection on MGS usage
		if (exists($Gene2MGS->{$gn})){
			$locMGScnt{$Gene2MGS->{$gn}}++;
		}
	}
	print scalar(keys(%subG))." genes, " . scalar(keys(%locMGScnt)). " MGS\n";
	my @histoMGScnts = values %locMGScnt;
	histoMGS(\@histoMGScnts, "Possible Bins in sample");

	
	
	
	my $sd = $sm; #this is current sample
	my $sd2 = $sd;
#	my $writeLink = 1;
	if (exists(  $map{altNms}{$sd}  )){
		$sd2 = $map{altNms}{$sd}; $replN{$sd} = $sd2;
	}
	#print "SMMM: $sd $sd2 $replN{$sd}\n";
	#check if sample in map
		#print "map s: " .scalar(keys%map)."\n";

	unless (exists ($map{$sd2}) ) {
		print STDERR "Can't find map entry for $sd\n"; #die;
		return;
		die;
	}
	my @subGKs = keys %subG;
	die "regex failed: $subGKs[0]\n" unless ($subGKs[0] =~ m/^(.*)__/);
	#find out if other samples are in the same assmblGrp..
	my @subSds = ($sd2);
	my $cAssGrp = $map{$sd2}{AssGroup};
	if (exists($AGlist{$cAssGrp})){
		@subSds = @{$AGlist{$cAssGrp}};
	}
	

	
	#print "$map{$sm}{SeqTech}\t2:$map{$sd2}{SeqTech}\t3:$map{$subSds[0]}{SeqTech}\n";
	
	#print "YY @subSds : $sd2 $sd\n";#die;
	#go into each sample ($sd3) from assembly group ($sd), that an assembly might be associated to (across multiple assemblies in assmblGrp)
	foreach my $sd3 (@subSds){
		#print "Time A: " . timeNice(time - $sttime)  . "\n";
		my $locSpace = "$locTmpDir/$sd3.cons/"; 
		#my $locSpace = "$preConDir/$sd3.cons/"; 
		
		my %locFAA; my %locFNA;my%locCSP;
		my %locMGSgenes; #keep track of genes written for each MGS..
		my $cD = $map{$sd3}{wrdir}."/";
		#print "$cD\n";
		my $rename = 0;
		$rename = 1 if ($sd2 ne $sd3);
		#print "r:$rename $sd3  (from $sd2) ";
		my $metaGD = getAssemblPath($cD);
		if ($metaGD eq ""){die "assembly not available: $cD $sd3";}
		#get NT's
		#my $tar = $metaGD."genePred/genes.shrtHD.fna";
		
		#pre-calculated, as in old MGTK versions (pre 0.69):
		my $fastaf = "$cD/$lSNPdir/$lConsFNA";
		my $fastafAA = "$cD/$lSNPdir/$lConsFAA";
		my $fastafVCF = "$cD/$lSNPdir/$lConsVCF";
		my $locForceVCF2FNA=$forceVCF2FNA;
		
		if (exists($preCompSNPs{$sd3})){ #precomputation was done before.. no need to calc then again here
			print "Found precomputed files: $preCompSNPs{$sd3}{NT}\n";
			$fastaf=$preCompSNPs{$sd3}{NT};
			$fastafAA=$preCompSNPs{$sd3}{AA};
			$locForceVCF2FNA=0;
		}
		
		#need to recreate fna/faa on the fly?? -> or does user want this?
		if ( $locForceVCF2FNA  || (! fileGZe( $fastaf ) && fileGZe($fastafVCF))){
			system "mkdir -p $locSpace";
			print "Recreating consensus fasta files on the fly..\n";
			#store these in scratch, uncompressed (much faster)
			$fastaf = "$locSpace/$sd3.cons.genes.fna";
			$fastafAA = "$locSpace/$sd3.cons.prots.faa";
			createConsFastas($cD, $sd3, $fastaf, $fastafAA, 1, 0);
		}
		#print "$fastaf\n";
		unless (-e $fastaf){
			print "\n=====================================\nCan't find nt file $fastaf -> skip sample\n=====================================\n";
			#die;
			next;
		}
		#print "Time A1: " . timeNice(time - $sttime)  . "\n";
		#print "$fastaf\n";
		#read the assemble nt and AA genes from the sample
		my $FNA = readFasta($fastaf,1,"\\s");#,\%subG);
		#my %FNA = %{$hr};
		my $FAA2 = readFasta($fastafAA,0);#,"\\s",\%subG);#read full head string
		my %FAA ;#= {};
		my %depths;
		#my $abunHR = readTabbed($cD.$abundF);
		#print "Time B: " . timeNice(time - $sttime)  . "\n";

		#my %FAA = %{$hr};
		#convert FAA hd
		my %conspSc;#read conspecific strain score from SNP consensus call..
		foreach my $k(keys %{$FAA2}){
			#$k =~ m/^(\S+)\s.*CSP=([0-9\.]+)/;
			#requires vcf2fn v 0.25
			$k =~ m/^(\S+)\sD=([0-9\.]+)\s.*CSP=([0-9\.]+)/;

			my $tmp = $1;
			if (!defined $3){print "no CSP $tmp\n";
			} else { $conspSc{$tmp} = $3;}
			if (defined($2)){$depths{$tmp} = $2;}
			$FAA{$tmp} = $FAA2->{$k};
		}
		$FAA2 = {};
		print "Time C: " . timeNice(time - $sttime)  . "\n";

		#some stats on gene extractions..
		my $missGene=0; my $foundGene=0; my $SInum=0; my $conspGen=0;my $SNPresFail=0;
		my $doubleGenes=0; my $MGStoolowGskip=0;my $missAbundance=0;
		#stats on different ways to filter genes
		my $geneLost=0; my $conSpecFail=0; my $abundFail=0; my $doubleGsFail=0;
		
		
		#3rd part: genes were read and renamed.. now write them out already here to save mem overall
		#currently takes too long in large GCs..
		foreach my $SI (@specis){ #("MGS.128"){#
			#print "$SI ";
			my @COGprios1 = @{$COGprios->{$SI}};
			#die "Can't find $SI in COGprios!\n" unless (defined(@COGprios1));#exists($COGprios->{$SI}));
			next if (scalar(@COGprios1) == 0);
			#next if (exists($ConspecificMGS{$SI}));
			#print "MGSsmplConsp{$SI}{$sd3}\n";
			#if (exists($MGSsmplConsp{$SI}{$sd}) ||  exists($MGSsmplConsp{$SI}{$sd3} )){next;}#print " DIIIIIIIIIIIIIIIID\n\n"; next;}
			
			my $OCstr=""; my $OFstr = ""; my $OAstr = ""; my $OLstr = "";
			my $locConSpecGen=0; my $accAbu=0; my $Gtrials=0; my $LmissG=0; my $doubleCntL=0; 
			 my $LmuissAbu=0;
			#get actual gene & gene2assmblname
			my @genes2 = (); #stores semi-final list of genes
			my @abunGs = (); #abundance vector of genes
			my %curcgs ;
			my %linkStr; #temp storage for links to gene cat etc of catalogues genes
			#my $curGcnt=0;
			my $MGSgcnt = scalar(@COGprios1);
			
			# 1: decide which gene to use in case of multiple COGs -> this is no longer needed, each COG should be represented by the first listed gene per MGS (decided in prior routines)
			
			foreach my $cog (@COGprios1){ 
				next if ($cog eq "");
				
				my $tar = "";my @genes ; 
				$tar = $SIgenes->{$SI}{$cog};
				next unless (exists($locCl2G2{$tar}));
				@genes = @{$locCl2G2{$tar}};
				#my $bestCOGcnt=0; next unless (scalar(@genes) > $bestCOGcnt);$bestCOGcnt = scalar(@genes);
				

				my $curG = "";my $maxAB =0;my $bestAB=100000;
				#if (1 ){#@genes > 1){
				#$doubleGenes++ if (@genes > 1); #95%gene is represented by >1 gene in sample.. potentially conspecific
					#if several genes: select most abundant from $abunHR->{}
				my $gX ; my $nonZeroCnt=0;
				foreach $gX (  @genes ){
					$Gtrials++;
					next if ($gX eq "");
					if ( !exists($FAA{$gX})){ #(exists($gene2genes{$gX}) && !exists($FAA{$gene2genes{$gX}} )) && 
						$LmissG++; 
						#die "$gene2genes{$gX} $gX\n";
						next;
					}
					#my $gX2 = $gX;if (exists($gene2genes{$gX}) && exists($FAA{$gene2genes{$gX}}) ){$gX2 = $gene2genes{$gX} ;} 
					#my $abundLoc = $abunHR->{$gX};
					my $abundLoc = $depths{$gX};
					if (defined($abundLoc)){
						$nonZeroCnt++ if ( $abundLoc > 0 );
						$curG = $gX ;
						$bestAB = $abundLoc;
					} else {
						$LmuissAbu++;
					}
				}
				next if ($curG eq "");
				#first layer filtering: gene represented by too many multi occurring genes?
				$doubleCntL++ if ($nonZeroCnt>1);
				#second layer: conspecific SNPs detected??
				if ($curG ne "" ){
					if (exists($conspSc{$curG})){
						if ($conspSc{$curG} > $conspecificSpThr){#too many indicators that gene is from conspecific strain
							$locConSpecGen ++ ;  $curG = "";#deactivate COG repri altogether..
						}
					} else {print "Can't find CSP for  \"$curG\" \n";}
				}
				if ($nonZeroCnt != 1 ){ #either 0 (gene not present) or >1 (too many copies) is not wanted
					$curG = ""; #deactivate COG repri altogether..
				}  
				next if ($curG eq "");
				
				push (@genes2 , $curG); 
				$curcgs{$curG} = $cog;
				#$curGcnt++;
				push(@abunGs, $bestAB);
				$accAbu += $bestAB ;
				
				#write link file, but only needs to be done once.. this avoids doing this later when the cat file is written
				
					#header used later: $ng = "$sd3$SaSe$curcgs{$gX}";
					my $laterHd = "$sd3$SaSe$cog";
					 $linkStr{$curG} ="$laterHd\t$cog\t$tar\t".scalar @genes . "\t".join(",",@genes)."\n" ;
				

			}
		

			#conditions where the MGS is NOT registered for current sample:
			#1) no genes..
			next if (scalar(@genes2) == 0 );
			#2) double gene counts too high  #3) too many conspecific genes
			$doubleGenes += $doubleCntL;
			$conspGen+=$locConSpecGen;
			$missGene += $LmissG;
			$missAbundance += $LmuissAbu;
			if (($doubleCntL/$MGSgcnt) > $multiGeneSmplMax 
					|| ($locConSpecGen/$MGSgcnt) >= $conspGeneSmplMax 
					# || (($locConSpecGen+$doubleCntL)/$MGSgcnt) >= (($conspGeneSmplMax +$multiGeneSmplMax)/2)
			){
				#print "$doubleCntL/$MGSgcnt > $multiGeneSmplMax\n";
				push(@{$ConspecificMGS{$SI}}, "$sd3" ); 
				if (($doubleCntL/$MGSgcnt) > $multiGeneSmplMax){$doubleGsFail++;}
				elsif (($locConSpecGen/$MGSgcnt) >= $conspGeneSmplMax){$conSpecFail++;}
				#elsif ((($locConSpecGen+$doubleCntL)/$MGSgcnt) >= (($conspGeneSmplMax +$multiGeneSmplMax)/2)){$doubleGsFail++;$conSpecFail++;}
				next;
			}
			

			my $quan10 = quantileArray(0.1,@abunGs);
			my $quan90 = quantileArray(0.9,@abunGs);
			my @genes3=();
			#filter genes based on their abundance.. only want to include approx same abundant genes
			for (my $i=0;$i<scalar(@abunGs);$i++){
				if ($abunGs[$i] <= (0.9*$quan10) || $abunGs[$i] >= (1.1*$quan90)){
					$abundFail++; next;
				}
				push (@genes3, $genes2[$i]);
			}
			
			#print "SIZE: BF:" . scalar(@genes2) . " AF:" . scalar(@genes3) . ":: $quan10, $quan90, $accAbu, $curGcnt\n";
				
			#die;
			
			if (scalar(@genes3)< $MGStoolowGsThr){
				$MGStoolowGskip++;next;
			}
				
			#now write MGS into local temp storage for later tree building..
			my $locCnt=0;
			foreach my $gX (  @genes3 ){
				print STDERR "Could not find \"$gX\" gene\n" unless (exists($FAA{$gX}));
				my $strCpy = ""; $strCpy = $FAA{$gX};# if (exists($locFAA{$gX}));
				my $AAlen = length($strCpy);
				if ($AAlen == 0){$SNPresFail++; next;}
				my $num1 = $strCpy =~ tr/[\-Xx]//;
				if ($num1 >= ($AAlen-1)){ $SNPresFail++; next;}
				if ($locCnt >= $maxNGenes){ next;}
				
				#write gene out
				my $ng = "$sd3$SaSe$curcgs{$gX}"; #must contain 2 informations: 1)sampleID 2)COG 
				#die;
				$OFstr .= ">$ng\n$FNA->{$gX}\n"; #FNA
				$OAstr .= ">$ng\n$strCpy\n"; #FAA
				$locCnt++;
				#add to category for later..
				$OCstr .= "$SI\t$curcgs{$gX}\t$sd3\t$ng\n";
				#$SIcat{$SI}{$cog}{$sd3} = $ng;
				#$genesWrite{$SI}++;
				
				if ($writeLink){
					$OLstr .= $linkStr{$gX}
				}

			}#
			
			
			if ($OFstr eq "" || $locCnt < $MGStoolowGsThr){ #5 genes is really too little to be considered valid as good strain rep..
				$MGStoolowGskip++;
				#delete $locMGSgenes{$SI};
				next;
			}
			$locMGSgenes{$SI} = $locCnt;
			
			if (!exists($OAstrH{$SI})){#set up base strings
				$OAstrH{$SI} = "";$OFstrH{$SI} = "";$OLstrH{$SI} = "";$OCstrH{$SI} = "";
			}
			#save in tmp hash (faster than opening bunch of files..
			$OAstrH{$SI} .= $OAstr;$OFstrH{$SI} .= $OFstr;$OLstrH{$SI} .= $OLstr;$OCstrH{$SI} .= $OCstr;
			$SInum ++ if ($locCnt>0);
			$foundGene+=$locCnt;
			#clenup tmp
			system "rm -rf $locSpace";
		}
		print "Time D: " . timeNice(time - $sttime)  . "\n";

		my @genesPmgs = values %locMGSgenes; 	@genesPmgs = sort { $a <=> $b}  @genesPmgs;
		histoMGS(\@genesPmgs,"Detected Bin Genes:");
		
		print "$sd3 - Missed/MissAbund/lost/abundFilterFail/SNPresFail Gs: ${missGene}/${missAbundance}/${geneLost}/${abundFail}/$SNPresFail\tConspecGs/doublGs/consMGS/failcMGS: ${conspGen}/${doubleGenes}/$conSpecFail/$doubleGsFail\tFoundGs: $foundGene/". scalar %locFAA . "\tMGS/skipped MGS: ${SInum}/$MGStoolowGskip\t";
		print "GperMGS (median,mean): " . median(@genesPmgs) . "/". int(mean(@genesPmgs)+0.5);#int($foundGene/$SInum) if ($SInum);
		print "\n";
	}
}


