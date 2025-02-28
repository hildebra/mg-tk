#!/usr/bin/env perl
#perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree4.pl -fna /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFNAs.fna -aa /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFAAs.faa -cats /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//categories4ete.txt -outD /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2/testMSA/ -cores 12 -useEte 0 -NTfilt 0.8 -runIQtree 0 -calcDistMat 1 -continue 0
#perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree4.pl -fna /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFNAs.fna -aa /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFAAs.faa -cats /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//categories4ete.txt -outD /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2/tesssst/ -cores 12 -useEte 0 -NTfilt 0.8 -NonSynTree 0 -SynTree 0 -runRAxML 0 -runGubbins 0
#perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec//secScripts/phylo/buildTree4.pl -fna /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat///intra_phylo//MB2bin314//allFNAs.fna -aa /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat///intra_phylo//MB2bin314//allFAAs.faa -smplSep '\|' -cats /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat///intra_phylo//MB2bin314//all.cat -outD /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat///intra_phylo//MB2bin314/  -runIQtree 1 -runFastTree 0 -cores 20  -AAtree 0 -bootstrap 000 -NTfiltCount 300 -NTfilt 0.05 -NTfiltPerGene 0.7 -GenesPerSpecies 0.2 -runRaxMLng 0 -minOverlapMSA 2 -SynTree 0 -NonSynTree 0 -MSAprogram 2 -continue 1 -AutoModel 0 -iqFast 1 -superTree 0 -outgroup MB2bin720 -superCheck 1
# perl /g/scb/bork/luetge/pangenomics/speciation/dNdS/scripts/buildTree4_mod2.pl -fna /g/scb/bork/luetge/pangenomics/speciation/dNdS/fasta/allOrtho_freeze11_cluster_10.fna -aa /g/scb/bork/luetge/pangenomics/speciation/dNdS/fasta/allOrtho_freeze11_cluster_10.faa -cats /g/scb/bork/luetge/pangenomics/speciation/dNdS/catFiles/freeze11_cluster_10_categories4MSA.txt -outD /g/scb/bork/luetge/pangenomics/speciation/dNdS/outFiles/test/ -cores 12 -useEte 0 -NTfilt 0.8 -NonSynTree 0 -SynTree 0 -runRAxML 0 -runGubbins 0 -runLengthCheck 0 -runDNDS 0 -genesToPhylip 0 -continue 1 -runFastgear 0 -runFastGearPostProcessing 1 -clustername cluster_10

#ARGS: ./buildTree.pl -fna [FNA] -faa [FAA] -cat [categoryFile] -outD [outDir] -cores [CPUs] -useEte [1=ETE,0=this script] -NTfilt [filter]
#versions: ver 2 makes a link to nexus file formats, to be used in MrBayes and BEAST etc
#8.12.17: added mod3 from Mechthild
#2.1.20: rewrite of workflow to extend superTrees to superCheck
#version 5 added: hyphy fubar, R scripts for theta, guidance2
#2.1.25: v5.02: reduced threshold for including genes from MGS
#14.2.25: v5.03: added treshrink

use warnings;
use strict;
#use threads ('yield','stack_size' => 64*4096,'exit' => 'threads_only','stringify');
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(gzipopen systemW readFasta readFastHD writeFasta convertMSA2NXS quantile);
use Mods::phyloTools qw(MSA filterMSA getTreeLeafs calcDisPos2 runRaxML runRaxMLng runQItree 
			runFasttree fixHDs4Phylo getGenoGenes getFMG readFMGdir );

use Getopt::Long qw( GetOptions );
#use Mods::ext::TreeIO;
#use Mods::IO::MaybeXS qw(encode_json decode_json);
use Mods::IO::PP qw (decode_json);
#use JSON qw( decode_json ); 
use Data::Dumper;
use Mods::math qw (medianArray avgArray meanArray);


sub convertMultAli2NT;
sub mergeMSAs;
sub synPosOnly;
sub calcDisPos;#gets only the dissimilar positions of an MSA, as well as %id similarity
sub calcDisPos2;#de novo aligns pairwise via vsearch and calcs id (iddef 2)
sub calcDiffDNA;
sub selecAnalysis;
sub runFastgear;
sub mergePids; sub WattTheta;
sub singleGeneMSAprocess;
sub pruneTree;
sub prepGenoDirs;
sub createTreeOpt;
sub treePresent;

my $doPhym= 0;
my $version = 5.03;

my $pal2nal = getProgPaths("pal2nal"); #"perl /g/bork3/home/hildebra/bin/pal2nal.v14/pal2nal.pl";
#die $pal2nal;
my $fasta2phylip = getProgPaths("fasta2phylip_scr");
my $phymlBin = getProgPaths("phyml");
my $trimalBin = getProgPaths("trimal");
my $pigzBin  = getProgPaths("pigz");
my $trDist = getProgPaths("treeDistScr");
my $msaFbin = getProgPaths("MSAfix");
my $RpogenS = getProgPaths("pogenStats");
my $eteBin = getProgPaths("ete3");
my $stBin = getProgPaths("supertree",0);


my $gubbinsBin = "/g/bork3/home/hildebra/bin/gubbins/python/scripts/run_gubbins.py";
my $pamlBin = "";#getProgPaths("codeml"); #PAML: currently unused
#my $evoConda = getProgPaths("evoEnv");#systemW "source $evoConda";



my $fastgearBin = "/g/bork3/home/luetge/softs/fastGEARpackageLinux64bit/run_fastGEAR.sh";
my $matlabBin = "/g/bork3/home/luetge/softs/matlab/v901";
my $fastgearSummaryBin = "/g/bork3/home/luetge/softs/fastGearPostprocessingLinux64bit/run_collectRecombinationStatistics.sh";
my $fastgearReconstrBin = "/g/bork3/home/luetge/softs/fastGearPostprocessingLinux64bit/run_startAncestryReconstruction.sh";
my $fastgearReorderBin = "/g/bork3/home/luetge/softs/fastGearPostprocessingLinux64bit/run_reorderMultipleGenes.sh";
my $partiExt=".partition.RAXML";

#die "TODO $trimalBin\n";
#trimal -in /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/COG0185.faa -out /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/tst.fna -backtrans /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/inMSA0.fna -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon
#some runtim options...
#my $ncore = 20;#RAXML cores
my $ntFrac =0.2; my $ntFracGene = 0.1; 
my $GeneFracPSpec = 0.1; #replacement for ntFracGene, as works also with supertrees
my $clustalUse = 2; #do MSA with clustal (1) or msaprobs (0), mafft(2), guidance2(3), MUSCLE5 (4)
my $calcDistMat = 0; #distmat of either AA or NT (depending on MSA)
my $calcDistMatExt = 0; #distmat of other AA or NT (depending on MSA), e.g. running two times an MSA
my $calcDistMatExtGo = 0;
my $treeAutoModel=1; #iqtree: choose model automatically (a bit slower)
my $fracMaxGenesFilter = 0.2;
my $fracMaxGenes90pct = 0.25; #gene cats to keep, e.g. 25% of 90th percentile


my $ntCntTotal =0; my $bootStrap=0; my $subsetSmpls = -1;
my ($fnFna, $aaFna,$cogCats,$outD,$ncore,$Ete, $filt,$smplDef,$smplSep,$calcSyn,$calcNonSyn,
			$useAA4tree,$calcDNAdiff,$tmpD ) = ("","","","",1,0,0.8,1,"_",0,0,0,0,"");

my ($continue,$isAligned) = (0,0);#overwrite already existing files?
my $outgroup="";
my $fixHeaders = 0;
my ($doGubbins,$doCFML,$doRAXML,$doFastTree, $doIQTree,$doRAXMLng) = (0,0,0,0, 1, 0);#fastree as default tree builder

#check length of fasta to avoid frameshift
my $doLengthCheck=1;
my $doDNDS=0;
my $doTheta=0;
my $mapF = "";
my $doGenesToPh=0;
my $doFastGear=0;
my $doFastGearSummary=0;
my $postFilter = "";
my $clusterName="";
my $MSAreq = 1;
my $iqFast=0;
my $minOverlapMSA = 0;
my $maxGapPerCol = 1 ;
my $minPcId = 0;
my $doSuperTree =0;
my $doSuperCheck=0;#check if tree's of single genes behave "strange"
my $gzipInput =0; my $removeMSA = 0;
my $useTreeShrink =0;


#EDBUGGIN
my $reparseHyphyJson = 1;


## parameters for pogen stats
my $selGene=0; 		#run dnds just on subset (given by genesForDNDS) of genes 
my @genesExtra;	#list with selected genes just for dnds
my @model= (0,1,2);  # codeml models: 0 -> neutral selection,  (1,2,7,8) -> test for postive selection
my @omegas = (0.3, 1.3); #try different omega values to check for convergence	
my $repeatCounts=2;	#set how often each model should be repeated to check for convergence 
my $codemlOutD=""; 
my $treeFile="";
my $genoindir = "";
my $wildcardflag = "";#"/*\.fna";
my $subsetPopgenStats = "10,20,30,100,200,500"; #maybe add later to options..
my $MSAsubsD =""; #only needed to create subsets of MSAs for hyphy etc calcs


die "no input args!\n" if (@ARGV == 0 );


GetOptions(
	"genoInD=s" => \$genoindir, #provide a dir with complete genomes, will extract FGMs and build tree between genomes (NT/AA flag)
	"wildcardflag=s" => \$wildcardflag,
	"fna=s" => \$fnFna,
	"aa=s"      => \$aaFna,
	"cats=s"      => \$cogCats,
	"outD=s"      => \$outD,
	"tmpD=s" => \$tmpD,
	"cores=i" => \$ncore,
	"superTree=i" => \$doSuperTree,
	"superCheck=i" => \$doSuperCheck,
	"fixHeaders=i" => \$fixHeaders, ## fix the fasta headers, if too long or containing not allowed symbols (nwk reserved)
	"useEte=i"      => \$Ete,
	"NTfilt=f"      => \$filt,
	"NTfiltPerGene=f"      => \$ntFracGene,
	"GenesPerSpecies=f" => \$GeneFracPSpec,
	"fracMaxGenes90pct=f" => \$fracMaxGenes90pct,
	"NTfiltCount=i" => \$ntCntTotal,
	"smplDef=i"	=> \$smplDef, #is the genome somehow quantified with a delimiter (_) ?
	"smplSep=s" => \$smplSep, #set the delimiter
	"outgroup=s"	=> \$outgroup,
	"AAtree=i" => \$useAA4tree,
	"MSAprogram=i" => \$clustalUse, #(0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5
	"minOverlapMSA=i" => \$minOverlapMSA, #min overlap in MSA columns, in order to retain column
	"maxGapPerCol=f" =>\$maxGapPerCol, #same as minOverlapMSA, but for MSAfix and %of gaps allowed in a column
	"calcDistMat=i" => \$calcDistMat,
	"calcDistMatExt=i" => \$calcDistMatExt,
	"calcDiffDNA=i" => \$calcDNAdiff,
	"minPcId=f" => \$minPcId, #sequence is filtered from data, unless the average minPcId is >= $minPcId
	"SynTree=i"	=> \$calcSyn,
	"NonSynTree=i"	=> \$calcNonSyn,
	"continue=i" => \$continue,
	"bootstrap=i" => \$bootStrap,
	"subsetSmpls=i" => \$subsetSmpls,
	"postFilter=s" => \$postFilter, # "," sep list of zorro,guidance2,macse
	"rmMSA=i" => \$removeMSA, #to save diskspace
	"gzInput=i" => \$gzipInput, #to save diskspace
	"isAligned=i" => \$isAligned,
	"runRAxML=i" => \$doRAXML,
	"runRaxMLng=i" => \$doRAXMLng,
	"runFastTree=i" => \$doFastTree,
	"treeShrink=i" => \$useTreeShrink,
	"runIQtree=i" => \$doIQTree,
	"AutoModel=i" => \$treeAutoModel,
	"iqFast=i" => \$iqFast, #fast qiTree mode
	"runClonalFrameML=i" => \$doCFML,
	"runGubbins=i" => \$doGubbins,
	"runLengthCheck=i" => \$doLengthCheck,		#check that sequence length can be divided by 3
	"runDNDS=i" => \$doDNDS,			#run dNdS analysis
	"runTheta=i" => \$doTheta,
	"genesForDNDS=s{,}" => \@genesExtra,		#list with selected genes just for dnds
	"DNDSonSubset=i" => \$selGene,			#run dnds just on subset (given by genesForDNDS) of genes
	"codemlRepeats=i" => \$repeatCounts,
#	"treefileCodeml=s" => \$treeFile,
	"outDCodeml=s"=> \$codemlOutD,
	"genesToPhylip=i" => \$doGenesToPh,	
	"runFastgear=i" => \$doFastGear,
	"runFastGearPostProcessing=i" => \$doFastGearSummary,
	"map=s" =>\$mapF,
	"clustername=s" => \$clusterName,
) or die("Error in command line arguments\n");


print "BuildTree 5 script v$version\nOutDir: $outD\n";
##DEBUG
#die if ($useAA4tree == 0);





######### indir
if ($genoindir ne ""){
	my $genoindir2 = $genoindir;
	if (!-d $genoindir2){ $genoindir2 =~ s/[^\/]+$//;}
	if ($outD eq ""){$outD = $genoindir2."/phylo/";}
}

##### setup dirs
$codemlOutD = "$outD/codeml" if ($codemlOutD eq "");;	

$tmpD = $outD."/tmp/" if ($tmpD eq "");
my $treeD = "$outD/phylo/";#raxml, fasttree, phyml tree output dir

my $MsaD = "$outD/MSA/";
if ($removeMSA){
	$MsaD = "$tmpD/MSA_$clusterName/";
}

$MSAsubsD = "$MsaD/clnd/";


if ($subsetSmpls >0){
	$MsaD =~ s/\/$/_S$subsetSmpls\//;
	$treeD =~ s/\/$/_S$subsetSmpls\//;
}
######

if ($clustalUse == 0){print "Warning:  MSAprobs with trimal gives warnings (ignore them)\n";}
if ($doCFML && !$doRAXML){die "Need RaxML alignment, if Clonal fram is to be run..\n";}

if ($aaFna eq "" || $useAA4tree){	$calcSyn=0;$calcNonSyn=0;}
if ($filt <1){$ntFrac=$filt; print "Using filter with $ntFrac fraction of nts\n";}
if ($outgroup ne ""){print "Using outgroup $outgroup\n";}
if ($bootStrap>0){print "Using bootstrapping in tree building\n";}
#if (($calcDistMat || $calcDistMatExt) && $isAligned || !$clustalUse){die"Can't calc distance mat, unless clustalO is being used for MSA\n";}
#else {$ntCntTotal = $filt;}

$MSAreq = 0 if (!$doFastTree && !$doRAXML && !$doRAXMLng && !$doCFML && !$doGubbins && !$doIQTree);

system "mkdir -p $tmpD" unless (-d $tmpD);
my $cmd =""; my %usedGeneNms;


my $outD_clust = "";
if($clusterName eq ""){$outD_clust = "$outD/MSA_FG";}

#------------------------------------------
#sorting by COG, MSA & syn position extraction
if ($Ete){
	$cmd = "$eteBin build -n $fnFna -a $aaFna -w clustalo_default-none-none-none  -m sptree_raxml_all --cpu $ncore -o $outD/tree --clearall --nt-switch 0.0 --noimg  --tools-dir /g/bork3/home/hildebra/bin/ete/ext_apps-latest"; #--no-seq-checks
	$cmd .= " --cogs $cogCats" unless ($cogCats eq "");
	print "Running tree analysis ..";
	print $cmd."\n";
	die;
	system $cmd . "> $outD/tree/ETE.log";
	print " Done.\n$outD/tree\n";
	exit(0);
}


#general routine, starting with MSA
#followed by merge of MSA / NT -> AA conversion, MSA filter etc
#followed by tree building

if ($fixHeaders){
	if ($cogCats ne ""){die"implement fix hds for cats\naborting\n";}
	$aaFna=fixHDs4Phylo($aaFna);$fnFna = fixHDs4Phylo($fnFna); 
}


system "rm -fr $treeD $MsaD" if (!$continue);
system "mkdir -p  $MsaD" unless(-d "$MsaD");
system "mkdir -p  $treeD/" unless(-d "$treeD");
my $multAli = "$MsaD/MSAli.fna";
my $multAliSyn = $multAli.".syn.fna";
my $multAliNonSyn = $multAli.".nonsyn.fna";
my @theRealMSAs;
my $partiFile="";#partitioning for multi gene MSAs
my %specList; #list of species (without _COG00012 tag);
my %samples; 
my $MSAcat = "$MsaD/MSAcat.fna";

prepGenoDirs($genoindir);

#prep tree Options
my $tOhr = createTreeOpt($multAli,"allsites","",0,"");
my %Tree1 = %{$tOhr};
my $tOhrNSun = createTreeOpt($multAliNonSyn,"nonsyn","",0,$Tree1{nwk});
my $tOhrSyn = createTreeOpt($multAliSyn,"syn","",0,$Tree1{nwk});



#DEBUG
#mergePids("$outD/MSA/",40, "NT") ;die;
#my $tmp = "/g/bork5/hildebra/results/TEC2/v5/T2dphylo/rDNA2/fullGenomes/ini16S.fna";
#calcDisPos2($tmp,"$outD/MSA/percID_syn.txt",1); die;


my @MSAs; my @MSA_AA; my @MSAsSyn; my @MSAsNonSyn;#full MSAs and MSAs with syn / nonsyn pos only
my @MSrm; 
my %FAA ; my %FNA ; my @geneList; my @geneListF;
my $doMSA = 1;
my $treesDone =0; $treesDone=1 if (treePresent($tOhr) && treePresent($tOhrNSun) && treePresent($tOhrSyn));
my $calcMSA=1; $calcMSA=0 if ($treesDone || -e $multAli || !$continue);
if (!$treesDone){#cleanup, avoid checkpoints..
	system "rm -f $treeD/*";
}
$doMSA =0 if ($isAligned || (
			$continue && (-e $multAli || !$calcMSA) && (-e $multAliSyn ||!$calcSyn)&& (-e $multAliNonSyn ||!$calcNonSyn)) );  ## checks if MSA already exists
#die "$doMSA $calcMSA $treesDone $continue $multAli\n";
#die "$doDNDS\n";
#my @xx = keys %FAA; die "$xx[0] $xx[1]\n$FAA{HM29_COG0185}\n";
if ($isAligned){
	if (-e $fnFna){
		system "ln -s $fnFna $multAli";
		$useAA4tree = 0;
		print "Using NT sequences to build tree..\n\n";
	}
} elsif (!$doMSA && $cogCats ne ""){
	fillGeneList($cogCats);
} elsif ($doMSA && $cogCats ne ""){
	my $hr = readFasta($aaFna,1); my %FAA = %{$hr};
	#die "$aaFna\n";
	#die $FAA{"1214150.PRJNA171417_05DKX"}."\n";
	$hr = readFasta($fnFna,1); my %FNA = %{$hr};
	#die $FNA{"1214150.PRJNA171417_05DKX"}."\n";		
	print "ReadFasta\n";

	############# test length of fna sequences can be divided by 3 ##############################
		
	if($doLengthCheck){
		my $FNAseq; my $length; my $div;
		while (($FNAseq) = each (%FNA)){
			$length = length($FNA{$FNAseq});
			$div = $length/3;
			print "AA seq can not divided by 3 in $FNAseq\n" if($div =~ /\D/);
		}
	}

	############# test if enough seq in Sample to add to tree (avoids confusion in MSA) ##############################

	my $cnt = -1;  # line count in cats file
	my $ogrpCnt=0;#my %genCats; 
	my %totalNTs;#overall nt counts (not N)
	my %charCnts;#per gene NT counts
	my %maxNtCnt;# per gene max NTs observed
	my %meanNTcnt; #probably more sensible to use this re overpredicting gene length etc
	my %qtl90NTcnt;
	my @genesPerCat;
	my $geneTooShort = 0; #count genes with too little NTs in gene..
	my $geneTooLong = 0;
	
	my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");
	#open my $xI,"<$cogCats" or die "Can't open cogcats $cogCats\n";
	chomp(my @linesCats = <$xI>);
	close $xI;
	#first cleanup of cat file..
	my @linesCats2; my @linesCats3;
	foreach (@linesCats){ #check first some parameters..
		$cnt++; my @spl = split /\t/;
		if ($spl[0] =~ m/^#/){shift @spl;}
		@spl = grep !/^NA$/, @spl;#remove NAs
		my @spl2;
		#$genesPerCat[$cnt] = scalar(@spl) ;
		my @geneLs;
		$spl[0] =~ m/^(.*)$smplSep(.*)$/;	my $sp = $1;my $gene = $2;
		foreach my $seq (@spl){### $seq = genomeX_NOGY
			$seq =~ m/^(.*)$smplSep(.*)$/;	$sp = $1;
			die "can't find AA seq $seq\n" unless (exists ($FAA{$seq}));
			die "can't find fna seq $seq\n" if (!exists ($FNA{$seq}) && !$useAA4tree);
			#print "$MFAA{$curK}\n";			#my $ss = $FAA{$seq}; 			#filter per sequence 
			my $num1 = $FAA{$seq} =~ tr/[\-Xx]//;
			my $geneL = (length( $FAA{$seq})-$num1);#AA length
			$charCnts{$sp}{$seq} = $geneL;
			push(@geneLs, $geneL);
		}
		my $qtl = quantile(0.9,@geneLs);#values(%{$charCnts{$sp}}));
		#print "Q$qtl $gene @geneLs\n";
		$qtl90NTcnt{$gene} = $qtl;#
		foreach my $seq (@spl){
			$seq =~ m/^(.*)$smplSep(.*)$/;	my $sp = $1;
			#quantile(0.8,values(%{$charCnts{$sp}}));
			if ( $charCnts{$sp}{$seq} >= ($qtl90NTcnt{$gene}  * $ntFracGene)){
				push(@spl2, $seq);
				$geneTooLong++;
			} else {
				$geneTooShort++;
			}
		}
		push(@linesCats2,\@spl2);
		#has to work with what is actually there, not what could have been..
		$genesPerCat[$cnt] = scalar(@spl2);
		#die;
	}
	#die;
	my $GenesQtl90 = quantile(0.9,@genesPerCat);
	my $GenesQtl50 = quantile(0.5,@genesPerCat);
	$cnt=-1;
	foreach my $aRef (@linesCats2){ #remove genes with just too few genes..
		$cnt++; my @spl = @{$aRef};
		if (@spl >= (($GenesQtl90 * $fracMaxGenes90pct) ) ){ #$GenesQtl50 || 
			push(@linesCats3,\@spl);
		}
		#print @spl . " ";
	}
	
	print "\n\n-----------------  Prefilter  ------------------\n";
	print "Remaining gene cats: ". scalar(@linesCats3) . "/" . scalar(@linesCats)."; Removed $geneTooShort/$geneTooLong genes < $ntFracGene qtl90 gene length\n";
	print "Warning:: Size linesCats3:: " . @linesCats3 . " linesCats2:: " .@linesCats2 ."\n$GenesQtl50 || $GenesQtl90 * $fracMaxGenes90pct\n" if (@linesCats3 < 20);
	@linesCats2 = (); #make space..
	$cnt=-1;
	foreach my $aRef (@linesCats3){
		$cnt++; my @spl = @{$aRef};
		$spl[0] =~ m/^(.*)$smplSep(.*)$/;#my @spl3 =($1,$2);
		#my @geneLgt; 
		my $gene = $2;
		foreach my $seq (@spl){
			$seq =~ m/^(.*)$smplSep(.*)$/;#my @spl3 =($1,$2);
			my $sp = $1;
			die "Wrong gene in $seq, expected $gene!\n" if ($2 ne $gene);
			$specList{$sp} ++;
			#my $seq2 = $seq;
			if (!exists($maxNtCnt{$gene})){
				$maxNtCnt{$gene} = $charCnts{$sp}{$seq};
			} elsif (
				$charCnts{$sp}{$seq} > $maxNtCnt{$gene}){ $maxNtCnt{$gene} = $charCnts{$sp}{$seq};
			}
			#push(@geneLgt,$charCnts{$sp}{$seq});
			$meanNTcnt{$gene}+=$charCnts{$sp}{$seq};
		}
		$meanNTcnt{$gene} /= scalar(@spl);
		# sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
		#$qtl90NTcnt{$gene} = (sort { $a <=> $b }( @geneLgt)) [ int($#geneLgt*0.9) ];
		#print "DEB: $gene $meanNTcnt{$gene} $qtl90NTcnt{$gene}\n";
		#second round, do some prefiltering already, now that we have qtl and mean gene size
		foreach my $seq (@spl){
			$seq =~ m/^(.*)$smplSep(.*)$/;#my @spl3 =($1,$2);
			my $sp = $1;
			#first check if gene gets removed
			#next if ( ($charCnts{$sp}{$seq} < $maxNtCnt{$gene} ) * $ntFracGene);
			#next if ( $charCnts{$sp}{$seq} < ($qtl90NTcnt{$gene} * $ntFracGene));
			$totalNTs{$sp} += $charCnts{$sp}{$seq};
		}
	}
	#die "@genesPerCat\n$GenesQtl90\n".$fracMaxGenes90pct*$GenesQtl90."\n" ;
	my @specs = keys %specList;
	#print "specs:: @specs\n";
	die "No species left after filtering!!\n" if (@specs == 0);
	
	
	my $maxGenes=0; my $maxNtCntTotal=0; #my @allNTcnts;
	foreach my $sp (@specs){
		if ($specList{$sp}>$maxGenes){$maxGenes = $specList{$sp};}
		if (!exists($totalNTs{$sp})){$totalNTs{$sp}=0;}
		if ($maxNtCntTotal< $totalNTs{$sp}){$maxNtCntTotal = $totalNTs{$sp};}
		#push(@allNTcnts,$totalNTs{$sp});
	}
	my $qtl90NTcntAll =  quantile(0.9,values(%totalNTs));#@allNTcnts);#(sort { $a <=> $b }( @allNTcnts)) [ int($#allNTcnts*0.9) ];
	my $qtl95Genes = quantile(0.95,values(%specList));
	my $qtl90Genes = quantile(0.9,values(%specList));

	my %smplsRmvd; my $tooFewGenes=0;my $tooFewNTs=0;my $tooFewNTs2=0; my $specsRemain = 0;
	#print "Samples removed due to low gene presence:\n";
	my $OGfnd=0;
	foreach my $sp (@specs){
		my $isOG=0;  if ($outgroup ne "" && $outgroup eq $sp){$isOG = 1;$OGfnd++;}
		
		my $NTfilter = 0; $NTfilter =1 if ( $totalNTs{$sp} < ($qtl90NTcntAll * $ntFrac));
		my $factor = 3 ; #AA counts..
		
		my $NTfilter2 =  0;$NTfilter2 = 1 if ( ($totalNTs{$sp}/$factor) <  ($ntCntTotal) ); #totalNTs are actually AA from FAA{}, $ntCntTotal are in nt's
		my $NTlengFilt = 0; $NTlengFilt =1 if ($specList{$sp} <  ($qtl90Genes * $GeneFracPSpec) );

		if (!$isOG && ($NTlengFilt || $NTfilter || $NTfilter2) ){
			$smplsRmvd{$sp}=1;
			$tooFewNTs++ if ($NTfilter);
			$tooFewNTs2++ if ($NTfilter2);
			$tooFewGenes++ if ($NTlengFilt);
			
			#print " $sp:$specList{$sp}:$totalNTs{$sp}; ";
		} else {
			$specsRemain ++;
		}
	}
	
	if ($OGfnd == 0 && $outgroup ne ""){
		die "could not find outgroups in sequence set!\n$outgroup\n$fnFna\n";
	}
	#############################################################################################

	
	print "Per species: MaxGenes: $maxGenes, Qtl90Genes: $qtl90Genes, MaxAA: $maxNtCntTotal, Qtl90 NTs: $qtl90NTcntAll\n";
	print "Species/Smpls removed: <NTs($ntFrac,$ntCntTotal):$tooFewNTs,$tooFewNTs2 ; <genes($GeneFracPSpec):$tooFewGenes\n";
	print "Remaining Smpls/Strains: $specsRemain/".scalar(keys%specList)."\n";
	print "------------------------------------------------\n";
	#die "$maxGenes\n";
	@linesCats = (); #empty array


	#die;
	$cnt=-1; #line counter
	foreach my $aRef (@linesCats3){#go over each gene category, building MSA for each
		$cnt++; my @spl = @{$aRef};
		if (@spl ==0){print "No categories in cat file line $cnt\n";next;}
		if ($spl[0] =~ m/^#/){shift @spl;}
		$spl[0] =~ m/^(.*)$smplSep(.*)$/;
		my $gene = $2;
		my @spl2 = ($1,$2);#split /$smplSep/,$spl[0] ;	
		#die "@spl\n";		
		my $ogrGenes = "";
		if ($outgroup ne ""){
			foreach my $seq (@spl){
				#my @spl3 = split /$smplSep/,$seq ; 
				$seq =~ m/^(.*)$smplSep(.*)$/;#my @spl3 =($1,$2);
				if ($1 eq $outgroup){$ogrGenes = $seq; $ogrpCnt ++ ;last;}
			}
		}
		
		#die "$spl2[0]\t$spl2[1]\n";
		

		die "Double gene name in tree build pre-concat: $spl2[0] $spl2[1]\n" if (exists($usedGeneNms{$spl2[1]}));
		$usedGeneNms{$spl2[1]} = 1;
		#die "@spl\n";
		my $tmpInMSA = "$tmpD/inMSA$cnt.faa";
		my $tmpInMSAnt = "$tmpD/inMSA$cnt.fna";
		my $tmpOutMSA2 = "$MsaD/$spl2[1].$cnt.faa";
		my $tmpOutMSA = "$MsaD/$spl2[1].$cnt.fna";
		open O,">$tmpInMSA" or die "Can;t open tmp faa file for MSA: $tmpInMSA\n";
		open O2,">$tmpInMSAnt" or die "Can;t open tmp fna file for MSA: $tmpInMSAnt\n";
		my $seqType = "AA";my $seqTypeOth = "NT";
		my $seqLength = 0; my $numSeq =0;
		
		#1st: collate sequences
		#do here already per gene length check .. probably better for alignment
		foreach my $seq (@spl){### $seq = genomeX_NOGY
			$seq =~ m/^(.*)$smplSep(.*)$/;#my @spl2 =($1,$2);
			my $sp = $1;
			next if (exists($smplsRmvd{$1}));
			if ($specList{$sp} <  ($qtl90Genes * $GeneFracPSpec) ){die "buildTree: GeneFracPSpec maxGenes shouldn't be here!\n";}
			my $seq2 = $seq;
			#just for this singular case applying..
			#next if ( $charCnts{$sp}{$seq} < ($qtl90NTcnt{$gene}  * $ntFracGene));  #maxNtCnt{$gene}
			if (!$smplDef){#create artificial head tag
				#TODO.. don't need it now for tec2, since no good NCBI taxid currently...
			}
			$samples{$1} = 1; #$genCats{$spl2[1]} = 1; 
			$FAA{$seq} =~ s/\*//g if (!$clustalUse);
			print O ">$seq2\n$FAA{$seq}\n";
			if (!$useAA4tree){
				$FNA{$seq} =~ s/-//g;
				print O2 ">$seq2\n$FNA{$seq}\n";
				$seqLength += length($FNA{$seq});
			} else {
				$seqLength += length($FAA{$seq});
			}
			$numSeq++;
		}
		close O;close O2;
		#done, samples are in O2
		if ($numSeq <= 3){ #actually pretty useless, no tree can be built from this, so just rm this one...
			system "rm -f $tmpInMSA $tmpInMSAnt";#
			next;
		}

		$seqLength /= @spl;
		#print "$tmpInMSA,$tmpOutMSA2\n";
		$tmpOutMSA2 = MSA($tmpInMSA,$tmpOutMSA2,$ncore,$clustalUse,$continue,$numSeq);
		$tmpOutMSA2 = filterMSA($tmpInMSA,$tmpOutMSA2,$ncore,$postFilter,$useAA4tree);
		if ($tmpOutMSA2 eq ""){ #filter failed...
			print "skipped protein, too many bad positions\n";
			system "rm -f $tmpOutMSA2 $tmpOutMSA $tmpInMSA $tmpInMSAnt";
			next;
		}
		
		#dist mat related
		my $tmpDMatOth = "$MsaD/${seqTypeOth}_clustalo_percID_${cnt}_".int($seqLength).".txt";
		my $tmpDMat = "$MsaD/${seqType}_clustalo_percID_${cnt}_".int($seqLength).".txt";
		my $inFastaOth = $tmpInMSAnt;
		my $percIDhr; my $avgID; my $pIDsmplhr;
		if ($calcDistMat){ #for dmat: calc each gene spearately and merge scores later
			($avgID,$pIDsmplhr,$percIDhr)  = calcDisPos2($tmpInMSA,$tmpDMat,0,$ncore,$tmpD);
		#			if ($calcDistMatExt && -e $inFastaOth){
			if (-e $inFastaOth){ 
				($avgID,$pIDsmplhr,$percIDhr) = calcDisPos2($inFastaOth,$tmpDMatOth,1,$ncore,$tmpD);
				$calcDistMatExtGo = $avgID;
			} else { $calcDistMatExtGo = 0;}
		}
		
		#die "@MSAs\n";
		if (!$useAA4tree){
			#this part now is all concerned about NT level things..
			convertMultAli2NT($tmpOutMSA2,$tmpInMSAnt,$tmpOutMSA);
			my ($tmpOutMSAsyn,$tmpOutMSAnonsyn) = synPosOnly($tmpOutMSA,$tmpOutMSA2,0,$ogrGenes,$calcSyn,$calcNonSyn);
			#this will not affect 4-fold only etc..
			$cmd = "$msaFbin -i $tmpOutMSA  -maskLowID -maskBorderGap -rmGapColsGreater ".$maxGapPerCol." -minGoodPosFrac 0.6\n";
			systemW $cmd;
			push (@MSAs,$tmpOutMSA);
			push (@MSAsSyn,$tmpOutMSAsyn) if ($tmpOutMSAsyn ne "");
			push (@MSAsNonSyn,$tmpOutMSAnonsyn) if ($tmpOutMSAnonsyn ne "");
			#die "@MSAs\n";
		} else {
			push (@MSA_AA,$tmpOutMSA2);
		}
		system "rm -f $tmpInMSA $tmpInMSAnt";# $tmpOutMSA2";
		push (@MSrm,$tmpOutMSA2,$tmpOutMSA);
		#die "$MSrm[1]\n";
		print "$cnt "; 
	}
	
	my $mergPIDtag = "_merge";
	mergePids("$MsaD/",$cnt, "AA",$mergPIDtag) if ($calcDistMat); #merge different percIDs
	mergePids("$MsaD/",$cnt, "NT",$mergPIDtag) if ($calcDistMat); #merge different percIDs
	
	#could be used to filter genes further, but not for now
	#pogenStatsFilter();
	if ($outgroup ne ""){print "Found $ogrpCnt of $cnt outgroup sequences\n";}
} elsif ($doMSA) {#no marker way, single gene
	my $r1; my $r2;
	#,$r1,$r2)
	$multAli = singleGeneMSAprocess($multAli)#;,\@MSAs,\@MSA_AA);
	#@MSAs = @{$r1};	@MSA_AA = @{$r2};
}

#die "@MSA_AA\n\n";
if ($calcMSA && (@MSAs == 0 && @MSA_AA == 0) ){
	my $IQtreef= "$treeD/IQtree_allsites.treefile";
	system "mkdir -p $treeD"; system "touch $IQtreef";
	print "No MSAs generated\n"; exit 0;
}

#prep final MSA file that is correct NT or AA and is merged
if (!$useAA4tree) {
	$calcSyn=0;$calcNonSyn=0;
	if ($cogCats eq ""){ #single gene case
		my ($hr,$OK) = readFasta($multAli,1); writeFasta($hr,$multAli);#complicated way to shorted headers of infile
	}
	mergeMSAs(\@MSAs,\%samples,$multAli,0,0);
	mergeMSAs(\@MSAsSyn,\%samples,$multAliSyn,1,0) if ($calcSyn);
	mergeMSAs(\@MSAsNonSyn,\%samples,$multAliNonSyn,1,0) if ($calcNonSyn);
	@theRealMSAs = @MSAs;

} else {#useAA4tree
	mergeMSAs(\@MSA_AA,\%samples,$multAli,0,1); #sames files as in @MSrm
	@theRealMSAs = @MSA_AA;
}

#phylip conversion??
if ( $doGenesToPh){ 
	my $phylipD = "$outD/phylip/";
	system "mkdir -p  $phylipD" unless(-d $phylipD || (!$doGenesToPh ));
	foreach my $MSAfn (@MSAs){
		#my $MSAfn = $tmpOutMSA;
		system "rm -f $phylipD/$MSAfn.ph*\n";
		my $cmd2 = "$fasta2phylip -c 50 $MSAfn > $MSAfn.ph\n";
		systemW $cmd2;
		push(@geneList, $MSAfn);
	}
}


#die;

#-------------------------------------------
#Tree prep phase (MSA clean up, conversion, 4fold sites etc)
#convert fasta again




#-------------------------------------------
#Supertrees, gubbins etc
#-------------------------------------------

if ($doSuperTree || $doSuperCheck){#can be for 2 reasons: 1) build actual super tree 2) quality control
	my @treeCol;
	for (my $i=0;$i<@theRealMSAs;$i++){
		print "===============>  Subtree $i  <===============\n";
		my $tOhrST = createTreeOpt($theRealMSAs[$i],"allsites",$i,1,"");
		my $trRetH = treeAtHeart($tOhrST);
		push(@treeCol,${$trRetH}{IQtreeout}.".treefile");
		if ($calcSyn){
		} 
		if ($calcNonSyn){
		}
	}
	if ($doSuperTree){
		my $outST = "$treeD/IQtree_allsites.treefile";
		my $specFile = "$treeD/IQtree_allsites.species";
		open OU,">$specFile";
		print OU join "\n",keys (%specList);
		close OU;
		my $cmd = "$stBin -s $specFile -F -o - @treeCol  | grep '\\[F01\\]' | cut -f2 -d' ' > $outST"; #-F
		#die "ST:\n$cmd\n";
		systemW $cmd;
		die $outST."\n";
	} elsif ($doSuperCheck) {
		
	}
}


if ($calcDNAdiff){
	calcDiffDNA($multAli,"$MsaD/percID_2.txt");
}


my $phyloTree = "";
if ($doGubbins){
	my $outDG = "$outD/gubbins/"; 
	system "mkdir -p $outDG" unless (-d $outDG);
	$outDG .= "GD";
	if ($continue && -e $outDG.".final_tree.tre"&& -e $outDG.".summary_of_snp_distribution.vcf"){
		print "Gubbins result already exists in output folder, run will be skipped\n";
	} else {
		system "rm -rf $outDG\n";
		my $cmdG = "source activate py3k\n";
		$cmdG .= "$gubbinsBin --filter_percentage 50  --tree_builder hybrid --prefix $outDG --threads $ncore $multAli";
		if (0){$cmdG.=" --outgroup $outgroup";}
		$cmdG.="\n";
		$cmdG .= "source deactivate py3k\n";
		systemW $cmdG;
		#die $cmdG."\n";
		print "Gubbins run finished\n";
	}
}


#distamce matrix, this is fast
#system "$trimalBin -in $multAli -gt 0.1 -cons 100 -out /dev/null -sident 2> /dev/null > $outD/MSA/percID.txt\n";
if (0 && !$useAA4tree){ #this is outdated
	calcDisPos($multAli,"$MsaD/percID.txt",1) unless(-e "$MsaD/percID.txt" && $continue);
	if ($calcSyn){
	#		system "$trimalBin -in $multAliSyn -gt 0.1 -cons 100 -out /dev/null -sident 2> /dev/null > $outD/MSA/percID_syn.txt\n";
	#	calcDisPos($multAliSyn,"$outD/MSA/percID_syn.txt",1);
	}
	if ($calcNonSyn){
		#system "$trimalBin -in $multAliNonSyn -gt 0.1 -cons 100 -out /dev/null -sident 2> /dev/null > $outD/MSA/percID_nonsyn.txt\n";
	#	calcDisPos($multAliNonSyn,"$outD/MSA/percID_nonsyn.txt",1);
	}
} elsif (0) {
	calcDisPos($multAli,"$MsaD/AA_percID.txt",1) unless(-e "$MsaD/percID.txt" && $continue);
}

#print "$theRealMSAs[0]\n\n\n";

#-------------------------------------------
#Tree building part with RaxML, IQtree, fasttree2, phyml
#-------------------------------------------

my $trRetH = treeAtHeart($tOhr );
if ($calcSyn){ #tree at syn pos
	treeAtHeart($tOhrSyn);
} 
if ($calcNonSyn){ #tree at non-syn pos
	treeAtHeart($tOhrNSun);
}
#system "rm -f $multAli.ph $multAliSyn.ph $multAliNonSyn.ph";

if ($useTreeShrink){
	my $trShr = getProgPaths("treeshrink");
	my $cmd = "$trShr -i $outD -t ${$trRetH}{IQtreeout}.treefile -q 0.05  -O TS. -f";
	print $cmd; 
	die;
}

#die "$distTree_scr -d -a --dist-output $raxD/distance.syn.txt $raxD/RXML_sym.nwk\n";


#### post tree building #######
#### compute dNdS and other popgen values ######

pogenStatsFilter();

#might require $doGenesToPh for paml??
if($doDNDS){
	system "mkdir -p  $codemlOutD" unless(-d "$codemlOutD");
#die "@geneList";
	if($selGene){@geneList=@genesExtra;	}
	#my $tmpDir = $codemlOutD."/tmp/";;
	system "mkdir -p  $tmpD" unless(-d $tmpD);	
	my $treeFile = $Tree1{nwk} if ($treeFile eq "");
	selecAnalysis(\@geneList, $treeFile, $codemlOutD, $tmpD);   

}
#if ($doTheta){ #Watterman estimator (Theta)
#	if($selGene){@geneList=@genesExtra;	}
#	WattTheta(\@geneList,$MsaD,$codemlOutD);
#}

FastGear();

system "rm -rf $MsaD" if ($removeMSA);
if ($gzipInput){
	system "$pigzBin -p $ncore  $aaFna " unless ($aaFna =~ m/\.gz$/);
	system "$pigzBin -p $ncore $fnFna  " unless ($fnFna =~ m/\.gz$/);
	system "$pigzBin -p $ncore $cogCats" unless ($cogCats =~ m/\.gz$/);
}

system "rm -rf $tmpD";
	###################### ETE ######################3

print "All done: $outD \n\n";
exit(0);













##########################################################################################
##########################################################################################





sub treePresent{
	my ($hr) = @_;
	my %treeOpts = %{$hr};
	my $ret = 1;
	if ($doFastTree){
		$ret=0 unless ($continue && -e $treeOpts{fastTrOut});
	}
	if ($doIQTree){
		my $IQtree = "$treeOpts{IQtreeout}";
		$ret=0 unless ($continue && -e "$IQtree.treefile");
	}
	if ($doRAXMLng){
		$ret=0 unless ($continue && -e $treeOpts{RAXNGtreeout});
	}
	if ($doRAXML){
		$ret=0 unless (-e $treeOpts{RAXtreeout});
	}
	return $ret;
}



sub createTreeOpt{
	my ($multF,$siteTag,$tcnt,$silent,$consTree) = @_;
	#$siteTag="allsites";
	my $isSubTree = 0;
	my $outgroupL = $outgroup;
	$isSubTree = 1 if ($tcnt ne "");
	$outgroupL = "" if ($isSubTree);
	my $partiF=$multF.$partiExt;
	$partiF="" unless (-e $partiF);
	#object to transfer options to tree (and get them back..)
	my $BStag = ""; if ($bootStrap>0){$BStag="_BS$bootStrap";}
	my %treeOpts = (inMSA => $multF,
					IQtreeout => "$treeD/IQtree${tcnt}_${siteTag}",
					ncore => $ncore,
					tcnt => $tcnt,
					outgr => $outgroupL,
					bootStrap => $bootStrap,
					useAA => $useAA4tree,
					iqtreeFast => $iqFast,
					autoModel => $treeAutoModel,
					cont => $continue,
					silent => $silent,
					partition => $partiF,
					constraintTree => $consTree,
					isSubTree => $isSubTree,
					PhymTree => "$treeD/phyml${tcnt}_${siteTag}.nwk",
					fastTrOut => "$treeD/FASTTREE${tcnt}_${siteTag}.nwk",
					RAXtreeout => "$treeD/RXML${tcnt}_${siteTag}$BStag.nwk",
					RAXNGtreeout => "$treeD/RXng${tcnt}_${siteTag}$BStag.nwk",
					);
	return \%treeOpts;
}

#core routine to calculte (start) phylo reconstruction
sub treeAtHeart{
	my ($hr) = @_;
	my %treeOpts = %{$hr};
	my $consTree = $treeOpts{constraintTree}; my $multF = $treeOpts{inMSA};
	my $silent = $treeOpts{silent}; my $tcnt = $treeOpts{tcnt};
	if ($consTree ne ""){
		die "ref tree $consTree does not exist\n" unless (-e $consTree);
		
		my $hr = getTreeLeafs($consTree);
		my %nwLfs= %{$hr};
		
		my $cntMissTree=0;
		$hr = readFasta($multF,1,"input MSA to check for constraint tree");
		my %FNA = %{$hr};
		my $unpresent=0;
		foreach my $ge (keys %FNA){
			unless (exists($nwLfs{$ge})){$unpresent++; last;}
		}
		my @genomeList;
		if ($unpresent){
			my $cntMissTree=0;
			open O,">$multF.tmp" or die "can't open MSA out $multF.tmp\n";
			foreach my $genome (keys %FNA){
				#print "$genome\n";
				unless (exists($nwLfs{$genome})){$cntMissTree++; next;}
				my $seq = $FNA{$genome};
				print O ">$genome\n$seq\n";
				push (@genomeList, $genome);
			} 
			close O;
			if ($cntMissTree>0){print "Removed $cntMissTree extra sequences from MSA not in constraint tree\n";}
			system "rm -f $multF;mv $multF.tmp $multF";
		} else {
			@genomeList = keys %FNA;
		}
		
		#my $ar = readFastHD($multF);
		my $nwkPrune = $multF."pruned.nwk";
		if (scalar(@genomeList) > 3){
			pruneTree($consTree,\@genomeList,$nwkPrune);
			$treeOpts{constraintTree} = $nwkPrune;
		} else {
			$treeOpts{constraintTree} = "";
		}

	}
	#print "cons: $treeOpts{constraintTree}\n";
	#convert MSA to NEXUS
	#convertMSA2NXS($multAli,"$multAli.nxs");
	#format conversion for raxml..
	if ($doRAXML){
		my $f = $treeOpts{inMSA};
		my $tcmd = "rm -f $f.ph*; $fasta2phylip -c 50 $f > $f.ph\n";
		systemW $tcmd;#) {die "fasta2phylim failed:\n$tcmd\n";}
	}


	if ($doFastTree){
		unless ($continue && -e $treeOpts{fastTrOut}){
			runFasttree($treeOpts{inMSA},$treeOpts{fastTrOut},$treeOpts{useAA},$treeOpts{ncore});
		}
		$phyloTree = $treeOpts{fastTrOut};
	}
	if ($doIQTree){
		my $IQtree = "$treeOpts{IQtreeout}";
		unless ($continue && -e "$IQtree.treefile"){
			runQItree(\%treeOpts);
		}
		$phyloTree = "$IQtree.treefile";
	}
	if ($doRAXMLng){
		unless ($continue && -e $treeOpts{RAXNGtreeout}){
			runRaxMLng(\%treeOpts);
		}
		$phyloTree = $treeOpts{RAXNGtreeout};
	}
	if ($doRAXML){
		$treeOpts{inMSA} = "$multF.ph";
		if (!-e $treeOpts{inMSA}){ die "Can't find expected *.ph file: $multF.ph";}
		runRaxML(\%treeOpts);#"$multF.ph",$bootStrap,$outgroup,"$treeD/RXML_allsites$BStag.nwk",$ncore,$continue,!$useAA4tree);
		$phyloTree = $treeOpts{RAXtreeout};#"$treeD/RXML_$siteTag$BStag.nwk";
	}
	if ($doCFML && !$treeOpts{isSubTree}){
		my $outDG = "$outD/clonalFrameML/";
		system "mkdir -p $outDG" unless (-d $outDG);
		$outDG .= "CFML";
		my $CFMLbin = "/g/bork3/home/hildebra/bin/ClonalFrameML/src/./ClonalFrameML";
		my $cmd = "$CFMLbin $phyloTree $multF $outDG\n";
		die $cmd;
	}

	#phyml

	if ($doPhym){
		die "Phym is outdated, check code to reactivate..\n";
		my @thrs;
		my $tcmd = "";
		my $nwkFile = $treeOpts{PhymTree};#"$treeD/phyml${tcnt}_${siteTag}.nwk";
		$tcmd = "$phymlBin --quiet -m GTR --no_memory_check -d nt -f m -v e -o tlr --nclasses 4 -b 2 -a e -i $multF.ph > $nwkFile\n";
		push(@thrs, threads->create(sub{system $tcmd;}));
		for (my $t=0;$t<@thrs;$t++){
			my $state = $thrs[$t]->join();
			if ($state){die "Thread $t exited with state $state\nSomething went wrong with RaxML\n";}
		}
	}
	$treeOpts{nwk} = $phyloTree;
	return (\%treeOpts);
}



sub singleGeneMSAprocess($){
#$multAli,\@MSAs,\@MSA_AA
	my ($multAli) = @_;#,$MSAsar,$MSA_AAar) = @_;
	#my @MSAs = @{$MSAsar};
	#my @MSA_AA=@{$MSA_AAar};
	print "No gene categories given, assumming 1 gene / species in input\n";
	my $tmpInMSA = $aaFna;
	#my $tmpInMSAnt = $fnFna;
	my $tmpOutMSA2 = "$tmpD/outMSA.faa";
	#my $tmpOutMSAsyn = $multAliSyn;#"$tmpD/outMSA.syn.fna";
	#my $tmpOutMSAnonsyn = $multAliNonSyn;
	
	#print "$tmpInMSAnt\n";
	my $numFas;
	if (-e $aaFna){ $numFas = `grep -c '^>' $aaFna`;#just faster to use aa file..
	} else {$numFas = `grep -c '^>' $fnFna`;
	}
	chomp $numFas;
	#die "$fnFna $numFas\n"; 
	my $seqType = "AA";
	my $seqTypeOth = "NT";
	my $inFasta = $aaFna;
	my $inFastaOth = $fnFna;
	my $mainTypeIsAA=1;
	if ($aaFna eq ""){ #always choose AA as default alignment, nt is only fallback
		$inFasta = $fnFna ;		$inFastaOth = $aaFna;
		$seqType = "NT";$seqTypeOth = "AA";	$mainTypeIsAA = 0;
	}
	if ($numFas <= 1){print "Not enough Sequences\n"; exit(0);}
	if ($isAligned){
		system "cp $inFasta $tmpOutMSA2";
	} else{
		#MSA calculation
		$tmpOutMSA2 = MSA($inFasta,$tmpOutMSA2,$ncore,$clustalUse,$continue,$numFas);
	}

	if ($calcDistMat){
		if (-e $tmpInMSA){
			calcDisPos2($tmpInMSA,"$outD/MSA/${seqType}_vsearch_percID.txt",!$mainTypeIsAA,$ncore,$tmpD);
		}
		if (-e $inFastaOth){
			calcDisPos2($inFastaOth,"$outD/MSA/${seqTypeOth}_vsearch_percID.txt",$mainTypeIsAA,$ncore,$tmpD);
		}
	}

	if ($tmpInMSA ne "" && !$useAA4tree){
		convertMultAli2NT($tmpOutMSA2,$fnFna,$multAli);
		$cmd = "$msaFbin -i $multAli  -maskLowID -maskBorderGap -rmGapColsGreater ".$maxGapPerCol." -minGoodPosFrac 0.6\n";
		systemW $cmd;
		($multAliSyn, $multAliNonSyn) = synPosOnly($multAli,$tmpOutMSA2,0,"",$calcSyn,$calcNonSyn);

		#system "rm $tmpInMSA $fnFna $tmpOutMSA2";
		system "rm -f $tmpOutMSA2";
	} else {
		system "mv $tmpOutMSA2 $multAli";
	}
	#$multAli = $tmpOutMSA; $multAliSyn = $tmpOutMSAsyn;
	
	print "OUTPUT:: $multAli\n";
	
	return $multAli;#,\@MSAs,\@MSA_AA);

}

#### fastgear ##

sub FastGear{

	if($doFastGear){
#		open I,"<$cogCats" or die "Can't open cogcats $cogCats\n"; 
		my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");

		while (<$xI>){
			my $cnt3 =0;
			chomp; my @splF = split /\t/;
			@splF = grep !/^NA$/, @splF;#remove NAs
			if (@splF ==0){print "No categories in cat file line $cnt3\n";next;}
			$splF[0] =~ m/^(.*)$smplSep(.*)$/;					
			my @splF2 = ($1,$2);#split /$smplSep/,$spl[0] ;	
			#die "@spl\n";
			push(@geneListF, $splF2[1]);
			$cnt3 ++;
		}
		close $xI;

		my $MsaDF1 = "$outD/MSA";
		#my $MsaDF2 = "$outD/MSA_FG";
		my $MsaDF2 = "$outD_clust/MSA_FG";
		my $MsaDF3 = ""; #added by Falk, debug, TODO!!
		system "mkdir $MsaDF3";		
		system "cp -r $MsaDF1 $MsaDF2";
		

		foreach my $geneF (@geneListF){
			my $outFG = "$outD/fastGear/fastGear_Results/$geneF";
			system "mkdir -p  $outFG" unless(-d "$outFG");
			my $outFileFG = "$outFG/${geneF}_res.mat";
			system "cat $MsaDF2/$geneF.*.fna | sed 's/_.*\$//' > $MsaDF2/$geneF.fna";
			my $FGparFile ="/g/bork3/home/luetge/softs/fastGEARpackageLinux64bit/fG_input_specs.txt";
			runFastgear($geneF, $outFileFG, $MsaDF2, $FGparFile);
		}
		system "rm -fr $outD_clust";
		#die;
	}
		


	### postprocessing fastgear output ##
		
	if($doFastGearSummary){
		my $FGDataD = "$outD/fastGear/";
		my $summaryD = "$FGDataD/fastGear_Summaries";
		system "mkdir -p  $summaryD" unless(-d "$summaryD");
		my $resultD = "$FGDataD/fastGear_Results";
		die "no fastgear results found\n" unless(-d $resultD);
			
		#die "$treeFileFG\n";
		my $allNamesFile ="$summaryD/allNamesFromTop.txt";
		my $treeNamesFile ="$summaryD/subtreeNamesFromTop.txt";
		
		# get a list with all genomes in tree
		if(! -e $allNamesFile |! -e $treeNamesFile){
			my @genomeListFG;		
			my @FNAheader_all = `grep '^>' $multAli`;
			#die "@FNAheader_all\n";
			foreach my $genome (@FNAheader_all){
				my ($genome2) = $genome =~ m/>(.*)?/;
				push (@genomeListFG, $genome2);
				}
			open T1,">$allNamesFile";
			print T1 join("\n",@genomeListFG);
			close T1;	
			#die;
			open T2,">$treeNamesFile";
			print T2 join("\n",@genomeListFG);
			close T2;
			#die "@genomeListFG\n";
		}

		# reorder files -> required for further steps?
		my $reorder_cmd = "$fastgearReorderBin $matlabBin $FGDataD fastGear_ allNamesFromTop.txt both";
		#die "$reorder_cmd\n";
		system $reorder_cmd; 
		
		## collect Recombination statistics #

		my $SRC_cmd = "$fastgearSummaryBin $matlabBin $FGDataD fastGear_";
		system $SRC_cmd; 
		print "fastgear collect recombination statistics finished";
		my $FG_sumOut = "$summaryD/fastGear__recSummaries.txt";
		if(-e $FG_sumOut){system "rm -rf $resultD";}
		#die;
	}

}

sub mergePids($ $ $ $){
	my ($dir,$max,$seqType,$tag) = @_;
	#$outD/MSA/${seqType}_clustalo_percID_$cnt.txt
	#my $seqType = "AA";  
	opendir(DIR, $dir) or die $!;
    my @subfls  = grep { /$seqType.*_percID_.*\.txt/ } readdir(DIR); #_percID_.*\.txt
	close DIR;
	return if (@subfls == 0);
	#die "@subfls\n$dir\n";
	my %bigMat; my %bigCnt;
	for (my $j=0;$j<$max;$j++){
		my $disM = $subfls[$j]; #"$outD/MSA/${seqType}_clustalo_percID_$j.txt";
		$disM =~ m/_(\d+)\.txt$/;
		my $curL = $1;
		die "can't find distance matrix $dir$disM\n" unless (-e "$dir/$disM");
		my $cc=-2;
		my @IDS; my %cL;
		#print "$disM\n";
		open I ,"<$dir/$disM";
		while (my $line = <I>){
			$cc++;
			next if ($cc==-1);
			chomp $line;
			my @spl = split /\s+/,$line;
			my $id = shift @spl;
			push (@IDS,$id);
			$cL{$cc} = \@spl;
		}
		close I;
		system "rm -f $dir/$disM";
		#matrix in mem, now relate to actual dist matrix
		for (my $j=0;$j<@IDS;$j++){
			$IDS[$j] =~ m/([^_]+)_/;
			my $id1 = $1;
			for (my $k=$j;$k<@IDS;$k++){
				my $curPID = $cL{$j}[$k];
				$IDS[$k] =~ m/([^_]+)_/;
				my $id2 = $1;
				#print "$id1 $id2 $curPID\n";;
				$bigMat{$id1}{$id2} += $curPID * $curL;
				$bigCnt{$id1}{$id2} += $curL;
			}
		}
		#die;
	}
	my $oMat = "$dir/${seqType}${tag}_percID.txt";
	open O,">$oMat" or die "Can't open out dis mat $oMat\n";
	my @IDS = sort(keys %bigMat);
	print O "\t".join("\t",@IDS)."\n";
	for (my $j=0;$j<@IDS;$j++){
		my $id1 = $IDS[$j];
		print O $id1;
		for (my $k=0;$k<@IDS;$k++){
			my $id2 = $IDS[$k];
			if (exists($bigMat{$id1}{$id2} )){
				print O "\t".$bigMat{$id1}{$id2} / $bigCnt{$id1}{$id2};
			} elsif (exists($bigMat{$id2}{$id1} ) ) {
				print O "\t".$bigMat{$id2}{$id1} / $bigCnt{$id2}{$id1};
			} else {
				print O "\tNA"
				#die "$id1 $id2\n";
			}
		}
		print O "\n";
	}
	
	close O;
	#die $oMat;
}




#just get positions different between alignments, and their relative position
sub calcDiffDNA($ $){
	my ($MSA,$opID) = @_;
	my $kr = readFasta($MSA);
	my %MS = %{$kr};
	my $isNT = 1;
	my %diffArs;my %perID;
	print "Calculating distance matrix..\n";
	foreach my $k1 (keys %MS){
		$MS{$k1} = uc ($MS{$k1});
	}
	foreach my $k1 (keys %MS){
		my $ss1 = $MS{$k1};
		foreach my $k2 (keys %MS){
			next if ($k2 eq $k1);
			my $ss2 = $MS{$k2};
			my $mask = $ss1 ^ $ss2;
			my $diff=0;
			my$N2=($ss2 =~ tr/[-]//);
			my$N1=($ss1 =~ tr/[-]//);
			if ($isNT){
				$N1+=($ss1 =~ tr/[N]//);$N2+=($ss2 =~ tr/[N]//);
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			} else {
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "N" ||$s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "N" ||$s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			}
			my $nonDiff = ($mask =~ tr/[\0]//);
			$nonDiff -= $N1;
			$perID{$k1}{$k2}= $nonDiff/($diff+$nonDiff)*100;
		}
	}
	open O,">$opID" or die "Cant open out perc ID file $opID\n";
	my @smpls = keys %MS;
	print O "percID\t".join("\t",@smpls)."";
	foreach my $k1 (@smpls){
		print O "\n$k1\t";
		foreach my $k2 (@smpls){
			if ($k1 eq $k2){print O "\t100";
			} else {
				print O "\t".$perID{$k1}{$k2};
			}
		}
	}
	close O;
	my @NTdiffs = sort {$a <=> $b} (keys %diffArs);
	#die "Cant open out perc ID file $opID\n";
	my $MSAredF = $MSA; $MSAredF =~ s/\.[^\.]+$//;
	my $NSAposF = $MSAredF; $MSAredF.=".reduced.fna";
	$NSAposF .= ".reduced.pos";
	open O2,">$NSAposF" or die "Can't open reduced MSA position file $NSAposF\n";
	foreach my $i (@NTdiffs){
		print O2 "$i\n";
	} close O2;

	open O,">$MSAredF" or die "Can't open reduced MSA file $MSAredF\n";
	foreach my $k1 (@smpls){
		my $seq1 = $MS{$k1};
		my $seq = "";my $pos="";
		foreach my $i (@NTdiffs){
			$seq.=substr($seq1,$i,1);
			$pos .="$i\n";
		}
		print O ">$k1\n$seq\n";

	}
	close O;
	print "Dont calculating Distance matrix\n";
	#die "done\n";
}

sub calcDisPos($ $ $){
	my ($MSA,$opID, $isNT) = @_;
	#$cmd = $clustaloBin." -i $MSA -o $MSA.tmp --outfmt=fasta --percent-id --use-kimura --distmat-out $opID --threads=$ncore --force --full\n";
	#$cmd .= "rm -f $MSA.tmp\n";
	
	#too slow
	my $kr = readFasta($MSA);
	my %MS = %{$kr};
	my %diffArs;my %perID;
	print "Calculating distance matrix..\n";
	foreach my $k1 (keys %MS){
		$MS{$k1} = uc ($MS{$k1});
	}
	foreach my $k1 (keys %MS){
		my $ss1 = $MS{$k1};
		foreach my $k2 (keys %MS){
			next if ($k2 eq $k1);
			my $ss2 = $MS{$k2};
			my $mask = $ss1 ^ $ss2;
			my $diff=0;
			my$N2=($ss2 =~ tr/[-]//);
			my$N1=($ss1 =~ tr/[-]//);
			if ($isNT){
				$N1+=($ss1 =~ tr/[N]//);$N2+=($ss2 =~ tr/[N]//);
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			} else {
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "X" ||$s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "X" ||$s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			}
			my $nonDiff = ($mask =~ tr/[\0]//);
			$nonDiff -= $N1;
			$perID{$k1}{$k2}= $nonDiff/($diff+$nonDiff)*100;
		}
	}
	open O,">$opID" or die "Cant open out perc ID file $opID\n";
	my @smpls = keys %MS;
	print O "percID\t".join("\t",@smpls)."";
	foreach my $k1 (@smpls){
		print O "\n$k1\t";
		foreach my $k2 (@smpls){
			if ($k1 eq $k2){print O "\t100";
			} else {
				print O "\t".$perID{$k1}{$k2};
			}
		}
	}
	close O;
	my @NTdiffs = sort {$a <=> $b} (keys %diffArs);
	#die "Cant open out perc ID file $opID\n";
	my $MSAredF = $MSA; $MSAredF =~ s/\.[^\.]+$//;$MSAredF.=".reduced.fna";
	open O,">$MSAredF" or die "Can't open reduced MSA file $MSAredF\n";
	foreach my $k1 (@smpls){
		my $seq1 = $MS{$k1};
		my $seq = "";
		foreach my $i (@NTdiffs){
			$seq.=substr($seq1,$i,1);
		}
		print O ">$k1\n$seq\n";
	}
	close O;
	print "Dont calculating Distance matrix\n";
	#die "done\n";
}

sub mergeMSAs($ $ $ $){
	my ($MSAsAr,$samplesHr,$multAliF,$del,$isAA) = @_;
	my @MSAs = @{$MSAsAr}; my %samples = %{$samplesHr};
	my @smps = keys %samples;
	if (@smps == 0){#no cats file
		push(@MSAs ,$multAliF);
		return;
	}
	my %bigMSAFAAnxs;my %bigMSAFAA;foreach my $sm (@smps){$bigMSAFAA{$sm} ="";$bigMSAFAAnxs{$sm}="";}
	my @lengthsParts;
	foreach my $MSAf (@MSAs){
		#print $MSAf."\n"; 
		my $hit =0; my $miss =0; my $keyI=0;
		my $hr = readFasta($MSAf,1); my %MFAA = %{$hr};
		system "rm -f $MSAf" if ($del);
		my @Mkeys = keys %MFAA;
		next if (@Mkeys == 0);
		#die "$Mkeys[0]\n";
		my $smplSep1 =$smplSep; $smplSep1 =~ s/\\//g;
#		my @spl2 = split /$smplSep/,$Mkeys[0];
		while ($keyI < @Mkeys && $Mkeys[$keyI] !~ m/^(.*)$smplSep(.*)$/){
			print STDERR "Can't recognize key $Mkeys[$keyI]\n";
			$keyI++;
		}
		$Mkeys[$keyI] =~ m/^(.*)$smplSep(.*)$/;my @spl2 =($1,$2);
		my $gcat = $spl2[1];
		my $len = length( $MFAA{$Mkeys[$keyI]} );
		if ($len == 0 || !defined($MFAA{$Mkeys[$keyI]})){print STDERR "0 length sequence discovered\n";next ;}
		push(@lengthsParts,$len);
		foreach my $sm (keys %samples){
			my $curK = $sm.$smplSep1.$gcat; #print $curK. " ";
			if ( exists($MFAA{$curK}) && defined($MFAA{$curK})  ) {
				my $seq = $MFAA{$curK}; 
				die "Sequence lengths within MSA unequal: $len != ".length($seq)."\n" if (length($seq) != $len);
				$hit++;
				$bigMSAFAA{$sm} .= $seq;
				$seq =~ s/^(-+)/"?" x length($1)/e;
				$seq =~ s/(-+)$/"?" x length($1)/e;
				#die $seq;
				$bigMSAFAAnxs{$sm} .= $seq;
			} else {
				$bigMSAFAA{$sm} .= "-"x$len; $miss++;#print "nooooooo ";
				$bigMSAFAAnxs{$sm} .= "?"x$len; $miss++;
			}
		}
		
		#die "$hit - $miss\n";
	}
	#filter part - count "-" in each seq
	my $factor = 1; $factor = 3 if ($isAA);
	my @ksMSAFAA = keys %bigMSAFAA;
	my $iniSeqNum = @ksMSAFAA; my $remSeqNum = 0;
	my %charCnts; my $maxNtCnt=0;
	#my @usedPos; #is now implemented in C++ program..
	#simply count gaps and N's
	foreach my $kk (@ksMSAFAA){
		#my $strCpy = $bigMSAFAA{$kk};
		my $num1 = 0;
		if ($isAA){
			$num1 = $bigMSAFAA{$kk} =~ tr/[\-Xx]//;
		} else {
			$num1 = $bigMSAFAA{$kk} =~ tr/[\-Nn]//;
		}
		
		$charCnts{$kk} = (length($bigMSAFAA{$kk})-$num1);
		#print "$kk GGGG  $charCnts{$kk} $num1\n";
		if ( $charCnts{$kk} > $maxNtCnt){
			$maxNtCnt = $charCnts{$kk};
		}
		next; #overlap implemented in C++
	}
	if ($maxNtCnt == 0){ #something really wrong
		print "No useable MSA positions.. not good.. aborting\n";
		die;
	}
	
	my $qtl90NTcnts = quantile(0.9,values(%charCnts));
	my $qtl50NTcnts = quantile(0.5,values(%charCnts));
	my $qtl25NTcnts = quantile(0.25,values(%charCnts));
	
	
	#final check on MSA's that enough data is present
	foreach my $kk (@ksMSAFAA){
		my $num1 = $charCnts{$kk};
		#print "$num1\n";
		if ( $maxNtCnt == 0 ||  ($num1 < ($qtl90NTcnts *$ntFrac) && $num1 < $qtl25NTcnts) || ($num1 < ($ntCntTotal/$factor) ) ){
			print "($num1 / $qtl90NTcnts ) < $ntFrac) && $num1 < $qtl25NTcnts || ($num1 < ($ntCntTotal/$factor) )\n";
			delete $bigMSAFAA{$kk}; delete $bigMSAFAAnxs{$kk}; $remSeqNum++; 
			if ($maxNtCnt == 0){
				print "$kk $num1  infinite $ntFrac \n";
			} else {
				print "$kk $num1  " . int($num1*1000 / ($maxNtCnt) )/1000 ." $ntFrac \n";
			}
		}
		#print "$num1  $kk \n";#$bigMSAFAA{$kk}\n\n"; last;
	}
	open O,">$multAliF" or die "Can't open MSA outfile $multAliF\n";
	open O2,">$multAliF.nxs" or die "Can't open MSA nexus outfile $multAliF.nxs\n";
	my @allKs = keys %bigMSAFAA;
	if (@allKs == 0){die "no genes for nexus output format.\nAborting\n";}
	print O2 "#NEXUS\nBegin data;\nDimensions ntax=".scalar(@allKs)." nchar=".length($bigMSAFAAnxs{$allKs[0]}).";\nFormat datatype=dna missing=? gap=-;\nMatrix\n";
	my $scnt=0;
	foreach my $kk (@allKs){
		print O ">$kk\n"; my $s1 = $bigMSAFAA{$kk};
		print O2 "\n$kk\t"; my $s2 = $bigMSAFAAnxs{$kk};
		#if (length($s1) !=  @usedPos){die "s1 and usedPos are not same length: ".length($s1)." : ".@usedPos."\n";}
		#if ($minOverlapMSA <=0){
		print O "$s1\n"; print O2 "$s2\n";
		#} else {
		#	my $usedP=0;
		#	for (my $i=0; $i<length($s1);$i++){
		#		if ($usedPos[$i] >= $minOverlapMSA){
		#			print O substr($s1,$i,1);
		#			print O2 substr($s2,$i,1);
		#			$usedP++;
		#		}
		#	}
		#	print O "\n"; print O2 "\n";
		#	print "Used $usedP positions of ". length($s1)  . " total positions.\n" if ($scnt==0);
		#	$scnt++;			
		#}
	}
	print O2 "\n;\nend;";

	close O;close O2;
	#die "$multAliF\n";
	
	#prepare partition file to record segment lengths..
	my $partiFile = $multAliF.$partiExt;
	open O,">$partiFile" or die "Can't open output partioning file $partiFile\n";
	my $lastP=0;
	my $TypeTag = "DNA";
	$TypeTag = "LG" if ($isAA); #this is the model to be used...
	for (my $i=0;$i<@lengthsParts;$i++){
		#DNA, part1 = 1-100
		print O "$TypeTag, part".($i+1) ." = ". ($lastP+1) ."-". ($lengthsParts[$i]+$lastP) ."\n";
		$lastP+=$lengthsParts[$i];
	}
	close O;
	
	print "Removed $remSeqNum of $iniSeqNum sequences\n";
}


sub convertMultAli2NT($ $ $){
	my ($inMSA,$NTs,$outMSA) = @_;
	my $tmpMSA=0;
	if ($inMSA eq $outMSA){$outMSA .= ".tmp"; $tmpMSA=1;}
	my $cmd = "";
	#"$pal2nal $inMSA $NTs -output fasta -nostderr -codontable 11 > $outMSA\n";
	
	#$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon\n";
	$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -ignorestopcodon  -gt 0.1 -cons 60 2>/dev/null\n";
	#die "$cmd\n$inMSA,$NTs,$outMSA\n";
	#my $hr1= readFasta($inMSA);
	#my %MSA = %{$hr1};
	#$hr1= readFasta($NTs);
	#my %NTs = %{$hr1};
	if ($tmpMSA){$cmd .= "rm -f $inMSA;mv $outMSA $inMSA;\n";}
	#print $cmd;
	#die "$cmd\n";
	die "Can't execute $cmd\n" if (system $cmd);
}

sub synPosOnlyAA($ $){#only leaves "constant" AA positions in MSA file.. 
#stupid, don't know if pal2nal can handle this.. prob not
	my ($inMSA,$outMSA) = @_;
	#print "Syn";
	my $hr = readFasta($inMSA,1); my %FNA = %{$hr};
	my @aSeq = keys %FNA;
	my $len = length ($FNA{$aSeq[0]});
	for (my $i=0; $i< $len; $i+=3){
		my $cod = substr $FNA{$aSeq[0]},$i,3;
		my $iniAA = "A";
		for (my $j=1;$j<@aSeq;$j++){
		}
	}
	#print " only\n";

}

sub synPosOnly{#now finished, version is cleaner
	my ($inMSA,$inAAMSA, $ffold, $outgroup, $doSyn, $doNSyn) = @_;
	
	my $outMSA = $inMSA; 	my $outMSAns = $inMSA;
	$outMSAns =~ s/\.fna/\.nonsyn\.fna/;
	$outMSA =~ s/\.fna/\.syn\.fna/;

	#print "Syn NT";
	my %convertor = (
    'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',    # Serine
    'TTC' => 'F', 'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L', 'TTG' => 'L',    # Leucine
    'TAC' => 'Y',  'TAT' => 'Y',    # Tyrosine
    'TAA' => '*', 'TAG' => '*', 'TGA' => '*',    # Stop
    'TGC' => 'C', 'TGT' => 'C',    # Cysteine   
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',    # Leucine
    'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',    # Proline
    'CAC' => 'H', 'CAT' => 'H',    # Histidine
    'CAA' => 'Q', 'CAG' => 'Q',    # Glutamine
    'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',    # Arginine
    'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',    # Threonine
    'AAC' => 'N','AAT' => 'N',    # Asparagine
    'AAA' => 'K', 'AAG' => 'K',    # Lysine
    'AGC' => 'S', 'AGT' => 'S',    # Serine
    'AGA' => 'R','AGG' => 'R',    # Arginine
    'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',    # Valine
    'GCA' => 'A','GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',    # Alanine
    'GAC' => 'D', 'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E', 'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G','GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',    # Glycine
    );
	my %ffd;
	if ($ffold){ #calc 4fold deg codons in advance to real data
		foreach my $k (keys %convertor){
			my $subk = $k; my $iniAA = $convertor{$subk} ;
			next if (exists($ffd{$iniAA}));
			my $cnt=0;
			foreach my $sNT ( ("A","T","G","C") ){
				
				substr ($subk,2,1) = $sNT;
				#print $subk ." " ;
				$cnt++ if ($convertor{$subk} eq $iniAA);
				
			}
#			if( $cnt ==4){ $ffd{$k} = 4;
#			} else {$ffd{$k} = 1;}
			if( $cnt ==4){ $ffd{$iniAA} = 4;
			} else {$ffd{$iniAA} = 1;}
			#die"\n$ffd{$iniAA}\n";
		}
	}

	#assumes correct 3 frame for all sequences in inMSA
	my $hr = readFasta($inMSA,1); my %FNA = %{$hr};
	#my %FAA;
	#if (0  || !$ffold){
	#	$hr = readFasta($inAAMSA); %FAA = %{$hr};
	#}
	#print "$inMSA\n$inAAMSA\n$outMSA\n";
	my @aSeq = keys %FNA; 
	my %outFNA;#syn
	my %outFNAns;#non syn
	for (my $j=0;$j<@aSeq;$j++){$outFNA{$aSeq[$j]}="";}
	my $len = length ($FNA{$aSeq[0]});
	my $nsyn=0;my $syn=0;
	for (my $i=0; $i< $len; $i+=3){ #goes over every position
		my $j =0;
		my $iniAA = "-";
		my $iniCodon ;
		while (1){ #check for first informative position
			$iniCodon = substr $FNA{$aSeq[$j]},$i,3;
			if ($iniCodon =~ m/---/ || $iniCodon =~ m/N/i){$j++; last if ($j >= @aSeq); next;}
			die "error: $iniCodon\n" if ($iniCodon =~ m/-/); #should not happen
			die "codon doesn't exist $iniCodon \n" unless (exists($convertor{$iniCodon}));
			$iniAA = $convertor{$iniCodon};#substr $FAA{$aSeq[0]},$i,1; 
			last;
		}
		#die "$iniAA\n";
		my $isSame = 1;my $ntSame = 1;
		next unless (!$ffold || $ffd{$iniAA} == 4);
	#print $i." $iniAA ";
		for (;$j<@aSeq;$j++){
			if ($aSeq[$j] eq $outgroup){next;}#print "HIT   $aSeq[$j] eq $outgroup";
			my $newCodon = substr $FNA{$aSeq[$j]},$i,3;
			my $newAA = "-";
			if ($newCodon !~ m/-/ && $newCodon =~ m/[ACTG]{3}/i){
				die "Unkown AA $newCodon\n" unless (exists $convertor{$newCodon} );
				$newAA = $convertor{$newCodon} ; # substr $FAA{$aSeq[$j]},$i,1;
			} elsif ($newCodon =~ m/---/ || $newCodon =~ m/[NWYRSKMDVHB]/i){
			} else {
				die "newCodon wrong $newCodon\n" ;
			}
			if ($iniAA ne $newAA && $newAA ne "-"){
				$isSame =0; $ntSame =0; last;
			}
			if ($iniCodon ne $newCodon){
				$ntSame=0;
			}
		}
		for (my $j=0;$j<@aSeq;$j++){
			my $curCod = substr $FNA{$aSeq[$j]},$i,3;
			if ($ntSame){#add nts to file
				$outFNA{$aSeq[$j]} .= $curCod;
				$outFNAns{$aSeq[$j]} .= $curCod;
			} elsif ($isSame){#add nts to file
				if ($ffold){
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},($i)+2,1;
				} else {
					$outFNA{$aSeq[$j]} .= $curCod;
				}
				#print substr $FNA{$aSeq[$j]},$i*3,3 . " ";
				$syn++;
			} else {
				$outFNAns{$aSeq[$j]} .= $curCod;
				$nsyn++;
			}
		}
	}
	#die $inMSA."\n";
	if ($doSyn){
		if ($syn ==0){
			$outMSA = "";
		} else {
			open O ,">$outMSA" or die "Can't open outMSA $outMSA\n";
			for (my $j=0;$j<@aSeq;$j++){
				print O ">$aSeq[$j]\n$outFNA{$aSeq[$j]}\n";
			}
			close O;
		}
	}
	if ($doNSyn){
		if ($nsyn ==0){
			$outMSAns = "";
		} else {
			open O ,">$outMSAns" or die "Can't open outMSA $outMSAns\n";
			for (my $j=0;$j<@aSeq;$j++){
				print O ">$aSeq[$j]\n$outFNAns{$aSeq[$j]}\n";
			}
			close O;
		}
	}
	$aSeq[0] =~ m/^.*$smplSep(.*)$/;
	#die "$outMSA\n";
	print "$1 ($syn / $nsyn) ".@aSeq." seqs \n";
	#print " only\n";
	#print "\n";
	return ($outMSA,$outMSAns);
}

sub codeml{
	my ($MSAfile2,$codemlOutDTmp,$gene,$nwkFile_gene2,$repeatCounts) = @_;
	my @omegaStart = @omegas;			
	my $codemlOutDFile = "$codemlOutD/${gene}_run2";
	system "mkdir -p  $codemlOutDFile" unless(-d $codemlOutDFile);

	chdir $codemlOutDTmp;
	my $modelName;
	for (my $mod=0; $mod < scalar(@model); $mod++){		
		$modelName = $model[$mod];
		my @repSel;
		for (my $rep=1; $rep <= $repeatCounts; $rep++) {

			open M0,">$codemlOutDTmp/${gene}_${rep}_${modelName}.c" or die "Can't open control file for codeml: $gene\n";
	
			print M0 "seqfile = $MSAfile2\n";
			print M0 "verbose = 2\n";
			print M0 "treefile = $nwkFile_gene2\n";
			print M0 "outfile = $codemlOutDTmp/codemlOut_${gene}_${rep}_${modelName}.txt\n";
			print M0 "aaDist = 0\n";
			print M0 "fix_blength = 0\n";
			print M0 "runmode = 0\n";
			print M0 "seqtype = 1\n";
			print M0 "CodonFreq = 2\n";
			print M0 "clock = 0\n";
			print M0 "model = 0\n";
			print M0 "NSsites = $modelName\n";
			print M0 "fix_omega = 0\n";
			print M0 "omega = $omegaStart[($rep-1)]\n";
			print M0 "cleandata = 0\n";
			print M0 "getSE = 0\n";	
			print M0 "icode = 0\n";
			print M0 "fix_kappa = 0\n";
			print M0 "kappa = 2\n";
			print M0 "Mgene = 0\n";
			print M0 "ncatG = 8\n";
			print M0 "RateAncestor = 0\n";
			print M0 "Small_Diff = 1e-6\n";
			print M0 "noisy = 0\n";
			

			close M0;

			## run codeml  
			$cmd = "$pamlBin $codemlOutDTmp/${gene}_${rep}_${modelName}.c\n";			
			die "$cmd\n";
			systemW $cmd; 
			print "finished codeml Model $modelName run $rep on $gene\n";
			#die;
			
			#get lnL and push to array
			open(CM, "<$codemlOutDTmp/codemlOut_${gene}_${rep}_${modelName}.txt" ) or die "could not find $!";
			while (my $line = <CM>) {
					if ($line =~ /^lnL*/) {
							$line =~ m/^.*\):\s*([-+]?[0-9]*\.?[0-9]+)\s.*$/;
					push @repSel, $1;
					last;
					}				
			}
			close CM;
			#system "rm $codemlOutDTmp/rst1";

		} 
		
		my $idxMax = 0;
			$repSel[$idxMax] > $repSel[$_] or $idxMax = $_ for 1 .. $#repSel; 
		my $repSelected = $idxMax+1;
		#die "@repSel\n$repSelected\n";
		system "cp $codemlOutDTmp/codemlOut_${gene}_${repSelected}_${modelName}.txt $codemlOutDFile/out_M$modelName.txt";
	}	
}


#starts my R script to calc popgenStats and filter out bad genes (in multi gene approach)
sub pogenStatsFilter{
	return [] unless ($doTheta);
	system "mkdir -p $MSAsubsD" unless (-d $MSAsubsD);
	system "mkdir -p $codemlOutD" unless (-d $codemlOutD);
	system "rm -f $codemlOutD/PopStats.*";
	my $condaAct = getProgPaths("CONDA");
	my $cmd = "";
	#$cmd .= "$condaAct\nconda activate r_env\n"; #this was a workaround since the packages didn't run correctly on R cluster and different env was needed..
	$cmd .= "$RpogenS $outD $mapF $codemlOutD $MSAsubsD $subsetPopgenStats\n";
	#$cmd .= "conda deactivate\n";
	#die "$cmd\n";
	my $ret = `$cmd`;
	$ret =~ m/Outliers: (.*);;/;
	my @spl = split /,/,$ret;
	return \@spl;
} 

sub hyphy{
	my ($MSAfile2,$codemlOutDTmp,$gene,$nwkFile_gene2,$log) = @_;
	my $hyphyBin=getProgPaths("hyphy");
	my $cmd = "";#"source activate hyphy\n";
	$cmd .= "$hyphyBin CPU=$ncore fubar --alignment $MSAfile2 --tree $nwkFile_gene2 > $log\n";
	$cmd .= "gzip -c $MSAfile2.FUBAR.json > $log.json.gz\n";
	$cmd .= "rm -f $MSAfile2.FUBAR.*\n";
	print $log."\n";
	systemW $cmd ;#if (!-e $log);
	#die $cmd ;#if (!-e $log);
	return $log;
}

sub pruneTree($ $ $){
	my ($nwkFile,$aR,$nwkFile_gene) = @_;
	my @genomeList = @{$aR};
	my $cmd_prune = "$eteBin mod -t $nwkFile --prune @genomeList --unroot -o $nwkFile_gene";
	systemW $cmd_prune . " > $nwkFile_gene";
}

sub fubarXML($){
	my $inp = $_[0];
	return "" if (!-e $inp || (-e "$inp.json.gz" && !$reparseHyphyJson));
	
	if (0){#python
		return "";
	}
	my $rDNm; my $rDSm;
	my $txt = `zcat $inp.json.gz`;
	my $str=decode_json($txt);
	my %JS = %{$str};
	my @XX = @{$JS{MLE}{"content"}{"0"}};
	#print @XX."  BB  @XX\n";
	my $negSel=0; my $posSel=0; 
	my @ds = (); my @dn = ();
	for (my $i=0;$i<@XX;$i++){
		my @YY = @{$XX[$i]};
		push @ds,$YY[0];
		push @dn,$YY[1];
		next if (@YY<4 || !defined($YY[4]));
		$negSel++ if ($YY[3]>0.9);
		$posSel++ if ($YY[4]>0.9);
	}
	my $dnX=$2;my $dsX = $1;
	my $rDNmed = medianArray(@dn);my $rDSmed = medianArray(@ds);
	$rDNm = meanArray(\@dn); $rDSm = meanArray(\@ds);
#	$txt = `cat $inp`;
#	$txt =~ m/\* synonymous rate =  ([\d\.]+)\n.*\* non-synonymous rate =  ([\d\.]+)/;
	#print "hy report dnds: $dnX $dsX; $rDNm  $rDSm; $rDNmed $rDSmed\n"; #this is actually probably wrong
	return $JS{input}{"number of sequences"} ."\t". $JS{input}{"number of sites"} . "\t$negSel\t$posSel\t$rDNm\t$rDSm";
}

sub coreHyPhy{
	my ($MSADir,$gene,$xtra,$nwkFile,$codemlOutDTmp,$logF) = @_;
	return if (-e $logF && -e "$logF.json.gz");
	my $runCodeML = 0;
	my @genomeList;
	opendir(DIR, $MSADir);
	my @MSAfile = grep(/$gene.*$xtra\.fna/,readdir(DIR));
	closedir(DIR);
	if (@MSAfile == 0){
		#die "$gene.*$xtra\.fna";
		return;
	}
	#die "@MSAfile\n";
	
	my $MSAfile2 = "$MSADir/$MSAfile[0]";
	my $MSAfile3 = "$codemlOutDTmp/tmp.$MSAfile[0]";
	#get entries in nwk to match these up as well..
	my $hr = getTreeLeafs($nwkFile);
	my %nwLfs= %{$hr};
	
	#die "@nwLfs\n";
	my $cntMissTree=0;
	$hr = readFasta($MSAfile2,1,"input MSA for selection analysis");
	my %FNA = %{$hr};
	#		print "$MSAfile2\n";

	open O,">$MSAfile3" or die "can't open MSA out $MSAfile3\n";
	foreach my $genome (keys %FNA){
		#print "$genome\n";
		my ($genome2) = $genome =~ m/(.*)?$smplSep/;
		next if ($genome2 =~ m/$outgroup/);
		unless (exists($nwLfs{$genome2})){$cntMissTree++; next;}
		my $seq = $FNA{$genome};
		$seq =~ s/T[AG][GA]$//; #remove stop codon
		my $x=0;
		foreach my $sto (("TGA","TAG","TAA")){
			while ( 1 ){
				$x = index($seq,$sto,$x);
				last if ($x < 0);
				if ($x % 3 ==0){
					print "HIT\n";
					substr($seq,$x,3) = "NNN";
				}
				$x++;
			}
		}
		print O ">$genome2\n$seq\n";
		push (@genomeList, $genome2);
	} 
	close O;
	if ($cntMissTree>0){print "Removed from MSA $cntMissTree sequences due to not being present in tree\n";}
	next if(scalar @genomeList <3);
			
	my $nwkFile_gene = "$codemlOutDTmp/subtree_$gene.nwk";

	####################### Prepare tree ################################## 
	pruneTree($nwkFile,\@genomeList,$nwkFile_gene);
	#$thrs[$thrCnt]->join();
	if ($runCodeML){
		codeml($MSAfile3,$codemlOutDTmp,$gene,$nwkFile_gene,$repeatCounts) #$thrs[$thrCnt] = threads->create( );
	} else {
		hyphy($MSAfile3,$codemlOutDTmp,$gene,$nwkFile_gene,$logF);
	}
	system "rm -f $MSAfile3";
	system "rm -f $nwkFile_gene";
}

sub selecAnalysis($ $ $ $ $){
	my ($geneListRef, $nwkFile, $codemlOutD, $codemlOutDTmp) = @_;
	my @geneListFin = @{$geneListRef};
	#system "source activate hyphy" unless ($runCodeML);
	my $jsonExrScr = getProgPaths("fubarJson_scr");
	my $stdJSONheader = "Gene\tNseqs\tNsites\tnegSel\tposSel\tdn\tds\tdn2\tds2\n";
	#my @thrs;my $thrCnt=0;
	#for ($thrCnt=0;$thrCnt<$ncore;$thrCnt++){		$thrs[$thrCnt] = threads->create(sub{my $x=0;});	}
	#$thrCnt=0;
	#first go through $MsaD
	system "rm -f $codemlOutD/hyphy.fubar*" if ($reparseHyphyJson);
	my $logF1 =  "$codemlOutD/hyphy.fubar.txt";
	if (!-e $logF1 ){
		my %logs;
		foreach my $gene (@geneListFin){
			my $logF =  "$codemlOutD/$gene.hyphy.fubar.log";
			$logs{$gene} = "$logF";
			coreHyPhy($MsaD,$gene,"",$nwkFile,$codemlOutDTmp,$logF);
		}
		my $sumTxt=$stdJSONheader;
		print "reading jsons..";
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");
			$sumTxt .= "$gene\t$summary";
			#print "$sumTxt\n";
		}
		print "Done reading\n";
		open O,">$logF1";		print O $sumTxt;		close O;
	}
	#add by hand the unique seqs subset..
	$logF1 =  "$codemlOutD/hyphy.fubar.unID.txt";
	if (!-e $logF1 ){
		my %logs;
		foreach my $gene (@geneListFin){
			my $logF =  "$codemlOutD/$gene.hyphy.fubar.unID.log";
			$logs{$gene} = $logF;
			coreHyPhy($MSAsubsD,$gene,"\\.uInd",$nwkFile,$codemlOutDTmp,$logF);
		}
		my $sumTxt=$stdJSONheader;
		print "reading jsons unID..";
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");

			$sumTxt .= "$gene\t$summary";
			#system "rm -f $logs{$gene}*";
		}
		print "Done reading\n";
		open O,">$logF1";		print O $sumTxt;		close O;
	}
	#now go through subsets..
	
	foreach my $subs (split /,/,$subsetPopgenStats){
		my %logs;
		#print "$subs\n";
		my $logF1 =  "$codemlOutD/hyphy.fubar.s$subs.txt";
		next if (-e $logF1 );
		foreach my $gene (@geneListFin){
			my $logF =  "$codemlOutD/$gene.hyphy.s$subs.fubar.log";
			$logs{$gene} = $logF;
			#COG0008.0.uInd.s20.fna
			coreHyPhy($MSAsubsD,$gene,"uInd\\.s$subs",$nwkFile,$codemlOutDTmp,$logF);
		}
		my $sumTxt=$stdJSONheader;
		print "reading jsons s$subs..";
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");next if ($summary eq "");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");
			$sumTxt .= "$gene\t$summary";
			#system "rm -f $logs{$gene}*";
		}
		print "Done reading\n";
		open O,">$logF1";		print O $sumTxt;		close O;
	}
	system "rm -fr $codemlOutDTmp" if (-e $codemlOutDTmp);
	
}


#\@geneList,$MsaD,$codemlOutD
sub WattTheta{
	#hyphy /path/to/WattetrsonTheta.txt --alignment /path/to/alignment
	my ($geneListRef, $MSADir, $codemlOutD) = @_;
	my @geneListFin = @{$geneListRef};
	my $runCodeML = 0;
	#system "source activate hyphy" unless ($runCodeML);
	my $hyphyBin=getProgPaths("hyphy");
	#die "XX\n";
	my %logs ;
	my $logF1 =  "$codemlOutD/hyphy.Theta.log";
	return if (-e $logF1);
	my $otxt = "Gene\tNseqs\tNsites\tSegSites\tWattTheta\n";
	foreach my $gene (@geneListFin){
		my $cnt = 0;
		my $logF =  "$codemlOutD/$gene.hyphy.Theta.log";
#		print "$logF\n";
		$logs{$gene} = $logF;
		#next if (-e $logF);
		my @genomeList;
				
		opendir(DIR, $MSADir);
		my @MSAfile = grep(/$gene.*\.fna/,readdir(DIR));
		closedir(DIR);
		#print "@MSAfile\t$gene\t$MSADir\n";
		next if (@MSAfile == 0);
		my $MSAfile2 = "$MSADir/$MSAfile[0]";
		my $hyphyBin=getProgPaths("hyphy");
		my $cmd = "";#"source activate hyphy\n";
		$cmd .= "$hyphyBin CPU=$ncore /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/phylo/WattetrsonTheta.hyphy --alignment $MSAfile2 ";#> $logF\n";
		my $txt = `$cmd`;
#		print $txt."\n";
		$txt =~ m/Sequences          = (\d+)\nSites              = (\d+)\nSegregating Sites  = (\d+)\n.*Watterson.s theta  = ([\d\.]+)/;
		$otxt.="$gene\t$1\t$2\t$3\t$4\n";
	}
	#foreach my $g (keys %logs){	my $txt = `cat $logs{$g}`;	}
	open O,">$logF1" or die "Can't open log out $logF1\n";
	print O $otxt;
	close O;
}
sub fillGeneList{
	my ($cogCats) = @_;
#	open my $xI,"<$cogCats" or die "Can't open cogcats $cogCats\n";
	my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");
	chomp(my @linesCats = <$xI>);
	close $xI;
	foreach (@linesCats){
		chomp; my @spl = split /\t/;
		@spl = grep !/^NA$/, @spl;#remove NAs
		if ($spl[0] =~ m/^#/){shift @spl;}
		$spl[0] =~ m/^(.*)$smplSep(.*)$/;
		my @spl2 = ($1,$2);#split /$smplSep/,$spl[0] ;	
		push(@geneList, $spl2[1]);				
	}
}



### Fastgear -> test for recombination 
sub runFastgear($ $ $ $){
	my ($geneFG, $outFile, $inD, $parFile) = @_;
	$cmd = "$fastgearBin $matlabBin $inD/$geneFG.fna $outFile $parFile";
	#die "$cmd\n";
	system $cmd; 
	print "fastgear on $geneFG finished";
	#die;
}




sub prepGenoDirs($){
	my ($genoindir) = @_;
	# -genoInD '/hpc-home/hildebra/grp/data/DB/Genomes/Ecoli_ref/*.fna'
	if ($genoindir eq ""){return;}
	$aaFna="$outD/autoFAA.faa";$fnFna="$outD/autoFNA.fna";$cogCats="$outD/autoCAT.cat";
	my $SaSe = "_"; my %catT;
	$smplSep = "_";
	
	#return(); #DEBUG
	if (-d $genoindir){
		if ($wildcardflag ne ""){
			$genoindir .= $wildcardflag;#"/*\.fna";
		} else {
			print "Warning: please give wildcard for genoIndir.\nAssumming *.fna now..\n";
			$genoindir .= "/*\.fna";
		}
	}
	#print "$genoindir\n\n";
	my @sfiles = glob($genoindir);
	if (scalar(@sfiles) == 0){die "$genoindir contains no files!\n";}
	
	#die "@sfiles\n";
	print "External Genomes: $genoindir\nProcessing ".scalar(@sfiles)." genomes\n";
	my %FNAfmg; my %FAAfmg;my %MGSFMG;
	
	my $cnt=0;
	foreach my $tarG(@sfiles){
		next if ($tarG  =~ /\*/ || $tarG  =~ /\.genes\./);
	#my $tarG = "$tarDir/$genoN";
		my $tag = $tarG;
		$tag =~ s/.*\///;
		$tag =~ s/\.[^\.]*$//;
		print "$tarG\n";
		my ($genes,$prots) = getGenoGenes($tarG,0,$ncore);
		my $FMGdir = getFMG("",$prots,$genes,$ncore);
		my ($hrN,$hrA,$hrC) = readFMGdir( "$FMGdir",$tag ,".");
		$cnt++;
		my %FAA=%{$hrA};
		my %FNA=%{$hrN};
		my %COGcat=%{$hrC};
		#add to categories..
		foreach my $k (keys %FAA){
			$FAAfmg{$k}=$FAA{$k};
			$FNAfmg{$k}=$FNA{$k};
			die "unkown key $k \n" unless (exists($COGcat{$k}));
			$MGSFMG{$tag}{$COGcat{$k}} = $k;
		}
		#print " $tag ";
		#die;
		#last if ($cnt>5);#DEBUG
	}
	
	
	open OA,">$aaFna"  or die "Can't open faa out file $aaFna\n"; 
	open ON,">$fnFna"  or die "Can't open faa out file $fnFna\n"; 
	foreach my $mg (keys %MGSFMG){
		foreach my $cog (keys %{$MGSFMG{$mg}}){
			my $ng = "$mg$SaSe$cog";
	#		print ON ">$ng\n$FNAfmg{$MGSFMG{$mg}{$cog}}\n";
			die "$MGSFMG{$mg}{$cog}\n" unless (exists( $FAAfmg{$MGSFMG{$mg}{$cog}} ));
			print OA ">$ng\n$FAAfmg{$MGSFMG{$mg}{$cog}}\n";
			print ON ">$ng\n$FNAfmg{$MGSFMG{$mg}{$cog}}\n";
			push(@{$catT{$cog}},"$ng");
		}
	}
	close OA; close ON;

	open OC,">$cogCats" or die "Can't open cat file $cogCats\n";
	foreach my $cg (keys %catT){
		print OC join("\t",@{$catT{$cg}})."\n";
	}
	close OC;
	
	
	#die "$cogCats\n";
	
}
