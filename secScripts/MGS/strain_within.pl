#!/usr/bin/perl
#This will build a phylo tree for each MGS (between samples)
#./strain_within.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//MB2.clusters.ext.can.Rhcl.mgs 10 /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//between_phylo/prunned.nwk 1 1
#./strain_within.pl /g/bork3/home/hildebra/data/SNP/GCs/alienGC2 /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Binning/MetaBat/MB2.clusters.obs 12 
use warnings;
use strict;

use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw(fileGZe fileGZs readClstrRev systemW median mean readMapS readFasta getAssemblPath);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystem2 qsubSystemJobAlive qsubSystemWaitMaxJobs);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::geneCat qw(readGene2tax createGene2MGS);
my $bts = getProgPaths("buildTree_scr");
my $neiTree = getProgPaths("neighborTree");
sub extractFNAFAA2genes;
sub histoMGS;
sub readGenesSample_Singl;

die "Not enough args!\n" unless (@ARGV > 1);

#v.14: reworked massively how many genes get included
#v.15: included lessons learned from MGS.pl v0.21 
#v.16 added familyVar and groupStabilityVars arguments for stability calculations
#v.17: considerations to improve speed of intial fna/faa extractions..
my $version = 0.16;

#input args..
my $GCd = "";#$ARGV[0];
my $MGSfile = "";#$ARGV[1];
my $numCores = 1;#$ARGV[2];
my $maxCores = -1;
my $onlySubmit =0;#extract genes anew?
my $reSubmit=0;
my $treeFile = "";
my $doSubmit=0;
my $subMode="";
my $maxNGenes = 300;
my $presortGenes = 1600;
my $checkMaxNumJobs = 600;
my $useGTDBmg = "GTDB";
my $redoSubmissionData = 0;
my $deepRepair = 0;
my $rmMSA = 1; #argument passed to buildTree5.pl 
my $contTests = ""; my $discTests = ""; #stat tests to be given to strain_within_2.2.pl
my $familyVar = ""; my $groupStabilityVars = "";

#define local files..
my $lSNPdir="SNP"; 
my $lConsFNA = "genes.shrtHD.SNPc.MPI.fna.gz";
my $lConsFAA = "proteins.shrtHD.SNPc.MPI.faa.gz";


#$treeFile = $ARGV[3] if (@ARGV > 3);$onlySubmit = $ARGV[4] if (@ARGV > 4);
#$doSubmit = $ARGV[5] if (@ARGV > 5);$subMode = $ARGV[6] if (@ARGV > 6);


GetOptions(
	"GCd=s"          => \$GCd,
	"MGS=s"          => \$MGSfile,
	"submit=i"       => \$doSubmit,
	"onlySubmit=i"   => \$onlySubmit, #submit only jobs, or also recreate input fna/faa files? (can take days)
	"reSubmit=i"     => \$reSubmit, #for all MGS: resubmit tree phylo building
	"cores=i"        => \$numCores, #not used any longer..
	"maxCores=i"     => \$maxCores, #superseedes -cores, will dynamically allocate num cores based on input file size, if defined
	"MGSphylo=s"     => \$treeFile,
	"subMode=s"      => \$subMode,
	"presortGenes=i" => \$presortGenes, #how many potential genes to include, of the original MGS (receovered will vary strongly  between samples)
	"maxGenes=i"     => \$maxNGenes, #how many genes to try to include? -> will be decided on each samples
	"MGset=s"        => \$useGTDBmg,
	"redoSubmissionData=i" => \$redoSubmissionData,  #for all MGS: will resubmit phylo and rebuild the fna/faa files..
	"deepRepair=i"   => \$deepRepair, #for missing MGS phylos: will resubmit phylo and rebuild fna/faa 
	"rmMSA=i"        => \$rmMSA, #remove MSA, to save diskspace
	"ContTests=s"      => \$contTests, #continous stat tests to be handed to next step (just a passthrough)
	"DiscTests=s"      => \$discTests, #discrete stat tests to be handed to next step (just a passthrough)
	"familyVar=s"      => \$familyVar, #column name in metadata containing family id
	"groupStabilityVars=s"      => \$groupStabilityVars, #column names of categories used for calculation of resilience and persistence
);


my $subsSmpl = -1;
my $tmpD = getProgPaths("globalTmpDir",0);
$GCd =~ m/.*\/([^\/]+\/?)/;
$tmpD .= "/$1/";
#die $tmpD."\n";

my $takeAll = 0;
my $conspecificSpThr = 0.1;
my $multiCpyThr = 0.2;
my $mode = "MGS";
my %genesWrite; #keep stats/track

$mode = "FMG" if ($MGSfile eq "");
my $useSuperTree = 0;
if ($mode eq "FMG"){$takeAll = 0;}

if ($takeAll){$maxNGenes = -1;$mode="MGSall"; $doSubmit = 0;}


my $bindir = $MGSfile;
$bindir =~ s/[^\/]+$//; 
my $outD =  $bindir;#"$GCd/$mode/intra_phylo/";
$outD .= "/intra_phylo/";

print "\n!! WARNING !!: RESUBMISSION mode selected (will resubmit MSA + phylos even for already completed MGS) !!\n\n" if ($reSubmit);
print "\n!! WARNING !!: REDOSUBMISSIONDATA mode selected (will redo and resubmit MSA + phylos even for already completed MGS) !!\n\n" if ($redoSubmissionData);



print "============= Strains from MGS v$version =============\n";
print "Creating within species strains for ${mode}s in $GCd\n";
print "GC dir: $GCd\nIn Cluster: $MGSfile\n# cores: $numCores\n";
#print "Ref tree: $treeFile\n";
print "Using tree $treeFile to create automatically outgroups\n" if ($treeFile ne "");
print "Outdir: $outD\n";
print "MGs: $useGTDBmg\n";
print "Using $presortGenes genes from each MGS for location\n";
print "Deep repariing remaining submission files\n" if ($deepRepair);
if ($takeAll){print "**************** Take all genes MGS mode\n";}
else {print "Using first $maxNGenes genes found per sample\n";}
print "==============================================\n";
if ($onlySubmit){print "Only submission mode\n";
} else {
	print "Creation of strain genes, old data might be deleted!\nDo you want to continue? (10s wait, use Ctrl-c to abort)\n"; #sleep 10;
}
#die;


#STONES
system "mkdir -p $outD/stones/" unless (-d "$outD/stones/");
my $inputChk = "$outD/stones/0.fileChk.sto";

#set up some base paths specific to pipeline..
my $FNAstdof = "allFNAs.fna"; my $FAAstdof = "allFAAs.faa";
my $LINKstdof = "link2GC.txt"; my $CATstdof = "all.cat";
my $fnaSNPf = "/SNP/genes.shrtHD.SNPc.MPI.fna.gz";
my $aaSNPf = "/SNP/proteins.shrtHD.SNPc.MPI.faa.gz";

my $xtraGuids = 

my $QSBoptHR = emptyQsubOpt($doSubmit,"",$subMode);
my %QSBopt = %{$QSBoptHR};

system "mkdir -p $outD" unless (-d $outD);

if (!-e $MGSfile.".srt"){
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
#die; 

my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;chomp $mapF;
#$mapF = $GCd."LOGandSUB/inmap.txt" if ($mapF eq "");
my ($hr1,$hr2) = readMapS($mapF,-1);
my $hr3; my $hr4;
my %map = %{$hr1}; my %AsGrps = %{$hr2};
my %AGlist;
#get all samples in assembly group, but only last in mapgroup
my @samples = @{$map{opt}{smpl_order}};
my $fileAbsent = 0; my $inputChkd = 0;
$inputChkd =1 if (-e "$inputChk");
foreach my $smpl (@samples){ # just check that files are there..
#and fill %AGlist .. so always let run..
	#last if ($inputChkd);
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
	#DEBUG
	#next;
	
	#check if SNP file is present
	next if ($onlySubmit && -e $inputChk);
	my $cD = $map{$smpl}{wrdir}."/";
	#my $tarF = $cD."/SNP/genes.shrtHD.SNPc.MPI.fna.gz";
	my $tarF = $cD."/$lSNPdir/$lConsFNA";
	my $tarF2 = $cD."/$lSNPdir/$lConsFAA";
	if (!-e $tarF || !-e $tarF2){
		print "Can't find SNP  file: $cD\n" ;
		$fileAbsent = 1;
	}
}
unless (-e "$inputChk"){
	if ($fileAbsent){
		print "Not all required input present\n" ;
	} else {
		print "All samples have SNP calls\n";
	}
	system "touch $inputChk" ;
}
#foreach (sort keys %AGlist) {   print "$_ : @{$AGlist{$_}}\n";}die;

my %replN; 
#my %allFNA; my %allFAA; #big hash with all genes in @allGenes
my %gene2genes;
my %cl2gene2;
my %FNAref; my %FAAref;
#my %SIcat;


print "\n\n----------------------------------------------------\nPart I:: extracting relevant core MGS genes (SNP consensus called) from original assemblies\n----------------------------------------------------\n\n";



#this file actually often only contains a fraction of associated genes..
my $gene2taxF = "$GCd/FMG/gene2specI.txt";
if ($useGTDBmg eq "GTDB"){
	$gene2taxF = "$GCd/GTDBmg/gene2specI.txt";
}

#this file has much more info.. better..
if ($mode eq "MGS" || $mode eq "MGSall"){
	#print "$MGSfile\n";
	$gene2taxF = createGene2MGS($MGSfile,$GCd);
	print "Using MGS from $MGSfile, adding eggNOG in: $gene2taxF\n";
} 
#die;
#print "$gene2taxF\n";


#these might be very limited number of genes here..
($hr1,$hr2,$hr3,$hr4) = readGene2tax($gene2taxF,$presortGenes);#$maxNGenes);
my %SIgenes=%{$hr1}; my %Gene2COG=%{$hr2}; my %Gene2MGS = %{$hr3}; my %COGprios = %{$hr4};
my %SIdirs;
my %SIgenes_OG; #later reads in SIgenes again, but no restriction to length 
my @specis = sort(keys(%SIgenes));
#sort specis by numbers
my %sis; foreach (@specis){m/(\d+)$/; $sis{$_}=$1;}
@specis = sort {$sis{$a} <=> $sis{$b} } keys %sis;
#die "@specis\n";
my $cnt=0; my $SaSe = "|"; my $dirsArePrepped = 1; my $allCatFileE = 1;
foreach my $SI (@specis){ #loop creates per specI file structure to run buildTreeScript on..
	#PART I: create fasta files required by tree
	my $outD2 = "$outD/$SI/";
	$SIdirs{$SI} = $outD2;
	#print "$outD2\n";
	$dirsArePrepped =0 unless (-d $outD2);
	$allCatFileE = 0 if (!-e  "$SIdirs{$SI}/$CATstdof");
	if (-d $outD2 && $onlySubmit == 0){#don't delete folders if we want to submit a job later..
		system "rm -rf $outD2/*";
	}
	system "mkdir -p $outD2" unless (-d $outD2);
	
}
 #my $refFAA = "$GCd/compl.incompl.95.prot.faa";	my $hr = readFasta($refFAA,1,"\\s",\%Gene2COG);  %FAAref = %{$hr};	die "read ". scalar(keys %FAAref)." genes\nTher are ".scalar(keys %Gene2COG)." norm genes\n";



if ($dirsArePrepped == 0 || $onlySubmit == 0){
	print "Preparing base strain alignments, per MGS\nThis might take a good while..\n";
	($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",0,\%Gene2COG);my %cl2gene = %{$hr2}; $hr1 = {};
	#stores alt names (wihtout M4__ at end)
	#read binning based on SpecI's 
	my @allGenes;
	my $goodGene =0; my $badGene =0;
	#in this process we can also check for multi genes
	foreach my $gene (keys %Gene2COG){
		my $geneStr = $cl2gene{$gene};
		$geneStr =~ s/>//g;
		my @genegenes = split /,/,$geneStr;
		#$cl2gene2{$gene} = \@genegenes;
		my %tmpGen; #my %gnCnt;
		#finding the sample might be slow in the start.. but will pay off long term
		foreach my $sg (@genegenes){
			my @spls = split /__/,$sg;
			#die "{$gene}{$spls[0]}\n";
			push(@{$tmpGen{$spls[0]}}, $sg);
			#$gnCnt{$spls[0]} ++;
		}
		#my $totSmpl = scalar(keys(%gnCnt));
		my $badGn =0;
		foreach my $sm(keys %tmpGen){
			$badGn++ if (scalar(@{$tmpGen{$sm}}) > 1);
			#if ($gnCnt{$sm} > 1){$badGn++;}
		}
		#fraction of samples where gene was double. Genes were this is obviously too high need to be excluded (multi copy)
		#my $frac = $badGn/$totSmpl;
		my $frac = $badGn/@genegenes;
		#print "$frac = $badGn/$totSmpl \n";
		if ($frac > $multiCpyThr){
			$badGene++;
			next;
		}
		$goodGene ++;
		foreach my $sm(keys %tmpGen){
			next if (scalar(@{$tmpGen{$sm}}) > 1); #don't copy this gene in this sample .. no decission which is the right one..
			#here gene is added in
			@{$cl2gene2{$sm}{$gene}} = @{$tmpGen{$sm}};
		}
	}
	print "Found $goodGene / ".($badGene + $goodGene)." genes useable (single copy)\n";
	
	undef %cl2gene; #lessen mem


	#and extract the corresponding fna/ faa from every other dir.. main single core work
	#this will also determine how many genes per MGS are now extracted..
	extractFNAFAA2genes(\%cl2gene2);#@allGenes);
	
	%cl2gene2 = (); #no longer needed, delete
	
	
	print "\nGene extraction & redistribution finished, ready to proceed to phylogeny jobs\n";
	@allGenes = ();
	#	exit(0);

	
	#write logs to found genes etc.
	foreach my $SI (@specis){
		my $outD2 = $SIdirs{$SI}; my $llogF="$outD2/geneFnd.log";
		open LL ,">$llogF" or die "can't open local gene log $llogF\n";
		if (exists( $genesWrite{$SI} )) {
			print LL "Total genes write $SI: $genesWrite{$SI}\n";
		} else {
			print LL "Total genes write $SI: 0\n";
		}
		close LL;
	}
}
my $geneCatLoaded=0;
if (1 && (!$allCatFileE || $deepRepair || !$dirsArePrepped || $onlySubmit == 0 || $redoSubmissionData == 1)){
	#also read reference gene seqs (for outgroup)
	my $refFNA = ""; my $refFAA = "";
	if ($mode eq "MGS" || $mode eq "MGSall"){
		print "Reading reference genecat genes, to create outgroup sequences\n";
		$refFNA = "$GCd/compl.incompl.95.fna"; $refFAA = "$GCd/compl.incompl.95.prot.faa";
		$geneCatLoaded=1;
	} elsif ($mode eq "FMG"){
		print "reading FMG ref genes..";
		$refFNA = "$GCd/FMG/COG*.fna"; $refFAA = "$GCd/FMG/COG*.faa";
	}
	($hr1,$hr2,$hr3,$hr4) = readGene2tax($gene2taxF);
	%SIgenes_OG=%{$hr1}; my %Gene2COG_OG=%{$hr2}; 
	my $hr = readFasta($refFAA,1,"\\s",\%Gene2COG_OG); %FAAref = %{$hr};
	$hr = readFasta($refFNA,1,"\\s",\%Gene2COG_OG); %FNAref = %{$hr};
	print "read ". scalar(keys %FNAref)." genes\n";
	print "done\n";
	

}



print "\n\n----------------------------------------------------\nPart II:: resort .cat files, submit intraStrain phylogeny\n----------------------------------------------------\n\n";


die "Tree for outgroup specified, but file not found:$treeFile\nAborting..\n" if  ($treeFile ne "" && !-e $treeFile);

#die;
#go through every SpecI;
$cnt=0; my $lcnt=0; my @jobs;
foreach my $SI (@specis){ #loop creates per specI file structure to run buildTreeScript on..
	qsubSystemWaitMaxJobs($checkMaxNumJobs);
	#print "$SI  XX ";
	$lcnt++;
	#print "$SI\n";
	#next unless ($lcnt>10);
	#PART I: create fasta files required by tree
	my $outD2 = $SIdirs{$SI};
	my $treeStone = "$outD2/treeDone.sto";
	#next if (-e "$outD2/phylo/IQtree_allsites.treefile");
	#print "$outD2\n";
	my $multiSmpl=0;
	my $FNAtf = "$outD2/$FNAstdof"; my $FAAtf = "$outD2/$FAAstdof";
	my $CATtf = "$outD2/$CATstdof"; #my $Linkf = "$outD2/$LINKstdof";
	my $IQtreef= "$outD2/phylo/IQtree_allsites.treefile";
	
	#job done..
	next if (-e $treeStone && -e $IQtreef && !$reSubmit && !$redoSubmissionData);
	print "${SI}::"; 
	
	my $outgS = "";my $OG = "";
	if (-e "$outD2/data.log"){$OG = `cat $outD2/data.log`; chomp $OG; $OG=~s/^OG://;}
	my $contPhylo = 1; $contPhylo = 0 if ($reSubmit || $redoSubmissionData);
	
	#main command to build within species strain tree.. missing outgroup so far ($outgS)
	
	my $inputFNAsize = fileGZs($FNAtf) / (1024 * 1024); #size in MB
	if ( $inputFNAsize  > 50 ){ #only if FNA is > X mb
		$QSBoptHR->{useLongQueue} = 1 ;
	}
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
	my $totMem = int($inputFNAsize /250*40)+10;$totMem = 20 if ($totMem < 20);
	my $numCoreL = $numCores;	
	if ($maxCores >0){ #scale cores according to used memory size
		$numCoreL = int($maxCores * sqrt($totMem) / sqrt(250));
		$numCoreL = 4 if ($numCoreL < 4);		$numCoreL = $maxCores if ($numCoreL > $maxCores);
	}

	my $Tcmd= "$bts -fna $FNAtf -aa $FAAtf -smplSep '\\$SaSe' -cats $CATtf -outD $outD2  -runIQtree 1 -runFastTree 0 -cores $numCoreL  "; 
	$Tcmd .= "-AAtree 0 -bootstrap 0 -NTfiltCount 400 -NTfilt 0.07 -NTfiltPerGene 0.5 -GenesPerSpecies 0.1 -runRaxMLng 0 -minOverlapMSA 2 ";
	$Tcmd .= "-subsetSmpls $subsSmpl -fracMaxGenes90pct 0.7 "; #concentrate on almost complete gene groups.. can yield more samples overall and speeds up calc..
	$Tcmd .= "-rmMSA $rmMSA -gzInput 1 "; #save more diskspace..
	$Tcmd .= "-SynTree 0 -NonSynTree 0 -MSAprogram 4 -continue $contPhylo -AutoModel 0 -iqFast 1 -superTree $useSuperTree ";
	$Tcmd .= "-runDNDS 0 -runTheta 0 -tmpD $tmpD/$SI/ -map $mapF ";
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
	} else {
		#print OC $SIcatLoc{$cog}{$smpl};				print OC "\t".$SIcatLoc{$cog}{$smpl};		my $ng = "$OG$SaSe$cog";
		#print "Reconstructing tmp cat file.";
		open ICT,"<$CATtf" or die "Can't open cat file $CATtf\n";
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
		print "$SI:: $catLines cat lines, $cntItems items: $CATtf\n";
	}
	
	
	#my @curCogs = sort keys %{$SIcat{$SI}};
	my @curCogs = sort keys %SIcatLoc;
	if (scalar(@curCogs) < 10){
		die "$SI error: \@curCogs is empty or < 10 members\n";
	}
	#print "COGs: $curCogs[0] $curCogs[1]\n";
	
	#include outgroup?
	if ($treeFile ne ""){
		my $call = "$neiTree $treeFile $SI";
		print "$call\n";
		#print "$call";
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
		unless (exists($SIgenes{$OG})){
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
		#my $tar = $SIgenes{$SI}{$cog};
		#my $new=1; my %smplsSeen; my $multiSmpl1 = 0;
		#print OL "$cog\t$tar\t".scalar @genes . "\t".join(",",@genes)."\n";
		if ($OG ne "" && exists($SIgenes_OG{$OG}{$cog})){#deal with outgroup..
			#die "can't find gene $SIgenes{$OG}{$cog}" unless (exists($FNAref{$SIgenes{$OG}{$cog}}));
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
	#system "echo \"$OG\" > $outD2/numG_". scalar(keys %{$SIgenes{$SI}})."_smpls_$multiSmpl";
	system "echo \"OG:$OG\" > $outD2/data.log";
	$multiSmpl = scalar(keys %uniqSmpls);
	if ($multiSmpl>2){
		print "$SI: multiSmpls:\t$multiSmpl\tpotential genes: ". scalar(@curCogs) ."\tcores:$numCoreL mem:$totMem\n";
	} else {print "\n$SI: too few samples\n";next;}
	
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
qsubSystemJobAlive( \@jobs,\%QSBopt );
print "\nAll done for $cnt Bins\nRun strain_within_2.pl for summary stats:\n";

my $outDX =  $MGSfile;#"$GCd/$mode/intra_phylo/";
$outDX =~ s/[^\/]+$//;
my $MGSabundance = "$GCd/Anno/Tax/GTDBmg_MGS/specI.mat";
$MGSabundance = "$bindir/Annotation/Abundance/MGS.matL7.txt";

my $strain2Scr = getProgPaths("MGS_strain2_scr");

my $nxtCmd = "$strain2Scr -GCd $GCd -FMGdir $outD -MGSmatrix $MGSabundance -map $mapF -cores 4 -Hcores $maxCores -reSubmit 0 -DiscTests \"$discTests\" -ContTests \"$contTests\" -familyVar \"$familyVar\" -groupStabilityVars \"$groupStabilityVars\" \n"; #$GCd/MB2.clusters.ext.can.Rhcl.matL0.txt
	my ($dep,$qcmd) = qsubSystem($outD."strainAnalysis2.sh",$nxtCmd,1,"60G","2StrainSub","","",1,[],$QSBoptHR);
print "\n". $nxtCmd."\n";

exit(0);

 

#########################################################################################
#########################################################################################





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
	my ($OFstrHr, $OAstrHr, $OCstrHr, $OLstrHr, $writeLink) = @_;
	my $wrMGS=0;
	foreach my $SI (keys(%{$OFstrHr})){
		next if ($OFstrHr->{$SI} eq "");
		#handle file paths..
		my $outD2 = $SIdirs{$SI};
		my $FNAtf = "$outD2/$FNAstdof"; my $FAAtf = "$outD2/$FAAstdof";my $Linkf = "$outD2/$LINKstdof";
		my $CATtf = "$outD2/$CATstdof.tmp";
		open OF,">>$FNAtf" or die "Can't append NT file $FNAtf\n";
		open OA,">>$FAAtf" or die "Can't append AA file $FAAtf\n";
		open OL,">>$Linkf" or die "Can't append link file $Linkf\n" if ($writeLink);
		open OC,">>$CATtf" or die "Can't append to CAT file $CATtf\n";#this is only a temp file, that needs to be rewritten later..
		
		#writing strings out..
		print OF $OFstrHr->{$SI}; print OA $OAstrHr->{$SI};print OC $OCstrHr->{$SI} ;print OL $OLstrHr->{$SI};
		close OF; close OA;close OC;#tmp OC file
		close OL if ($writeLink);
		$OCstrHr->{$SI} = ""; $OFstrHr->{$SI} = ""; $OAstrHr->{$SI} = ""; $OLstrHr->{$SI} = "";
		$wrMGS++;
	}
	print "\nwrote for $wrMGS MGS data..\n";
}





#this routine hast to get genes out of each sample, that are needed
#and save them to be later written per specI
sub extractFNAFAA2genes{
	my ($cl2gene2) = @_;
	#my %cl2gene2 = %{$hr};
	my %perMGScnts;
	my %totGnes;
	#create gene to genes list
	foreach my $sm (keys %{$cl2gene2}){
		#my @locGenes;
		my $gnCnt=0;my $MGSgeneCnt=0;
		foreach my $gn (keys %{${$cl2gene2}{$sm}}){
			$totGnes{$gn} = 1;
			$gnCnt++;
			if (exists($Gene2MGS{$gn})){
				$perMGScnts{$Gene2MGS{$gn}}{$gn}=1;
				#print "1";
				$MGSgeneCnt++;
			}
		}
		#print "$sm  $gnCnt $MGSgeneCnt \n";
	}
	my @histoMGScnts ;#= values %perMGScnts;
	foreach my $MGS (keys %perMGScnts){push(@histoMGScnts,  scalar( keys( %{$perMGScnts{$MGS}} ) ));}
	#print "\n@histoMGScnts\n";
	print "Genes per MGS (prefiltering, N= ".  scalar(keys(%totGnes)) ." genes, " .scalar(keys(%perMGScnts)) . " MGS, avg " . int(0.5+scalar(keys(%totGnes))/scalar(keys(%perMGScnts))) . " genes/MGS):\n";
	#show_histogram(\@histoMGScnts,10);
	
	histoMGS(\@histoMGScnts,"Theorectical best Bin sizes: ");
	#some stats on genes/MGS
	my @srtdSmpls = sort (keys %{$cl2gene2});
	print "Extracting GC genes from " . scalar(@srtdSmpls). " dirs\n";

	
	#different way to go over genes..
	 my $smCnt=1;
	 #storage hash for raw fasta/faa/link files, needs to be written separately
	my $OCstrHR = {}; my $OFstrHR = {} ; my $OAstrHR = {} ; my $OLstrHR = {};
	#goes over every assembly group to extract SNP corrected genes that fall into each MGS
	my $writeLink = 1; my $appCnt=0;
	foreach my $sm (@srtdSmpls){
		#my @locGenes;
		my %subG; my %locMGScnt;
		my %locCl2G2 = %{$cl2gene2{$sm}};
		
		#%hash = map { $_ => 1 } @array;
		
		foreach my $gn (keys %locCl2G2){
			#put genes into hash to avoid duplicates..
			$subG{$_} = 1 foreach(@{$locCl2G2{$gn}}) ;
			#stats collection on MGS usage
			if (exists($Gene2MGS{$gn})){
				$locMGScnt{$Gene2MGS{$gn}}++;
			}
		}
		print "$smCnt/" . scalar(@srtdSmpls) ." $sm\t".scalar(keys(%subG))." genes, " . scalar(keys(%locMGScnt)). " MGS\n";
		my @histoMGScnts = values %locMGScnt;
		#foreach my $MGS (keys %locMGScnt){push(@histoMGScnts,  scalar( keys( %{$locMGScnt{$MGS}} ) ));}
		histoMGS(\@histoMGScnts, "Possible Bins in sample");
		readGenesSample_Singl(\%subG,$sm, $OFstrHR, $OAstrHR, $OCstrHR, $OLstrHR, $writeLink);
		$smCnt++; $appCnt++;
		if ($appCnt >= 50){
			appendWriteMGSgenes($OFstrHR, $OAstrHR, $OCstrHR, $OLstrHR, $writeLink);
			$appCnt=0;
		}
	}
	
	appendWriteMGSgenes($OFstrHR, $OAstrHR, $OCstrHR, $OLstrHR, $writeLink);
	$appCnt=0;
	#done at the point with gene extractions
	return;
}



sub readGenesSample_Singl{
	#go into curSpl dir and extract all marked gene reps.. 
	#write to correct format so they can be used in phylo later
	my ($subGHR,$sm, $OFstrHR, $OAstrHR, $OCstrHR, $OLstrHR, $writeLink) = @_;
	my %subG = %{$subGHR};#$_[0]};
	my $sd = $_[1]; #this is current sample
	my $sd2 = $sd;
#	my $writeLink = 1;
	#print "$sd ";
	if (exists(  $map{altNms}{$sd}  )){
		$sd2 = $map{altNms}{$sd}; $replN{$sd} = $sd2;
	}
	#check if sample in map
	unless (exists ($map{$sd2}) ) {
		print "Can't find map entry for $sd\n"; die;
	}
	my @subGKs = keys %subG;
	die "regex failed: $subGKs[0]\n" unless ($subGKs[0] =~ m/^(.*)__/);
	#find out if other samples are in the same assmblGrp..
	my @subSds = ($sd2);
	my $cAssGrp = $map{$sd2}{AssGroup};
	if (exists($AGlist{$cAssGrp})){
		@subSds = @{$AGlist{$cAssGrp}};
	}
	#print "YY @subSds : $sd2\n";
	#go into each sample, that an assembly might be associated to
	foreach my $sd3 (@subSds){
		my %locFAA; my %locFNA;my%locCSP;
		my %locMGSgenes; #keep track of genes written for each MGS..
		my $cD = $map{$sd3}{wrdir}."/";
		#print "$cD\n";
		my $rename = 0;
		$rename = 1 if ($sd2 ne $sd3);
		#print "r:$rename $sd3  ";
		my $metaGD = getAssemblPath($cD);
		if ($metaGD eq ""){die "assembly not available: $cD $sd3";}
		#get NT's
		#my $tar = $metaGD."genePred/genes.shrtHD.fna";
		my $fastaf = $cD.$fnaSNPf;
		unless (-e $fastaf){
			print "\n=====================================\nCan't find nt file $fastaf\n=====================================\n";
			next;
		}
		#print "$fastaf\n";
		#read the assemble nt and AA genes from the sample
		my $FNA = readFasta($fastaf,1,"\\s");#,\%subG);
		#my %FNA = %{$hr};
		$fastaf = $cD.$aaSNPf;
		my $FAA2 = readFasta($fastaf,0);#,"\\s",\%subG);#read full head string
		my $FAA = {};
		#my %FAA = %{$hr};
		#convert FAA hd
		my %conspSc;#read conspecific strain score from SNP consensus call..
		foreach my $k(keys %{$FAA2}){
			my $tmp = $k;
			$k =~ s/\s.*//;
			$FAA->{$k} = $FAA2->{$tmp};
			$tmp =~ m/ CSP=([0-9\.]+)/;
			if (!defined $1){print "no CSP $tmp\n";
			} else {$conspSc{$k} = $1;
			}
		}
		$FAA2 = {};
		my $geneLost=0;
		my @kks = keys %{$FNA};
		foreach my $ge (@subGKs){
			#die"YES $ge $fastaf" if ($ge =~ m/C1404_L=8071=_3/);
			if ( !exists($FNA->{$ge})){
				$geneLost++;
				#print "$ge ";
				next;
			}
			#die ;
			my $ge2= $ge;
			if ($rename){
				#die "${sd2}__/${sd3}__\n";
				unless ($ge2 =~ s/${sd}__/${sd3}__/){
					die "could not replace $sd with $sd3 in string $ge2\n";
				}
				push(@{$gene2genes{$ge}},$ge2); #only save subset..
			}
			
			$locFNA{$ge2} = $FNA->{$ge};
			$locFAA{$ge2} = $FAA->{$ge};
			$locCSP{$ge2} = $conspSc{$ge};
		}
		#fill /reset stat vector..
		foreach my $SI (@specis){$locMGSgenes{$SI} = 0;}
		
		
		#print scalar %locFAA . " genes found\n";
		$FAA = {}; $FNA = {};#%FAA = (); %FNA = (); 
		%conspSc = ();
		#some stats on gene extractions..
		my $missGene=0; my $foundGene=0; my $SInum=0; my $conspGen=0;
		
		#3rd part: genes were read and renamed.. now write them out already here to save mem overall
		foreach my $SI (@specis){
			my $OCstr=""; my $OFstr = ""; my $OAstr = ""; my $OLstr = "";
			
			my $locCnt=0; my $locConSpecGen=0;
			
			die "Can't find $SI in COGprios!\n" unless (exists($COGprios{$SI}));
			#get actual gene & gene2assmblname
			foreach my $cog (@{$COGprios{$SI}}){#keys %{$SIgenes{$SI}}){
				next unless (exists($SIgenes{$SI}{$cog}));
				my $tar = $SIgenes{$SI}{$cog};
				next unless (exists($cl2gene2{$sd}{$tar}));
				#print "yes ";
				my @genes = @{$cl2gene2{$sd}{$tar}};
				
				#write link file, but only needs to be done once.. this avoids doing this later when the cat file is written
				if ($writeLink){
					#$OLstr .= "$cog\t$tar\t".scalar @genes . "\t".join(",",@genes)."\n" ;
					$OLstr .= "$cog\t$tar\t".scalar @genes . "\t".join(",",@genes)."\n" ;
				}

				foreach my $gX ( @genes ){
					my @genes2 = ($gX);
					push (@genes2, @{$gene2genes{$gX}}) if (exists($gene2genes{$gX}));
					foreach my $g (@genes2){
						#not present for some reason.. (can happen in SNP call)
						if (!exists($locFAA{$g})){
							#if only present in nt, but not AA: skip..
							die "$g : $cog : $SI have NT but not AA\n" if (exists($locFNA{$g}));
							#print "miss: $g"; 
							#$skipGene{$SI}{$g}=1;
							$missGene++; next;
						}
						if (!exists($locCSP{$g})){print "Can't find $g CSP\n";}
						if ($locCSP{$g} > $conspecificSpThr){#too many indicators that gene is from conspecific strain
							$locConSpecGen ++ ; $conspGen++; next;
						}
						my $strCpy = ""; $strCpy = $locFAA{$g} if (exists($locFAA{$g}));
						my $AAlen = length($strCpy);
						if ($AAlen == 0){$missGene++; next;}
						my $num1 = $strCpy =~ tr/[\-Xx]//;
						if ($num1 >= ($AAlen-1)){ $missGene++; next;}
						if (exists($locMGSgenes{$SI}) && $locMGSgenes{$SI} >= $maxNGenes){next;}
						
						#write gene out
						my $ng = "$sd3$SaSe$cog"; #must contain 2 informations: 1)sampleID 2)COG 
						$OFstr .= ">$ng\n$locFNA{$g}\n";
						$OAstr .= ">$ng\n$strCpy\n";
						$locCnt++;
						#add to category for later..
						$OCstr .= "$SI\t$cog\t$sd3\t$ng\n";
						#$SIcat{$SI}{$cog}{$sd3} = $ng;
						$genesWrite{$SI}++;
						$locMGSgenes{$SI}++;
					}
				}
			}
			
			if ($locCnt == 0 || $OFstr eq ""){ #nothing to do here..
				next;
			}
			
			if ($locConSpecGen*0.05 >= ($locConSpecGen+$locCnt) ){
				$OCstr = ""; $OFstr = ""; $OAstr = ""; $OLstr = "";
				#print "skipping conspecific MGS ";
				next;
			} 
			
			if (!exists($OAstrHR->{$SI})){#set up base strings
				$OAstrHR->{$SI} = "";$OFstrHR->{$SI} = "";$OLstrHR->{$SI} = "";$OCstrHR->{$SI} = "";
			}
			#save in tmp hash (faster than opening bunch of files..
			$OAstrHR->{$SI} .= $OAstr;$OFstrHR->{$SI} .= $OFstr;$OLstrHR->{$SI} .= $OLstr;$OCstrHR->{$SI} .= $OCstr;
			$SInum ++ if ($locCnt>0);
			$foundGene+=$locCnt;
		}
		my @genesPmgs = values %locMGSgenes; 	@genesPmgs = sort { $a <=> $b}  @genesPmgs;
		histoMGS(\@genesPmgs,"Detected Bin Genes:");
		
		print "$sd3 - Missed(lost)Gs: $missGene($geneLost)\tConspec: $conspGen\tFoundGs: $foundGene/". scalar %locFAA . "\tMGS: $SInum\t";
		print "GperMGS (median,mean): " . median(@genesPmgs) . "/". int(mean(@genesPmgs)+0.5);#int($foundGene/$SInum) if ($SInum);
		print "\n";
	}
}


