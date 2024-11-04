#!/usr/bin/perl
#submits R scripts to work through phylos
#args: [GeneCat dir] [intra_phylo_dir] [abundance MGS] [mapping file[ [cores]
# perl strain_within_2.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat/intra_phylo/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat/MB2.clusters.ext.can.Rhcl.matL0.txt /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/drama4.map 1
use warnings;
use strict;

use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::TamocFunc qw ( getFileStr );

sub sumSummaries;
my $strainStatsR = getProgPaths("treeSubGrpsR");
my $RpogenS = getProgPaths("pogenStats");
my $MFdir = getProgPaths("MGSTKDir");
my $vizPhylos = getProgPaths("vizPhylosSign_R");
my $treewasRun_R = getProgPaths("treewasRun_R");
my $processTreewas_R = getProgPaths("processTreewas_R");

#.21: added $DiscTests $ContTests
#.22: R interface updated
#.23: Perl interface updated
#.24: added familyVar and groupStabilityVars arguments for stability calculations 
#.25: 3.4.24 added Kasia's scoary scripts
#.26: 11.7.24 replaced scoary scripts with treewas
my $version = 0.26;

my $rewriteRanalysis = 0; my $doSubmit = 1;

my $DiscTests =""; my $ContTests = "";
my $familyVar = ""; my $groupStabilityVars = "";
my $GCd = "";#$ARGV[0];
my $nCore = 4;#$ARGV[4];
my $nCoreHeavy = 32;
my $refMap = "";#$ARGV[3];#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/drama4.map";
#my $FMGpD = "$GCd/MGS/phylo";
my $FMGpD = "";#$ARGV[1];# if (@ARGV > 1);
my $abMatrix = "";#$ARGV[2];
#discrete and continous tests done on each MGS' strains.. based on column header in map file, comma separated
#$DiscTests = $ARGV[5]; $ContTests = $ARGV[6]; 
my $individualVar = "AssmblGrps";

GetOptions(
	"GCd=s"          => \$GCd,
	"map=s"          => \$refMap,
	"submit=i"       => \$doSubmit,
	"reSubmit=i"     => \$rewriteRanalysis,
	"cores=i"        => \$nCore,
	"Hcores=i"        => \$nCoreHeavy, #for heavy jobs (like scoary)
	"FMGdir=s"     => \$FMGpD,
	"MGSmatrix=s"      => \$abMatrix,
	"DiscTests=s" => \$DiscTests, #discrete categoried (from map) to be tested via perMANOVA 
	"ContTests=s"     => \$ContTests, #continous variables (from map) to be tested for strain phylo signal
	"familyVar=s"      => \$familyVar, #column name in metadata containing family id
	"groupStabilityVars=s"      => \$groupStabilityVars, #column names of categories used for calculation of resilience and persistence
	"individualVar=s"      => \$individualVar, #column name specifying individual IDs, AssmblGrps by default
);


my $cpDir = "";#$GCd/MGS/R_analysis/";
#my $outD
my $defTreeFile = "IQtree_allsites.treefile";
my $defTreeFileBase = "IQtree_allsites";
#die;
die "MGS phylo dir doesn't exist!\n$FMGpD\n" unless (-d $FMGpD);

my $QSBoptHR = emptyQsubOpt(1,"");
$QSBoptHR->{tmpSpace} = 0;
my $bts = getProgPaths("buildTree_scr");
my $SaSe = "|";
my $SaSe2 = "\|"; #for regex
my $numCores = 40; #used for phylos..
my $RsummaryTab = "$FMGpD/Rsummary.tab";


my %dirs;my %destDs; my %baseD;

print "Strain analysis v $version\n";



if ($rewriteRanalysis){ #faster to do once for all..
	print "\nWARNING:: Rewriting strain2 results!\n" ;
	#print "Removing old strain2 analysis..\n";
	system "rm -rf $FMGpD/*/within ;rm -f $RsummaryTab"
}


print "Reading dirs..\n";
opendir DIR, $FMGpD;
#loop over intra-phylo dir and check for file presence..
my %sizTrees;
while ( my $entry = readdir DIR ) {
    next if $entry eq '.' or $entry eq '..';
    next unless -d $FMGpD . '/' . $entry;
	next unless (-d "$FMGpD/$entry/phylo/");
	#my $destD = "$FMGpD/$entry/within/";
	#system "cp $destD/$entry.nwk $FMGpD/$entry/phylo/IQtree.treefile " if (-e "$destD/$entry.nwk");
	my $sizTree = -s "$FMGpD/$entry/phylo/$defTreeFile";
	next unless ($sizTree);
	#genuine MGS phylo dir-> store in %dirs %baseD
	$dirs{$entry} = "$FMGpD/$entry/phylo/"; 
	$baseD{$entry} = "$FMGpD/$entry";
	$sizTrees{$entry} = $sizTree;
}

closedir DIR;
print "Found ".scalar(keys %dirs)." dirs with calculated tree\n";
my $cnt=-1; my $curBatch=0; my $batchSize = 0; my $submitted=0;my @jobs;
my $MGstats = "$GCd/metagStats.txt";
$MGstats = "-1" unless (-e $MGstats);
my $treeAbsent = 0;
#my @k2d = sort keys %dirs;
my @k2d = sort { $sizTrees{$b} <=> $sizTrees{$a} } keys(%sizTrees);

my $cmd = "ulimit -s 20000\n";my $destD =""; my $wrHead=0; 
foreach my $d (@k2d){#loop over MGS intra-phylo dirs, submit R analysis
	$cnt++;
	$destD = $dirs{$d}; $destD =~ s/(.*)\/phylo/$1\/within/; 
	my $destBaseD = $dirs{$d}; $destBaseD =~ s/(.*)\/phylo/$1\//; 
	$destDs{$d} = $destD;
	#my $locTree = "$destD ../phylo/$defTreeFile"; #two args in one..
	if (!-e "$dirs{$d}/$defTreeFile"){
		$treeAbsent++;
		next;
	}
	system("rm -rf $destD/* ;") if ($rewriteRanalysis && -d $destD);
	#next; 
	next if ( #did script already finish analysis? -> skip dir
			#-e "$destBaseD/codeml/WithinStrainDiv.txt" && 
			
			 -e "$destD/$d.Ranalysis.log" 
			&& -e "$destD/$d.analysis.txt" 
			#&& -e "$destD/$defTreeFileBase.analysis.Rdata"
			#&& -e "$destD/$defTreeFileBase.analysis.txt"
			);
	#die "$destD\n";
	#next if ($cnt == 0);
	#die "$destD/$d.Ranalysis.log\n";
	systemW "rm -rf $destD/*";
	system "mkdir -p $destD" unless (-d $destD);
	my $tmp = `cat $1/data.log`;
	$tmp =~ m/OG:(.*)/; my $OG = $1;
	#system "cp $dirs{$d}/$defTreeFile $destD/$d.nwk";
	my $treeDef = "$dirs{$d}/$defTreeFile";
	my $BinN = 1000;
	if ($d =~ m/MB2bin(\d+)/){$BinN = $1;}
	if ($BinN<30){$nCore = 5} else {$nCore = 5;}
	
	$cmd .= "echo \"At tree $d\"\n";
	$wrHead=1 if ( $cnt == 0);
	$cmd .= "$strainStatsR --path $destD --tree ../phylo/$defTreeFile --taxN $d --outgroup $OG --map $refMap --metagStats $MGstats --abMat $abMatrix --ncore $nCore --siteMode 1 --MFDir $MFdir --wrColNms $wrHead --discPermTests \"$DiscTests\" --contPermTests \"$ContTests\" --familyCol \"$familyVar\" --groupStabilityVars \"$groupStabilityVars\" > $destD/$d.Ranalysis.log\n";
	$wrHead=0;
	if (0){#rerun popgen stats??
		$cmd .= "$RpogenS $destBaseD $refMap $destBaseD/codeml/ $destBaseD/MSA/clnd/ 10,20,30,100,200,500\n";
	}
	
	#print $cmd;
	#next;
	#system $cmd."\n";
	#$QSBoptHR->{useLongQueue} = 1;
	print "$d: "; 
	$curBatch++;
	if ($curBatch > $batchSize){
		my ($dep,$qcmd) = qsubSystem($destD."Ranalysis.sh",$cmd,$nCore,"20G","R$cnt","","",1,[],$QSBoptHR);
		#die " $destD\n";
		push(@jobs,$dep);
		$curBatch = 0; $cmd="";
		$submitted++;
	}
	#die;
	#last if ($cnt > 5);
}
if ($curBatch > 0){
	my ($dep,$qcmd) = qsubSystem($destD."Ranalysis.sh",$cmd,$nCore,"10G","R$cnt","","",1,[],$QSBoptHR);
	$curBatch = 0; $cmd="";
	push(@jobs,$dep);
}

if ($submitted>0){ #wait for all submitted R scripts, then continue in script
	
	#automate
	
	print "\n\nwaiting for R analysis to finish before subclustering step\n";
	qsubSystemJobAlive( \@jobs,$QSBoptHR );

} 

print "$treeAbsent phylos absent\n";

if (0){#get within strain nuc div
	my $countStrains=0;
	open O,">$FMGpD/withinStrain.tab";
	foreach my $d (@k2d){
		my $destBaseD = $dirs{$d}; $destBaseD =~ s/(.*)\/phylo/$1\//; 
		my $strainFile = "$destBaseD/within/IQtree_allsites.strains.txt";
		next unless (-e $strainFile);
		print O "$d\tNA\n";
		my $wiStF="$destBaseD/codeml/WithinStrainDiv.txt";
		open I,"<$wiStF" or die "can't find file $wiStF\n";
		while (<I>){ print O $_;}
		close I;
		my %seensStrains;
		open I,"<$strainFile" or die "can;t open strain file $strainFile\n";
		while (my $lin = <I>){
			chomp $lin; my @spl = split /\t/,$lin;
			#print $lin."\n"; # if (@spl < 2 || !defined($spl[1]));
			$seensStrains{$spl[1]} = 1;# unless ($seensStrains{$spl[1]} eq "NA");
		}
		close I;
		$countStrains += (scalar(%seensStrains)-1);
	}
	close O;
	print "Total Strains seen: $countStrains\n";
}

my %TS; 
#--------------------------------------------------------------
#summary of dnds
if (0){
	#fubar summaries
	sumSummaries("hyphy.fubar.txt","$FMGpD/fubar.tab");
	sumSummaries("hyphy.fubar.s10.txt","$FMGpD/fubar.s10.tab");
	sumSummaries("hyphy.fubar.s20.txt","$FMGpD/fubar.s20.tab");
	sumSummaries("hyphy.fubar.s30.txt","$FMGpD/fubar.s30.tab");
	sumSummaries("hyphy.fubar.s100.txt","$FMGpD/fubar.s100.tab");
	sumSummaries("hyphy.fubar.s200.txt","$FMGpD/fubar.s200.tab");
	sumSummaries("hyphy.fubar.s500.txt","$FMGpD/fubar.s500.tab");
	sumSummaries("hyphy.fubar.unID.txt","$FMGpD/fubar.unID.tab");
	#popstats
	sumSummaries("PopStats.txt","$FMGpD/PopStats.tab");
	sumSummaries("PopStats.10.txt","$FMGpD/PopStats.10.tab");
	sumSummaries("PopStats.20.txt","$FMGpD/PopStats.20.tab");
	sumSummaries("PopStats.30.txt","$FMGpD/PopStats.30.tab");
	sumSummaries("PopStats.100.txt","$FMGpD/PopStats.100.tab");
	sumSummaries("PopStats.200.txt","$FMGpD/PopStats.200.tab");
	sumSummaries("PopStats.500.txt","$FMGpD/PopStats.500.tab");
	#sumSummaries("PopStats.500.txt","$FMGpD/PopStats.unID.tab");
	sumSummaries("hyphy.Theta.log","$FMGpD/Theta.tab");
}

#die;
#summary  of R stats
#create summary tables


if (1 || !-e $RsummaryTab){
	system "rm -f $RsummaryTab";
	#reset output report file
	open O,">$RsummaryTab" or die $!;  close O;

	foreach my $d (@k2d){
		my $clsts = "$destDs{$d}/${d}.Ranalysis.log";
		if ($cpDir ne ""){
			system "mkdir -p $cpDir/$d" unless (-d "$cpDir/$d");
			system "cp $destDs{$d}/${d}* $cpDir/$d";
		}
		my $SCtrig=0;
		
		my $TXTreport = "$destDs{$d}/${d}.analysis.txt";

		if (-e $TXTreport){
			my $cmd = "cat $TXTreport >> $RsummaryTab;";
			system $cmd;
		}
		next;
	}
}





# functional enrichments of strains in conditions defined by user
my $funCmd = "";
my $treewasOut = "$FMGpD/GeneEnrich/";
my $treewasOutfile = "$treewasOut/treeWAS_results.csv";
my $summaryOutfile = "$treewasOut/treeWAS_results_functions.csv";
my $MGSd = $FMGpD; $MGSd =~ s/\/[^\/]+[\/]+$/\//;
$funCmd .= "mkdir -p $treewasOut\n";
$funCmd .= "#1st command: run treewas job\n";
$funCmd .= "$treewasRun_R --gene_cat_dir \"$GCd\" --n_threads $nCoreHeavy --metadata_vars \"$groupStabilityVars\" -o \"$treewasOut\" --mgs_dir \"$MGSd\" --metadata_file \"$refMap\" -r \"$MFdir\" -i \"$individualVar\" \n";
$funCmd .= "#2nd command: process results\n";
$funCmd .= "$processTreewas_R -i \"$treewasOutfile\" --gene_cat_dir \"$GCd\" --annot_files \"NOG,CZy,KGM\", --out_file \"$summaryOutfile\" --n_threads $nCoreHeavy -r \"$MFdir\" \n";


my ($dep,$qcmd) = qsubSystem($treewasOut."treeWAS.sh",$funCmd,$nCoreHeavy,"6G","treewas","","",1,[],$QSBoptHR);


#"$strainStatsR $destD ../phylo/$defTreeFile $d $OG $refMap $MGstats $abMatrix $nCore 1 $MFdir $wrHead $DiscTests $ContTests > $destD/$d.Ranalysis.log\n";
#my $taxFile = "$GCd/Anno/Tax/GTDBmg_MGS/specI.tax";
my $taxFile = "$FMGpD/../Annotation/MGS.GTDB.LCA.tax";
my $cmdPic = "$vizPhylos $RsummaryTab $taxFile $FMGpD phylo $MFdir $refMap -1\n";

print "Printing figures of most significant phylogenies\nThis might take several hours..\n";
print $cmdPic;
system $cmdPic;

#die "$FMGpD/Rsummary.tab";

print "\nFinished with strain postprocessing\n";
exit (0);

#--------------------------------------------------------------
#now starts the real work, take cluster files and make trees based on this subset.
#this part depends on pre created sets of fastas containing all genes for a given species
#then it will resort the cat file to only include samples deduced from R tree
foreach my $d (@k2d){
	my $clsts = "$destDs{$d}/${d}.cl_IndFam.txt";
	my %clusters;
	open I,"<$clsts" or die "can't open $clsts\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		push(@{$clusters{$spl[1]}}, $spl[0]);
	}
	close I;
	
	my $fnFile = "$baseD{$d}/allFNAs.fna";
	my $faFile = "$baseD{$d}/allFAAs.faa";
	my $catFile = "$baseD{$d}/all.cat";
	if (!-e $fnFile || !-e $catFile){
		die "Can't find input files:\n$catFile\n$fnFile\n";
	}
	my %genePres;
	open I,"<$catFile" or die "Can't open the category file $catFile\n";
	while (<I>){
		chomp;
		my @spl = split /\t/;
		foreach my $gn (@spl){
			$gn =~ m/^(.*)$SaSe2(.*)$/;
			$genePres{$2}{$1} = 1; #{COG}{SMPL}
		}
	}
	close I;
	my $cnt=0;
	#sort out which genes go into same tree..
	foreach my $clN (keys %clusters){
		my @smpls = @{$clusters{$clN}};
		my $nOD = "$baseD{$d}/Cl_$clN/";
		system "mkdir -p $nOD" unless (-d $nOD);
		my $clOF = "$nOD/cl$clN.cat";
		#TODO: add 3 smpls from other cluster as outliers..
		my $outg = "";
		open O,">$clOF";
		#each line one cog
		foreach my $cog (keys %genePres){
			my @cp;
			foreach my $s (@smpls){
				if (exists($genePres{$cog}{$s})){ #gene exists for sample, put in cat file
					push (@cp,"$s$SaSe$cog");
				}
			}
			if (@cp > 0 ){
				print O join("\t",@cp)."\n";
			}
		}
		close O;
		#cat file created, now submit job
		my $Tcmd= "$bts -fna $fnFile -aa $faFile -smplSep '\\$SaSe' -cats $clOF -outD $nOD -runIQtree 1 -iqFast 1 -runFastTree 0 -cores $numCores  ";
		$Tcmd .= "-AAtree 0 -bootstrap 000 -NTfiltCount 3000 -NTfilt 0.1 -NTfiltPerGene 0.6 -runRaxMLng 0 -minOverlapMSA 2 -MSAprogram 2 -AutoModel 0 \n";
		#die "$cmd\n" if ($cnt ==10);
		$QSBoptHR->{useLongQueue} = 1;
		my ($dep,$qcmd) = qsubSystem($nOD."subtree_$clN.sh",$Tcmd,$numCores,"1G","FT$cnt","","",1,[],$QSBoptHR);
		$cnt ++;
		if ($outg ne ""){
			open O,">$nOD/TODO";
			close O;
		}
	}
	
}






exit(0); #don't do FST for now..
#--------------------------------------------------------------
#calculate FST between countries / cities..
#depends on Alex's script for calculating FST
my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;
#$mapF = $GCd."LOGandSUB/inmap.txt" if ($mapF eq "");
my ($hr1,$hr2) = readMapS($mapF,-1,"Country");
my %map = %{$hr1}; my %AsGrps = %{$hr2};
my $FSTbin = getProgPaths("FSTpy");
foreach my $d (@k2d){
	my $clsts = "$destDs{$d}/${d}.cl_IndFam.txt";
	my $cntrFcnt=0;
	my %cntr;
	open I,"<$clsts" or die "can't open $clsts\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		if (exists($map{$spl[0]}{"Country"})){
			my $tcnt=$map{$spl[0]}{"Country"};
			push(@{$cntr{$tcnt}},$spl[0]) ;
			$cntrFcnt++;
		} else {
			print "not found: $spl[0]\n";
		}
	}
	close I;
	print "$d : $cntrFcnt country\n";
	#start FST script from alex
	my $FSTd = "$baseD{$d}/FST/";
	system "mkdir -p $FSTd" unless (-d $FSTd);
	my @presCntrs = sort keys %cntr;
	open O, ">$FSTd/FST.cntry.guide" or die $!;
	foreach my $cn (@presCntrs){
		print O "$cn\t".join(",",@{$cntr{$cn}}) . "\n";
	}
	close O;
	my $refpop = "USA";
	if (!exists($cntr{USA})){$refpop = $presCntrs[0];}
	my $Fcmd = "$FSTbin -i $baseD{$d}/MSA/MSAli.fna -b 100 -p $FSTd/FST.cntry.guide --refpop $refpop -o  $FSTd/FST.cntry.test\n";
	die "$Fcmd\n";
}





sub sumSummaries($ $){
	my ($inF,$outF) = @_;
	open O1,">$outF" or die "Can't open $outF\n";
	my $first=1;
	foreach my $d (@k2d){
		#next;
		my $destD = $dirs{$d}; $destD =~ s/(.*)\/phylo/$1\/codeml/; 
		if (-e "$destD/$inF"){
			open I,"<$destD/$inF";my $tmp = <I>; 
			print O1 "\t$tmp" if ($first);while (<I>){print O1 "$d\t$_";	}close I; 
		} else {
			print "missing: $destD/$inF\n";
		}
		$first=0;
	}
	close O1; 
}




