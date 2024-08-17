#!/usr/bin/perl
#submits R scripts to work through phylos
#args: [GeneCat dir] [intra_phylo_dir] [abundance MGS] [mapping file[ [cores]
# perl strain_within_2.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat/intra_phylo/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat/MB2.clusters.ext.can.Rhcl.matL0.txt /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/drama4.map 1
use warnings;
use strict;


use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::geneCat qw(readGene2tax createGene2MGS);
sub sumSummaries;
my $treeSubGrpsR = getProgPaths("treeSubGrpsR");
my $RpogenS = getProgPaths("pogenStats");


my $GCd = $ARGV[0];
my $nCore = $ARGV[4];
my $refMap = $ARGV[3];#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/drama4.map";
my $FMGpD = "$GCd/MGS/phylo";
$FMGpD = $ARGV[1];# if (@ARGV > 1);
my $abMatrix = $ARGV[2];
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
my $numCores = 40;

my %dirs;my %destDs; my %baseD;

print "Reading dirs..\n";

opendir DIR, $FMGpD;
#loop over intra-phylo dir and check for file presence..
while ( my $entry = readdir DIR ) {
    next if $entry eq '.' or $entry eq '..';
    next unless -d $FMGpD . '/' . $entry;
	next unless (-d "$FMGpD/$entry/phylo/");
	#my $destD = "$FMGpD/$entry/within/";
	#system "cp $destD/$entry.nwk $FMGpD/$entry/phylo/IQtree.treefile " if (-e "$destD/$entry.nwk");
	next unless (-e "$FMGpD/$entry/phylo/$defTreeFile");
	#genuine MGS phylo dir-> store in %dirs %baseD
	$dirs{$entry} = "$FMGpD/$entry/phylo/"; 
	$baseD{$entry} = "$FMGpD/$entry";
}
closedir DIR;
print "Found ".scalar(keys %dirs)." dirs with calculated tree\n";
my $cnt=-1;
my $MGstats = "$GCd/metagStats.txt";
$MGstats = "-1" unless (-e $MGstats);
my @k2d = sort keys %dirs;
foreach my $d (@k2d){#loop over MGS intra-phylo dirs, submit R analysis
	my $destD = $dirs{$d}; $destD =~ s/(.*)\/phylo/$1\/within/; 
	my $destBaseD = $dirs{$d}; $destBaseD =~ s/(.*)\/phylo/$1\//; 
	$destDs{$d} = $destD;
	
	#next; 
	next if ( #did script already finish analysis? -> skip dir
			-e "$destBaseD/codeml/WithinStrainDiv.txt" && 
			-e "$destD/$d.Ranalysis.log" && 
			-e "$destD/$defTreeFileBase.analysis.Rdata");
	#die "$destD\n";
	
	
	$cnt++;
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
	if ($BinN<30){$nCore = 5} else {$nCore = 1;}
	my $cmd = ". /g/bork3/home/hildebra/env/mini4/etc/profile.d/conda.sh\nconda activate r_env\n";
	$cmd .= "$treeSubGrpsR $destD ../phylo/$defTreeFile $d $OG $refMap $MGstats $abMatrix $nCore > $destD/$d.Ranalysis.log\n";
	
	if (0){#rerun popgen stats??
		$cmd .= "$RpogenS $destBaseD $refMap $destBaseD/codeml/ $destBaseD/MSA/clnd/ 10,20,30,100,200,500\n";
	}
	
	#print $cmd;
	#next;
	#system $cmd."\n";
	#$QSBoptHR->{useLongQueue} = 1;
	print "$d: "; 
	my ($dep,$qcmd) = qsubSystem($destD."Ranalysis.sh",$cmd,$nCore,"10G","R$cnt","","",1,[],$QSBoptHR);
	#die " $destD\n";
	#last if ($cnt > 5);
}

if ($cnt>0){
	
	#automate
	
	die "\n\nwaiting for R analysis to finish before subclustering step\n";
} 


$FMGpD = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/MGS/";

if (1){#get within strain nuc div
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
foreach my $d (@k2d){
	my $clsts = "$destDs{$d}/${d}.Ranalysis.log";
	if ($cpDir ne ""){
		system "mkdir -p $cpDir/$d" unless (-d "$cpDir/$d");
		system "cp $destDs{$d}/${d}* $cpDir/$d";
	}
	my $SCtrig=0;
	open I,"<$clsts" or die "can't open R summary $clsts\n";
	while (my $line = <I>){
		chomp $line; 
		if ($line =~ m/QTL0.1,0.5,0.75,0.9: ([\d\.e\-]+),([\d\.e\-]+),([\d\.e\-]+),([\d\.e\-]+)/){
			$TS{$d}{qtl01} = $1; $TS{$d}{qtl05} = $2; $TS{$d}{qtl09} = $4;
		#Remaining days: 577.5 / 577.5 
		} elsif ($line =~ m/Remaining days: ([\d\.e\-NA]+) \/ ([\d\.e\-NA]+)/){
			$TS{$d}{SurvRemainDays} = $1;$TS{$d}{SurvRemainDMax} = $2;
		} elsif ($line =~ m/survival Abundance p:  ([\d\.e\-NA]+)/){
			$TS{$d}{SurvAbund} = $1;
		} elsif ($line =~ m/RegSurvFit: Abundance: ([\d\.e\-NA]+) ([\d\.e\-NA]+) Age: ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){
			$TS{$d}{SurvRegAbundZ} = $1; $TS{$d}{SurvRegAbundP} = $2; $TS{$d}{SurvRegAgeZ} = $3;$TS{$d}{SurvRegAgeP} = $4;
		#Microdiversity:1.22969837587007
		} elsif ($line =~ m/Microdiversity:([\d\.e\-NA]+)/){
			$TS{$d}{Microdiversity} = $1;
		#Within fam trans v2 (0,V,H): same: 0,0,0 diff: 2,0,1
		} elsif ($line =~ m/Within fam trans v2 \(0,V,H\): same: (\d+),(\d+),(\d+) diff: (\d+),(\d+),(\d+)/){
			$TS{$d}{WiFam2Yes} = $1;$TS{$d}{WiFam2YesV} = $2;$TS{$d}{WiFam2YesH} = $3;
			$TS{$d}{WiFam2No} = $4;$TS{$d}{WiFam2NoV} = $5;$TS{$d}{WiFam2NoH} = $6;
			
		} elsif ($line =~ m/Within fam trans v3 \(VMC,VCM,VFC,VCF\): (\d+),(\d+),(\d+),(\d+)/){
			$TS{$d}{VM2C} = $1;$TS{$d}{VC2M} = $2;$TS{$d}{VF2C} = $3;$TS{$d}{VC2F} = $4;
		
		} elsif ($line =~ m/Within fam trans v4 \(Sc,Sv,Dc,Dv\): (\d+),(\d+),(\d+),(\d+)/){
			$TS{$d}{SamStrCae} = $1;$TS{$d}{SamStrVag} = $2;$TS{$d}{DifStrCae} = $3;$TS{$d}{DifStrVag} = $4;
		} elsif ($line =~ m/Within fam trans age v5 \(mean,Csec,Vag\): ([\d\.e\-NAa]+),([\d\.e\-NAa]+),([\d\.e\-NAa]+)/){
			$TS{$d}{meanTransfAge} = $1;$TS{$d}{meanTransfAgeCaes} = $2;$TS{$d}{meanTransfAgeVagi} = $3;
		
		} elsif ($line =~ m/Within fam dis\((\d+),(\d+)\): ([\d\.e\-NA]+);([\d\.e\-NA]+);([\d\.e\-NA]+)  Between fam dis:([\d\.e\-NA]+) pval:([\d\.e\-NA]+) total:(\d+)/){
			$TS{$d}{WiFamYes} = $1;$TS{$d}{WiFamNo} = $2; $TS{$d}{WiFamYesD} = $3;$TS{$d}{WiFamNoD} = $4;
			$TS{$d}{WiFamYesNoD} = $5; $TS{$d}{BeFamD} = $6; $TS{$d}{BeFamPval} = $7;$TS{$d}{WiBeFamN}=$8;
		} elsif ($line =~ m/QTL1 0.1,0.5,0.75,0.9: ([\d\.e\-NA]+),([\d\.e\-NA]+),([\d\.e\-NA]+),([\d\.e\-NA]+)/){
			$TS{$d}{qtl101} = $1; $TS{$d}{qtl105} = $2; $TS{$d}{qtl109} = $4;
		} elsif ($line =~ m/Colless: ([\d\.e\-NA]+) sackin: ([\d\.e\-NA]+)/){
			$TS{$d}{Colless} = $1; $TS{$d}{sackin} = $2; 
		} elsif ($line =~ m/Nodes < 2e-5: (\d+) and >  (\d+) /){
			$TS{$d}{NodesLess2e5} = $1; $TS{$d}{NodesMore2e5} = $2; 
		} elsif ($line =~ m/Nodes < 0.0001: (\d+) and >  (\d+) /){
			$TS{$d}{NodesLess1e4} = $1; $TS{$d}{NodesMore1e4} = $2; 
		} elsif ($line =~ m/Nodes < 0.001: (\d+) and >  (\d+) /){
			$TS{$d}{NodesLess1e3} = $1; $TS{$d}{NodesMore1e3} = $2; 
		} elsif ($line =~ m/Nodes < 0.01: (\d+) and >  (\d+) /){
			$TS{$d}{NodesLess1e2} = $1; $TS{$d}{NodesMore1e2} = $2; 
		} elsif ($line =~ m/Nodes < mean: (\d+) and >  (\d+) /){
			$TS{$d}{NodesLessMean} = $1; $TS{$d}{NodesMoreMean} = $2; 
		} elsif ($line =~ m/Mean\/Median tip length: ([\d\.e\-NA]+) \/ ([\d\.e\-NA]+)/){
			$TS{$d}{meanTipL} = $1; $TS{$d}{medianTipL} = $2; 
		} elsif ($line =~ m/Mean\/Median internal length: ([\d\.e\-NA]+) \/ ([\d\.e\-NA]+)/){
			$TS{$d}{meanIntL} = $1; $TS{$d}{medianIntL} = $2; 
		} elsif ($line =~ m/QTL4 0.1,0.5,0.75,0.9: ([\d\.e\-NA]+),([\d\.e\-NA]+),([\d\.e\-NA]+),([\d\.e\-NA]+)/){
			$TS{$d}{qtl401} = $1; $TS{$d}{qtl405} = $2; $TS{$d}{qtl409} = $4;
		} elsif ($line =~/Stable monophyl strain dis: ([\d\.e-]+);([\d\.e-]+) ; SC=(\d+) DC=(\d+) others=(\d+)/){ #same cluster based on closest neighbor def
			$TS{$d}{monoPdis} = $1;$TS{$d}{monoPdis2} = $2; $TS{$d}{SCp} = $3; $TS{$d}{DCp} = $4; $TS{$d}{TCp} = $5;
		} elsif ($line =~/Stable monophyl 4-fol strain dis: ([\d\.e\-NA]+)/){ #same cluster based on closest neighbor def
			$TS{$d}{monoPdis4F} = $1;
		
		} elsif ($line =~/Stable conccurent samples \(weighted\):(\d+)\(([\d\.e\-NA]+),([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+),([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurS} = $1;$TS{$d}{ConcurSTime} = $2; $TS{$d}{ConcurSTimKid} = $3; $TS{$d}{ConcDiff} = $4; $TS{$d}{ConcDiffTime} = $5;$TS{$d}{ConcDiffTimeKid} = $6;
		} elsif ($line =~/Stable conccurent samples deliv Vag:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurSvag} = $1;$TS{$d}{ConcurSTimevag} = $2; $TS{$d}{ConcDiffvag} = $3; $TS{$d}{ConcDiffTimevag} = $4;
		} elsif ($line =~/Stable conccurent samples deliv Csec:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurScae} = $1;$TS{$d}{ConcurSTimecae} = $2; $TS{$d}{ConcDiffcae} = $3; $TS{$d}{ConcDiffTimecae} = $4;
		} elsif ($line =~/Stable conccurent samples ABx:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurSabx} = $1;$TS{$d}{ConcurSTimeabx} = $2; $TS{$d}{ConcDiffabx} = $3; $TS{$d}{ConcDiffTimeabx} = $4;
		} elsif ($line =~/Stable conccurent samples not ABx:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurSnab} = $1;$TS{$d}{ConcurSTimenab} = $2; $TS{$d}{ConcDiffnab} = $3; $TS{$d}{ConcDiffTimenab} = $4;

		} elsif ($line =~/Stable conccurent samples not ABx Kids:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurSnabKid} = $1;$TS{$d}{ConcurSTimenabKid} = $2; $TS{$d}{ConcDiffnabKid} = $3; $TS{$d}{ConcDiffTimenabKid} = $4;
		} elsif ($line =~/Stable conccurent samples not ABx Adults:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurSnabAdu} = $1;$TS{$d}{ConcurSTimenabAdu} = $2; $TS{$d}{ConcDiffnabAdu} = $3; $TS{$d}{ConcDiffTimenabAdu} = $4;
		} elsif ($line =~/Stable conccurent samples ABx Kids:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurSabxKid} = $1;$TS{$d}{ConcurSTimeabxKid} = $2; $TS{$d}{ConcDiffabxKid} = $3; $TS{$d}{ConcDiffTimeabxKid} = $4;
		} elsif ($line =~/Stable conccurent samples ABx Adults:(\d+)\(([\d\.e\-NA]+)\), different:(\d+)\(([\d\.e\-NA]+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{ConcurSabxAdu} = $1;$TS{$d}{ConcurSTimeabxAdu} = $2; $TS{$d}{ConcDiffabxAdu} = $3; $TS{$d}{ConcDiffTimeabxAdu} = $4;



		} elsif ($line =~/Stable samples overall \(baby\):(\d+)\((\d+)\), different:(\d+)\((\d+)\)/){ #same cluster based on closest neighbor def
			$TS{$d}{SCpN} = $1;$TS{$d}{SCpNkid} = $2; $TS{$d}{DCpN} = $3; $TS{$d}{DCpNkid} = $4; 

		} elsif ($line =~/Stable monophyl 1-fol strain dis: ([\d\.e\-NA]+)/){ #same cluster based on closest neighbor def
			$TS{$d}{monoPdis1F} = $1;
		} elsif ($line =~/SurvivalCheck SC:(\d+) DC:(\d+) TC:(\d+)/){#same cluster based on strict distance cutoffs (within tree)
		
			if ($SCtrig==0){
				$SCtrig=1;
				$TS{$d}{SC1} = $1; $TS{$d}{DC1} = $2; $TS{$d}{TC1} = $3;
			} else {
				$TS{$d}{SC2} = $1; $TS{$d}{DC2} = $2; $TS{$d}{TC2} = $3;
			}
		} elsif ($line =~ m/([\d\.e]+) tips in tree/){
			$TS{$d}{treeTips} = $1; 
		} elsif ($line =~ m/On individual level tree: med_root_dis:([\d\.e]+) IQR_root_dis:([\d\.e]+)/){
			$TS{$d}{medRtD} = $1; $TS{$d}{iqrRtD} = $2; 
#Mantel test on geo
#Testing strain diffs baby vs adult P=1e-04 R=0.304865487746621
		} elsif ($line =~ m/Testing country one individual \/ family P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{CntrFPTP} = $1; $TS{$d}{CntrFPTR} = $2; $TS{$d}{CntrFPTRadj} = $3;$TS{$d}{CntrFPTF} = $4; 
		} elsif ($line =~ m/Testing country one individual \/ family \| age P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{CntrFPTPbA} = $1; $TS{$d}{CntrFPTRbA} = $2; $TS{$d}{CntrFPTRadjbA} = $3; $TS{$d}{CntrFPTFbA} = $4; 
			
		}elsif ($line =~ m/Testing strain diffs baby vs adult P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{BabyAdult_P} = $1; $TS{$d}{BabyAdult_R} = $2; $TS{$d}{BabyAdult_Radj} = $3;$TS{$d}{BabyAdult_F} = $4; 
		}elsif ($line =~ m/Testing strain diffs baby vs adult \| country P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{BabyAdult_PbC} = $1; $TS{$d}{BabyAdult_RbC} = $2; $TS{$d}{BabyAdult_RadjbC} = $3; $TS{$d}{BabyAdult_FbC} = $4; 

		}elsif ($line =~ m/Testing strain diffs Gender \| country P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{Gender_P} = $1; $TS{$d}{Gender_R} = $2; $TS{$d}{Gender_Radj} = $3; $TS{$d}{Gender_F} = $4; 
		}elsif ($line =~ m/Testing strain diffs DeliveryMode P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{Deliv_P} = $1; $TS{$d}{Deliv_R} = $2; $TS{$d}{Deliv_Radj} = $3;  $TS{$d}{Deliv_F} = $4; 
		}elsif ($line =~ m/Testing strain diffs DeliveryMode \| country P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{Deliv_PbC} = $1; $TS{$d}{Deliv_RbC} = $2; $TS{$d}{Deliv_RadjbC} = $3;  $TS{$d}{Deliv_FbC} = $4; 
		}elsif ($line =~ m/Testing country babies \(individuals\) P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{CntrBPTP} = $1; $TS{$d}{CntrBPTR} = $2; $TS{$d}{CntrBPTRadj} = $3;$TS{$d}{CntrBPTF} = $4; 
		}elsif ($line =~ m/Testing country adults \(individuals\) P=([\d\.e-]+) R=([\d\.e-]+) Radj=([\d\.e-]+) F=([\d\.e-]+)/){
			$TS{$d}{CntrAPTP} = $1; $TS{$d}{CntrAPTR} = $2; $TS{$d}{CntrAPTRadj} = $3; $TS{$d}{CntrAPTF} = $4; 
			
			
#Survival Stats..  #  coef: -21.03133 , XYES 
		}elsif ($line =~ m/survival  Country  p:  ([\d\.e-]+)  coef: ([\d\.e-]+) , (\S+)/){
			$TS{$d}{surv_cntry_P} = $1; 
		}elsif ($line =~ m/survival  DeliveryMode  p:  (\S+)  coef: (\S+) , (\S+)/){
			$TS{$d}{surv_Deliv_P} = $1; $TS{$d}{surv_Deliv_Coef} = $2; $TS{$d}{surv_Deliv_Term} = $3; 
		}elsif ($line =~ m/survival  Gender  p:  ([\d\.e-]+)  coef: ([\d\.e-]+) , (\S+)/){
			$TS{$d}{surv_gendr_P} = $1; 
		}elsif ($line =~ m/survival  Age  p:  (\S+)  coef: ([\d\.e-]+) , (\S+)/){
			$TS{$d}{surv_Age_P} = $1; 
		}elsif ($line =~ m/survival  AntibioticsSbjAny  p:  (\S+)  coef: (\S+) , (\S+)/){
			$TS{$d}{surv_Antib_P} = $1; $TS{$d}{surv_Antib_Coef} = $2; $TS{$d}{surv_Antib_Term} = $3; 
		}elsif ($line =~ m/survival  AntibioticsSbjAft  p:  (\S+)  coef: (\S+) , (\S+)/){
			$TS{$d}{surv_AntibL_P} = $1; $TS{$d}{surv_AntibL_Coef} = $2; $TS{$d}{surv_AntibL_Term} = $3; 
		}elsif ($line =~ m/survival  AgeCat  p:  ([\d\.e-]+)  coef: ([\d\.e-]+) , (\S+)/){
			$TS{$d}{surv_AgeCat_P} = $1; 
		}elsif ($line =~ m/Survival at 50,100,200,365:  ([\d\.e-]+), ([\d\.e-]+), ([\d\.e-]+), ([\d\.e-]+)/){
			$TS{$d}{surv_T50} = $1;$TS{$d}{surv_T100} = $2;$TS{$d}{surv_T200} = $3;$TS{$d}{surv_T365} = $4;
		}elsif ($line =~ m/Median Survival: ([\d\.e-]+) Robust Mean: ([\d\.e-]+)/){
			$TS{$d}{surv_Median} = $1; $TS{$d}{surv_Mean} = $1; 


		#Mantel tests..
		}elsif ($line =~ m/Mantel Test BMI \| country: P=([\d\.e-]+) R=([\d\.e-]+)/){
			$TS{$d}{BMI_P} = $1; $TS{$d}{BMI_R} = $2; 
		}elsif ($line =~ m/Mantel Test Age \| country: P=([\d\.e-]+) R=([\d\.e-]+)/){
			$TS{$d}{Age_P} = $1; $TS{$d}{Age_R} = $2; 
		}elsif ($line =~ m/Testing strain diffs Antibtiotics P=([\d\.e-]+) R=([\d\.e-]+)/){
			$TS{$d}{AB_straindiff_P} = $1; $TS{$d}{AB_straindiff_R} = $2; 
		}elsif ($line =~ m/Testing strain diffs Antibtiotics \| country P=([\d\.e-]+) R=([\d\.e-]+)/){
			$TS{$d}{AB_straindiff_PbC} = $1; $TS{$d}{AB_straindiff_RbC} = $2; 

#Distance category tests
		}elsif ($line =~ m/AllDisCorrelogP: ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){# ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){
			$TS{$d}{AllDisCorLogP1} = $1;$TS{$d}{AllDisCorLogP2} = $2;$TS{$d}{AllDisCorLogP3} = $3;$TS{$d}{AllDisCorLogP4} = $4;
			$TS{$d}{AllDisCorLogP5} = $5;$TS{$d}{AllDisCorLogP6} = $6;$TS{$d}{AllDisCorLogP7} = $7;
		}elsif ($line =~ m/AllDisCorrelogR: ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){# ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){
			$TS{$d}{AllDisCorLogR1} = $1;$TS{$d}{AllDisCorLogR2} = $2;$TS{$d}{AllDisCorLogR3} = $3;$TS{$d}{AllDisCorLogR4} = $4;
			$TS{$d}{AllDisCorLogR5} = $5;$TS{$d}{AllDisCorLogR6} = $6;$TS{$d}{AllDisCorLogR7} = $7;
		}elsif ($line =~ m/BabyDisCorrelogP: ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){# ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){
			$TS{$d}{BabyDisCorLogP1} = $1;$TS{$d}{BabyDisCorLogP2} = $2;$TS{$d}{BabyDisCorLogP3} = $3;$TS{$d}{BabyDisCorLogP4} = $4;
			$TS{$d}{BabyDisCorLogP5} = $5;$TS{$d}{BabyDisCorLogP6} = $6;$TS{$d}{BabyDisCorLogP7} = $7;
		}elsif ($line =~ m/BabyDisCorrelogR: ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){# ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){
			$TS{$d}{BabyDisCorLogR1} = $1;$TS{$d}{BabyDisCorLogR2} = $2;$TS{$d}{BabyDisCorLogR3} = $3;$TS{$d}{BabyDisCorLogR4} = $4;
			$TS{$d}{BabyDisCorLogR5} = $5;$TS{$d}{BabyDisCorLogR6} = $6;$TS{$d}{BabyDisCorLogR7} = $7;
		}elsif ($line =~ m/AdultDisCorrelogP: ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){# ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){
			$TS{$d}{AdulDisCorLogP1} = $1;$TS{$d}{AdulDisCorLogP2} = $2;$TS{$d}{AdulDisCorLogP3} = $3;$TS{$d}{AdulDisCorLogP4} = $4;
			$TS{$d}{AdulDisCorLogP5} = $5;$TS{$d}{AdulDisCorLogP6} = $6;$TS{$d}{AdulDisCorLogP7} = $7;
		}elsif ($line =~ m/AdultDisCorrelogR: ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){# ([\d\.e\-NA]+) ([\d\.e\-NA]+)/){
			$TS{$d}{AdulDisCorLogR1} = $1;$TS{$d}{AdulDisCorLogR2} = $2;$TS{$d}{AdulDisCorLogR3} = $3;$TS{$d}{AdulDisCorLogR4} = $4;
			$TS{$d}{AdulDisCorLogR5} = $5;$TS{$d}{AdulDisCorLogR6} = $6;$TS{$d}{AdulDisCorLogR7} = $7;
		} elsif ( $line =~ m/Mantel Test geodis all: P=([\d\.e-]+) R=([\d\.e-]+)/){
			$TS{$d}{CntrFMTP} = $1; $TS{$d}{CntrFMTR} = $2; 
		} elsif ( $line =~ m/Mantel Test geodis baby: P=([\d\.e-]+) R=([\d\.e-]+)/){
			$TS{$d}{CntrBMTP} = $1; $TS{$d}{CntrBMTR} = $2; 
		} elsif ( $line =~ m/Mantel Test geodis adult: P=([\d\.e-]+) R=([\d\.e-]+)/){
			$TS{$d}{CntrAMTP} = $1; $TS{$d}{CntrAMTR} = $2; 
			
		} elsif ( $line =~ m/Pvalues strain abu:([\d\.e-]+)\(([^\)]+)\), survival time:([\d\.e-]+)\(([^\)]+)\)/){
			$TS{$d}{strain_abu_hilo_p} = $1;$TS{$d}{strain_abu_hilo_cat} = $2; $TS{$d}{strain_surv_hilo_p} = $3;$TS{$d}{strain_surv_hilo_cat} = $4; 
		} elsif ( $line =~ m/ Most abundant, abu: ([\d\.e-]+) , surv time: ([\d\.e-]+)/){
			$TS{$d}{HiStrain_abu_hilo_p} = $1; $TS{$d}{HiStrain_surv_hilo_p} = $2; 
		} elsif($line =~ m/Stable in age cats:([\d,]+) Diff:([\d,]+)/){
			$TS{$d}{StabAgeCat} = $1; $TS{$d}{DiffAgeCat} = $2; 
		}elsif($line =~ m/Weighted stable in age cats:([\d,\.]+) Diff:([\d,\.]+)/){
			$TS{$d}{StabWeiAgeCat} = $1; $TS{$d}{DiffWeiAgeCat} = $2; 


		} elsif($line =~ m/Stable in age cats subset:([\d,]+) Diff:([\d,]+)/){
			$TS{$d}{StabAgeCatSub} = $1; $TS{$d}{DiffAgeCatSub} = $2; 
		}elsif($line =~ m/Weighted stable in age cats subset:([\d,\.]+) Diff:([\d,\.]+)/){
			$TS{$d}{StabWeiAgeCatSub} = $1; $TS{$d}{DiffWeiAgeCatSub} = $2; 


		} elsif($line =~ m/Stable in age cats not AB:([\d,]+) Diff:([\d,]+)/){
			$TS{$d}{StabAgeCatNab} = $1; $TS{$d}{DiffAgeCatNab} = $2; 
		}elsif($line =~ m/Weighted stable in age cats not AB:([\d,\.]+) Diff:([\d,\.]+)/){
			$TS{$d}{StabWeiAgeCatNab} = $1; $TS{$d}{DiffWeiAgeCatNab} = $2; 
		} elsif($line =~ m/Stable in age cats ABx:([\d,]+) Diff:([\d,]+)/){
			$TS{$d}{StabAgeCatAbx} = $1; $TS{$d}{DiffAgeCatAbx} = $2; 
		}elsif($line =~ m/Weighted stable in age cats ABx:([\d,\.]+) Diff:([\d,\.]+)/){
			$TS{$d}{StabWeiAgeCatAbx} = $1; $TS{$d}{DiffWeiAgeCatAbx} = $2; 



		}
	}
	close I;
}
my @hds = (
		"treeTips",
		"monoPdis","monoPdis2","monoPdis1F","monoPdis4F","SCp","DCp","TCp","SC1","DC1","TC1","SC2","DC2","TC2","medRtD","iqrRtD",
		"SCpN","SCpNkid","DCpN","DCpNkid","Microdiversity","ConcDiffTimeKid","ConcDiffTime","ConcDiff","ConcurSTimKid","ConcurSTime","ConcurS",
		"SamStrCae","SamStrVag","DifStrCae","DifStrVag","meanTransfAge","meanTransfAgeCaes","meanTransfAgeVagi",
		"ConcurSvag","ConcurSTimevag","ConcDiffvag","ConcDiffTimevag","ConcurScae","ConcurSTimecae","ConcDiffcae","ConcDiffTimecae",
		"ConcurSabx","ConcurSTimeabx","ConcDiffabx","ConcDiffTimeabx","ConcurSnab","ConcurSTimenab","ConcDiffnab","ConcDiffTimenab",
		"ConcurSabxKid","ConcurSTimeabxKid","ConcDiffabxKid","ConcDiffTimeabxKid","ConcurSnabKid","ConcurSTimenabKid","ConcDiffnabKid","ConcDiffTimenabKid",
		"ConcurSabxAdu","ConcurSTimeabxAdu","ConcDiffabxAdu","ConcDiffTimeabxAdu","ConcurSnabAdu","ConcurSTimenabAdu","ConcDiffnabAdu","ConcDiffTimenabAdu",
		"qtl01","qtl05","qtl09","qtl101","qtl105","qtl109","qtl401","qtl405","qtl409","meanTipL","medianTipL",
		"Colless","sackin","meanIntL","medianIntL","WiFamYes","WiFamNo","WiFamYesD","WiFamNoD","WiFamYesNoD","BeFamD","BeFamPval","WiBeFamN",   "WiFam2Yes","WiFam2YesV","WiFam2YesH","WiFam2No","WiFam2NoV","WiFam2NoH",
		"VM2C","VC2M","VF2C","VC2F",
		"SurvAbund","SurvRegAbundZ","SurvRegAbundP","SurvRegAgeZ","SurvRegAgeP",
		"NodesLess2e5","NodesMore2e5","NodesLess1e4","NodesMore1e4","NodesLess1e3","NodesMore1e3",
		"NodesLess1e2","NodesMore1e2","NodesLessMean","NodesMoreMean",
		"CntrFPTP","CntrFPTR","CntrFPTRadj","CntrFPTF","CntrFPTPbA","CntrFPTRbA","CntrFPTRadjbA","CntrFPTFbA",
		"CntrBPTP","CntrBPTR","CntrBPTRadj","CntrBPTF","CntrAPTP","CntrAPTR","CntrAPTRadj","CntrAPTF","CntrFMTP","CntrFMTR","CntrBMTP","CntrBMTR","CntrAMTP","CntrAMTR",
		"DiffAgeCat","StabAgeCat","StabWeiAgeCat","DiffWeiAgeCat",
		"StabAgeCatSub","DiffAgeCatSub","StabWeiAgeCatSub","DiffWeiAgeCatSub",
		"DiffAgeCatAbx","StabAgeCatAbx","StabWeiAgeCatAbx","DiffWeiAgeCatAbx",
		"DiffAgeCatNab","StabAgeCatNab","StabWeiAgeCatNab","DiffWeiAgeCatNab",
		#"strain_abu_hilo_p","strain_abu_hilo_cat","strain_surv_hilo_p","strain_surv_hilo_cat","HiStrain_abu_hilo_p","HiStrain_surv_hilo_p",
		"surv_cntry_P","surv_Deliv_P","surv_Deliv_Coef","surv_Deliv_Term","surv_Antib_P","surv_Antib_Coef","surv_Antib_Term","surv_AntibL_P","surv_AntibL_Coef","surv_AntibL_Term",
		"surv_gendr_P","surv_Age_P","surv_AgeCat_P","surv_T50","surv_T100","surv_T200","surv_T365",
		"surv_Mean","surv_Median",

		"AllDisCorLogP1","AllDisCorLogP2","AllDisCorLogP3","AllDisCorLogP4","AllDisCorLogP5","AllDisCorLogP6","AllDisCorLogP7",
		"AllDisCorLogR1","AllDisCorLogR2","AllDisCorLogR3","AllDisCorLogR4","AllDisCorLogR5","AllDisCorLogR6","AllDisCorLogR7",
		"BabyDisCorLogP1","BabyDisCorLogP2","BabyDisCorLogP3","BabyDisCorLogP4","BabyDisCorLogP5","BabyDisCorLogP6","BabyDisCorLogP7",
		"BabyDisCorLogR1","BabyDisCorLogR2","BabyDisCorLogR3","BabyDisCorLogR4","BabyDisCorLogR5","BabyDisCorLogR6","BabyDisCorLogR7",
		"AdulDisCorLogP1","AdulDisCorLogP2","AdulDisCorLogP3","AdulDisCorLogP4","AdulDisCorLogP5","AdulDisCorLogP6","AdulDisCorLogP7",
		"AdulDisCorLogR1","AdulDisCorLogR2","AdulDisCorLogR3","AdulDisCorLogR4","AdulDisCorLogR5","AdulDisCorLogR6","AdulDisCorLogR7",

		"Gender_P","Gender_R","Gender_Radj","Gender_F","Deliv_P","Deliv_R","Deliv_Radj","Deliv_F",
		"Deliv_PbC","Deliv_RbC","Deliv_RadjbC","Deliv_FbC",
		"BabyAdult_P","BabyAdult_R","BabyAdult_Radj","BabyAdult_F","BabyAdult_PbC","BabyAdult_RbC","BabyAdult_RadjbC","BabyAdult_FbC",
		"BMI_P","BMI_R","Age_P","Age_R","AB_straindiff_P","AB_straindiff_R","AB_straindiff_PbC","AB_straindiff_RbC");
open O,">$FMGpD/Rsummary.tab";
print O "\t".join("\t", @hds)."\n";
foreach my $d (@k2d){
	print O "$d";
	foreach my $kk (@hds){
		if (!exists($TS{$d}{$kk}) || !defined($TS{$d}{$kk})){
			#print "can't find $d $kk\n"; 
			print O "\tNA";
		} else {
			my $val = $TS{$d}{$kk}; #$val =~s/-+$//;
			print O "\t$val";
		}
	}
	print O "\n";
}
close O;


die "$FMGpD/Rsummary.tab";


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




