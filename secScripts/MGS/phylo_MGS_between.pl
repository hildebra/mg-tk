#!/usr/bin/perl
#script to get a set of marker genes from each FMG (40 MG), extract them, and build phylo tree
#relatively simple, since can use genes directly from gene cat, no need to get SNP called genes
#can also include reference genomes to include in tree
#perl /hpc-home/hildebra/dev/Perl/MATAF3//secScripts/MGS/phylo_MGS_between.pl -GCd /ei/projects/3/3c24aae4-5ce2-4156-a31a-82d4602c2176/data/GC_PDD1/ -MGS /ei/projects/3/3c24aae4-5ce2-4156-a31a-82d4602c2176/data/GC_PDD1//Binning//MB2.clusters.ext.can.Rhcl.filt -c 10 -outD /ei/projects/3/3c24aae4-5ce2-4156-a31a-82d4602c2176/data/GC_PDD1//Binning//customRefs/ -refGenos '/hpc-home/hildebra/geneCats/Chicken2/Cultured_genomes/99_ani_dRep/*.fasta'

use warnings;
use strict;
use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::phyloTools qw(calcDisPos2 getGenoGenes getFMG readFMGdir);
use Mods::geneCat qw(calculate_spearman_correlation read_matrix correlation checkAntiOcc);

my $bts = getProgPaths("buildTree_scr");
my $vizTree = getProgPaths("vizBtwPhylo_R");


if (@ARGV < 2){
	die "Not enough input args: use \n./phylo_MGS.pl -GCd [path to GC] -MGS [MGS file]\n";
}

my $wait4job = 0;#wait till tree is done? only neccessary in pipeline...
my $numCores = 20;
my $GCd = "";#$ARGV[0];
my $btout = "";#"$GCd/MGS/phylo/";#main output dir
my $MGSfile = "";#$ARGV[1];
$numCores = 1;#$ARGV[2] if (@ARGV > 2);
my $addRefGenos="";
my $fastphylo = 0;
my $MSAprog = 4; #4:MUSCLE5, 2:mafft
my $xtraMessageInSH = "";
my $mem = 120; #memory request in GB
#$btout = $ARGV[3] if (@ARGV > 3);
#$wait4job = $ARGV[4] if (@ARGV > 4);

#options to pipeline..
GetOptions(
	"GCd=s"      => \$GCd,
	"MGS=s" => \$MGSfile,
	#"submit=i" => \$doSubmit,
	"outD=s" => \$btout,
	"c|cores=i" => \$numCores,
	"wait2finish=i" => \$wait4job,
	"refGenos=s" => \$addRefGenos,
	"fast=i" => \$fastphylo,
	"MSAprogram=i" => \$MSAprog,
	"xtraMsg=s" => \$xtraMessageInSH,
	"mem=i" => \$mem,
);

die "No input args\n" if ($GCd eq "" || $MGSfile eq "");
$btout = "$GCd/MGS/phylo/" if ($btout eq "");#main output dir


#main objects to store dna/cats
my %FAAfmg; my %FNAfmg;my %catT; 
my %MGS; my %MGSFMG; my %dblList; my %totDbls;


#is there any ref genomes to add?
if ($addRefGenos ne ""){
	my $ncore = 1;
	my @sfiles = glob($addRefGenos);
	print "External Genomes: $addRefGenos\nProcessing ".scalar(@sfiles)." genomes\n";
	my $cnt=0;
	foreach my $tarG(@sfiles){
	#my $tarG = "$tarDir/$genoN";
		my $tag = $tarG;
		$tag =~ s/.*\///;
		$tag =~ s/\.[^\.]*$//;
		#print "$tarG\n";
		my ($genes,$prots) = getGenoGenes($tarG);
		my $FMGdir = getFMG("",$prots,$genes,$ncore);
		my ($hrN,$hrA,$hrC) = readFMGdir( "$FMGdir",$tag ,".");
		$cnt++;
		my %FAA=%{$hrA};
		my %COGcat=%{$hrC};
		#add to categories..
		foreach my $k (keys %FAA){
			$FAAfmg{$k}=$FAA{$k};
			die "unkown key $k \n" unless (exists($COGcat{$k}));
			$MGSFMG{$tag}{$COGcat{$k}} = $k;
		}
		#print " $tag ";
		#die;
		#last if ($cnt>5);#DEBUG
	}
}

#read FMG designation
my %FMG2COG;
open I,"<$GCd/FMG.subset.cats" or die "Can't open $GCd/FMG.subset.cats\n";
while (<I>){
	chomp;
	my @spl = split /\t/;
	my @s2 = split /,/,$spl[2];
	foreach my $x (@s2){
		$FMG2COG{$x} = $spl[0];
	}
}
close I;

print "Found ". scalar(keys(%FMG2COG)) ." FMG genes in total gene cat\n";

#read MGS genes
my $mfcnt=0; my $mfdbl=0;
open I,"<$MGSfile" or die "Can't open MGS input\n";
while (<I>){
	chomp; my @spl = split /\t/;
	next if (@spl <2);
	my @genes = split /,/,$spl[1];
	#not needed here..
	#$MGS{$spl[0]} = \@genes;
	foreach my $x (@genes){
		if (exists($FMG2COG{$x})){
			if (exists($MGSFMG{$spl[0]}{$FMG2COG{$x}})){
				$mfdbl++;
				$dblList{$spl[0]}{$FMG2COG{$x}}{$x} = 1;
				$dblList{$spl[0]}{$FMG2COG{$x}}{$MGSFMG{$spl[0]}{$FMG2COG{$x}}} = 1;
				$totDbls{$x} = 1; $totDbls{$MGSFMG{$spl[0]}{$FMG2COG{$x}}} = 1;
			} else {
				$MGSFMG{$spl[0]}{$FMG2COG{$x}} = $x;
			}
			$mfcnt ++;
		}
	}
}
close I;
print "Assigned $mfcnt genes to MGS in ".scalar(keys(%MGSFMG))." MGS (". int(10*$mfcnt/scalar(keys(%MGSFMG)))/10 ." on average, $mfdbl double)\n";


#routine to do double checking etc and also find targets for merging gene clusters..
if (0){
	my $hr = read_matrix("$GCd/Matrix.FMG.mat","\t",\%totDbls); my %FMGm = %{$hr};
	$hr = readFasta("$GCd/FMG/COG*.fna"); %FNAfmg = %{$hr};
	my $sdir = "$GCd/MGS/dbls/"; system "mkdir -p $sdir" unless (-d $sdir);
	my $mergeList = "";
	#take care of double list
	print scalar(keys %dblList)." MGS in double list\n";
	my @hist; my @hist2;
	foreach my $mg (sort keys %dblList){
		my $dcnt=0;
		my @avgIDs;
		open O,">$sdir/$mg.fna" or die $!;
		#$hist[scalar(keys (%{$dblList{$mg}}))]++;
		my @dblKeys = sort keys %{$dblList{$mg}};
		for (my $j=0;$j<@dblKeys;$j++){
			my $c  = $dblKeys[$j];
			open O2 ,">>$sdir/$mg.$c.fna" or die $!;
			foreach my $x (keys %{$dblList{$mg}{$c}}){
				print O ">$c.$x\n$FNAfmg{$x}\n";
				print O2 ">$c.$x\n$FNAfmg{$x}\n";
			}
			close O2;
			my $avgID = calcDisPos2("$sdir/$mg.$c.fna","",1,10,"/tmp/hildebra/test/");
			$avgIDs[$j] = $avgID;
			unlink("$sdir/$mg.$c.fna");
		}
		close O;
		my $avgID = calcDisPos2("$sdir/$mg.fna","$sdir/$mg.dist",1);
		$hist[$dcnt] ++;
		#get summary stat on MGS
		my $merge=0; my $couldMerge=0; my $nomerge=0;
		for (my $i=0;$i<@avgIDs;$i++){
			if ($avgIDs[$i] > 90){ #candidate for merging..
				my @potDG = keys %{$dblList{$mg}{$dblKeys[$i]}};
				$hr = checkAntiOcc(\@potDG,\%FMGm); my %mrgMat = %{$hr};
				foreach my $mk (keys %mrgMat){
					foreach my $mk2 (keys %mrgMat){
						if ( (exists($mrgMat{$mk}{$mk2}) && $mrgMat{$mk}{$mk2} > 0.9) ||
							(exists($mrgMat{$mk2}{$mk}) && $mrgMat{$mk2}{$mk} > 0.9)
						){
							$merge++ ;
							$mergeList .= "$mk\t$mk2\n";
						} else {
							$couldMerge++;
						}
					}
				}
				
			} else {
				$nomerge++ ;
			}
		}
		print "$mg ($merge/$couldMerge/$nomerge)\t@avgIDs\n";
		foreach my $ids (@avgIDs){
			$hist2[int($ids)]++;
		}
	}
	#print histogram
	for (my $h=0;$h<@hist2;$h++){
		last;
		print "$h:$hist2[$h] ";
	}
	open O,">$GCd/FMG/merges.FMG.txt" or die $!; print O $mergeList; close O;
	print "\n";
	die;
}


print "reading FMG ref genes..";
my $hr = readFasta("$GCd/FMG/COG*.faa"); %FAAfmg = (%FAAfmg,%{$hr});
print "done\n";

system "mkdir -p $btout" unless (-d $btout);

#open ON,">$btout/all.fna"; 
open OA,">$btout/all.faa"  or die "Can't open faa out file $btout/all.faa\n"; 
my $SaSe = "|";
foreach my $mg (keys %MGSFMG){
	foreach my $cog (keys %{$MGSFMG{$mg}}){
		my $ng = "$mg$SaSe$cog";
#		print ON ">$ng\n$FNAfmg{$MGSFMG{$mg}{$cog}}\n";
		die "$MGSFMG{$mg}{$cog}\n" unless (exists( $FAAfmg{$MGSFMG{$mg}{$cog}} ));
		print OA ">$ng\n$FAAfmg{$MGSFMG{$mg}{$cog}}\n";
		push(@{$catT{$cog}},"$ng");
	}
}
close OA;

open OC,">$btout/all.cats" or die "Can't open cat file $btout/all.cats\n";
foreach my $cg (keys %catT){
	print OC join("\t",@{$catT{$cg}})."\n";
}
close OC;
my $QSBoptHR = emptyQsubOpt(1,"");
$QSBoptHR->{useLongQueue} = 1;
my $treeFile = "$btout/phylo/IQtree_allsites.treefile";

my $cmd = "";
if (!-e $treeFile){
	print "Creating phylogeny for found specI's//\n";
	$cmd .= "$bts  -aa  $btout/all.faa -smplSep '\\$SaSe' -cats $btout/all.cats -outD $btout -runIQtree 1 -runFastTree 0 -runRaxMLng 0 -cores $numCores  -AAtree 1 -bootstrap 5000 -NTfiltCount 300 -NTfilt 0.1 -NTfiltPerGene 0.5 -minOverlapMSA 2 -MSAprogram $MSAprog -AutoModel 0 -iqFast 0 \n";
} else {
	print "Found already existing tree, skipping tree building\n";
	$cmd .= "#$bts  -aa  $btout/all.faa -smplSep '\\$SaSe' -cats $btout/all.cats -outD $btout -runIQtree 1 -runFastTree 0 -runRaxMLng 0 -cores $numCores  -AAtree 1 -bootstrap 5000 -NTfiltCount 300 -NTfilt 0.1 -NTfiltPerGene 0.5 -minOverlapMSA 2 -MSAprogram $MSAprog -AutoModel 0 -iqFast 0 \n";
}
$cmd .= "\n\n\n$xtraMessageInSH\n" if ($xtraMessageInSH ne "");

#add script for phylo visualization
$cmd .= "\n#visualize the newly created phylogeny\n";
my $abundMatrix = $MGSfile;  $abundMatrix =~ s/\/[^\/]+$/\//; $abundMatrix .= "Annotation/Abundance/MGS.matL7.txt";
$cmd .= "$vizTree $abundMatrix $treeFile $btout/phylo/IQtree_allsites.pdf \n";

#handle submission
my $scrNm = "btwFMGtree";
$scrNm = "btwCusFMGTree" if ($addRefGenos ne "");
my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
my ($dep,$qcmd) = qsubSystem($btout.$scrNm.".sh",$cmd,$numCores,int($mem/$numCores)."G",$scrNm,"","",1,[],$QSBoptHR);
#$cmd= "$bts  -aa  $btout/all.faa -smplSep '\\$SaSe' -cats $btout/all.cats -outD ${btout}_ST -runIQtree 1 -runFastTree 0 -runRaxMLng 0 -cores $numCores  -AAtree 1 -bootstrap 000 -NTfiltCount 300 -NTfilt 0.1 -NTfiltPerGene 0.5 -minOverlapMSA 2 -MSAprogram 2 -AutoModel 0 -iqFast 1 -superTree 1 \n";
#($dep,$qcmd) = qsubSystem($btout."btweenTreeST.sh",$cmd,$numCores,"1G","FMGstStree","","",1,[],$QSBoptHR);

if ($wait4job==1){
	qsubSystemJobAlive( [$dep],$QSBoptHR );
}
if ($wait4job==2){
	print "WAITID=$dep\n";
}
