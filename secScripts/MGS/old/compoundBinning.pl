#!/usr/bin/env perl
#uses a multi sample assembly to bin contigs with metabat
#see MGS.pl for gene catalog version of this
# perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec//helpers/MGS/compoundBinning.pl /g/scb/bork/hildebra/SNP/GCs/alienGC2/ /local/bork/hildebra/MB2test/ /g/scb/bork/hildebra/SNP/GCs/alienGC2/Canopy2/clusters.txt
# perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec//helpers/MGS/compoundBinning.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ /local/bork/hildebra/MB2test/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Canopy4_AC/clusters.txt
use warnings;
use strict;
use Data::Dumper;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(readMap  readClstrRev unzipFileARezip getAssemblPath systemW gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive);
use Mods::IO_Tamoc_progs qw(jgi_depth_cmd );
use Mods::TamocFunc qw ( getFileStr);
use Mods::Binning qw (runMetaBat createBinFAA runCheckM readMGS );

sub runMetaBat;
sub MB2assigns;
sub getGoodMBstats;
sub printL;
sub countUpBin;
sub filterMB2;
sub getConvergeWcanopy;
sub Rhclusts; sub evalRhcl_Bin;
sub invertIndex;

my $mb2Qual = getProgPaths("mb2qualCheck_scr");
my $rareBin = getProgPaths("rare");


#add this? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4748697/figure/fig-1/

if (@ARGV == 0){die "no input args given!\n";}

my $inD = $ARGV[0];
$inD.="/" unless($inD =~ m/\/$/);
my $doSubmit = 1;
my $numCore = 4;
my $canCore = 20;
my $treeCores=30;

#set up basic structures
my $QSBoptHR = emptyQsubOpt($doSubmit,"");
my %QSBopt = %{$QSBoptHR};
my $tmpD = $inD."/tmp/";
$tmpD = $ARGV[1] if (@ARGV > 1);
my $outD = $inD."/Binning/MetaBat/";
my $logDir = $outD."LOGandSUB/";
my $singleSample = 1;
print "Single Sample MetaBatting..\n";

system "mkdir -p $tmpD $outD $logDir";
#single sample
my $SmplNm = `cat $inD/mapping/done.sto`;
$SmplNm =~ s/-smd.bam\n?//;
my $jgiCov = $inD."mapping/$SmplNm-smd.bam.jgi.cov";my $jgiConn = $inD."mapping/$SmplNm-smd.bam.jgi.pairs.sparse";
my ($befZ,$aftZ) = unzipFileARezip( [$jgiCov,$jgiConn] );
$aftZ .= "\nrm -r $tmpD\n";	my $metaGD = getAssemblPath($inD);
my $scaffs = $metaGD."/scaffolds.fasta.filt";
$befZ .= jgi_depth_cmd($inD,$tmpD."covj",95,$numCore, $scaffs);
my $MBcmd =  runMetaBat($tmpD."covj.jgi.depth.txt",$outD,$SmplNm,$scaffs);
my $MBout = "$outD/$SmplNm";
systemW $befZ . $MBcmd . $aftZ . "\n" unless (-e $MBout && -e "$outD/MeBa.sto");	print "Single sample MetaBat2 done\n";

#qual check

my $postCmd = "";
$postCmd = "$mb2Qual $scaffs $MBout $tmpD\n" unless (-e "$MBout.cm" && -e "$MBout.assStat");
if ($MBcmd eq "" && $postCmd eq "") {exit(0);}#print "next "; next;}
print "Running checkM quality checks..\n";
systemW $postCmd;

exit (0);



