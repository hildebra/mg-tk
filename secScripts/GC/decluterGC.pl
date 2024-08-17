#!/usr/bin/env perl
#perl decluterGC.pl /g/bork3/home/hildebra/data/SNP/GCs/alienGC2 /local/hildebra/GCali/ 40
#script to look for 90% clutering genes from GC and check their abundance patter for clear anti-correlation
#this would mean that this gene is probably the same, but was clustered

use warnings;
use strict;
use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystem2 qsubSystemJobAlive);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::geneCat qw(readGene2tax createGene2MGS sortFNA);

die "not enough input args!\n" if (@ARGV<3);
my $GCd = $ARGV[0];
my $tmpD = $ARGV[1];
my $numCor = $ARGV[2];
my $local = 1;
my $totMem = 200;
$local = $ARGV[3] if (@ARGV > 3);
$totMem = $ARGV[4] if (@ARGV > 4);
my $finalStone = "";
$finalStone = $ARGV[5] if (@ARGV > 5);
my $cdhitBin = getProgPaths("cdhit");
my $rareBin = getProgPaths("rare");
my $pigzBin = getProgPaths("pigz");
my $mms2 = getProgPaths("mmseqs2");

my $clMode = 2; #0=cdhit, 2=mmseqs2

my $c90d = "$GCd/decluter/";
system "mkdir -p $c90d" unless (-d $c90d);

my $primFasta = "$GCd/compl.incompl.95.fna";
my $primFastaAA = "$GCd/compl.incompl.95.prot.faa";
my $ofna = "$tmpD/tmp.fna";
my $cltsv = "$tmpD/clust.txt";
system "mkdir -p $tmpD" unless (-d $tmpD);
print "=============================================\nDeclutering of gene catalog\n=============================================\n";
print "Phase I: clustering of gene cat and estimation of co-exclusion\n";
#choose very agressive parameters
my $cmd = "";
#DEBUG 
my $stone = "$c90d/declut.stone";
$cmd .= "mv $GCd/Matrix.mat.gz $c90d/Mat.predecl.mat.gz\n" if (-e "$GCd/Matrix.mat.gz");
#$cmd .= sortFNA($GCd,"compl.incompl.95",1,$tmpD,$numCor);
if ($clMode == 0){
	$cmd .= $cdhitBin."-est -i $primFasta -o $ofna -n 9 -G 1 -A 150 -r 0 -aS 0.4 -mask NX -aL 0.15 -d 0 -c 0.94 -g 1 -T $numCor -M ".int(($totMem+30)*1024)."\n";
	$cmd .= "mv $ofna* $c90d\n";
	$cmd .= "$rareBin decluter -i $GCd/Matrix.mat.gz -reference $ofna.clstr -o $tmpD -t $numCor -gz \n";
}elsif ($clMode == 2){ #mmseqs2 clustering 
	die "can't find input $primFastaAA" unless (-e $primFastaAA);
	my $MMdb = "$tmpD/compl.incompl.95.prot.faa.mms2.db";
	$cmd .= "touch $MMdb.stone\n";
	$cmd .= "echo \"Starting protein clustering at 95%\"\n";
	$cmd .= "$mms2 createdb  $primFastaAA $MMdb -v 1\n" unless (-e "$primFastaAA.mms2.db");
#	$cmd .= "$mms2 linclust $primFastaAA.mms2.db  $ofna $tmpD --cov-mode 2 -c 0.4 --min-seq-id 0.95 --threads $numCor\n";
	$cmd .= "$mms2 cluster $MMdb  $tmpD/primClus $tmpD --cov-mode 2 -c 0.3 --min-seq-id 0.95 --threads $numCor -v 2\n" unless (-e $cltsv);
	$cmd .= "$mms2 createtsv $MMdb $MMdb  $tmpD/primClus $cltsv --threads $numCor -v 2\n" unless (-e $cltsv);
#	$cmd .= "sort $tmpD/clus.tsv > $cltsv\n";
	$cmd .= "echo \"Starting gene matrix declutering\"\n";
	$cmd .= "$rareBin decluter2 -i $c90d/Mat.predecl.mat.gz -reference $cltsv -o $tmpD -t $numCor  -gz -pval 1e-5 \n";
	$cmd .= "rm -f $MMdb*\n";
	$cmd .= "echo \"Finished declutering\"\n";

} else {
	my $vsBin = getProgPaths("vsearch");
	$cmd .= "$vsBin --cluster_fast $primFasta --iddef 0 --usersort --mincols 150 --maxseqlength 90000 --consout $ofna --id 0.9 --strand plus --threads $numCor --uc $ofna.uc\n";
	die "$cmd\n";
}
$cmd .= "mv $tmpD/mat.decl.mat.gz $GCd/Matrix.mat.gz\n";
$cmd .= "mv $tmpD/concat.list $c90d\n";
$cmd .= "$pigzBin -p $numCor -c $GCd/compl.incompl.95.fna.clstr > $c90d/compl.incompl.95.fna.clstr.gz\n";
$cmd .= "$pigzBin -p $numCor -c $GCd/compl.incompl.95.fna.clstr.idx > $c90d/compl.incompl.95.fna.clstr.idx.gz\n";
$cmd .= "rm $GCd/compl.incompl.95.fna.clstr*\n";
$cmd .="touch $stone";

	#die "$cmd\n";


if ($clMode == 0){
}elsif ($clMode == 2){ 
}
#die "$cmd\n";
my @jdep=();
if (!-e "$stone"){
	if ($local){
		systemW $cmd;
	} else { #useful if only one core used in main routine..
		my $QSBoptHR = emptyQsubOpt(1,"");
		#$QSBoptHR->{useLongQueue} = 1;
		push(@{$QSBoptHR->{constraint}}, "sse4");
		$QSBoptHR->{useLongQueue} = 1;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		my ($dep,$qcmd) = qsubSystem($c90d."decluter.sh",$cmd,$numCor,int($totMem/$numCor)."G","declut","","",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;

		print "waiting for Phase I.. execute script again after clustering is finished\n"; #exit(0);
		push (@jdep, $dep);
		qsubSystemJobAlive( \@jdep,$QSBoptHR );
	}
}
if (!-e $stone){
	die "Could not find decluter step1 stone: $stone\n";
}

$cmd = "";
if (-e "$primFastaAA.mms2.stone"){
	die "decluter didn't work!\n $c90d/decluter.sh";
}
system "rm -f $primFastaAA.mms2";
#print "$tmpD/mat.decl.mat.gz\n$c90d/concat.list\n";
if (!-e $stone && (-e "$tmpD/mat.decl.mat.gz" || !-e "$c90d/concat.list")){
	$cmd = "mv $tmpD/mat.decl.mat.gz $GCd/Matrix.mat.gz\n";
	$cmd .= "mv $tmpD/concat.list $c90d\n";
	$cmd .= "mv $GCd/compl.incompl.95.fna.clstr* $c90d\n$pigzBin -p $numCor $c90d/compl.incompl.95.fna.clstr\n";
	$cmd .="touch $stone";
	systemW $cmd;
}
if (!-e "$stone"){
	die "decluter II didn't work!\n $c90d/decluter.sh";
}

#$ofna.".clstr";

print "Phase I complete\nPhase II: Rewriting clstr.idx\n";
#now comes just a bunch of rewriting..
my ($hr1,$hr2) = readClstrRev("$c90d/compl.incompl.95.fna.clstr.idx",0);my %cl2gene = %{$hr2}; $hr1 = {};

open I,"<$c90d/concat.list" or die $!;
while (my $line = <I>){
	chomp $line;
	#print "$line\n";
	my @spl = split /,/,$line;
	next if (@spl <= 1);
	my $tar = $spl[0];
	for (my $i=0;$i< @spl; $i++){
		die "can't find gene $spl[$i] XX $line\n" unless (exists($cl2gene{$spl[$i]}));
		next if ($i==0);

		#print "$cl2gene{$tar}\n";		print "$cl2gene{$spl[$i]}\n";
		#add to base gene
		$cl2gene{$tar} .= ",".$cl2gene{$spl[$i]};
		#and empty base gene
		$cl2gene{$spl[$i]} = "";
		#print "$cl2gene{$tar}\n";		print "$cl2gene{$spl[$i]}\n";
	}
	#die;
}
close I;
my @genes= sort(keys %cl2gene);
print "Rewriting into $GCd/compl.incompl.95.fna.clstr.idx\n";
open O,">$GCd/compl.incompl.95.fna.clstr.idx" or die $!;
foreach my $gene (@genes){
	print O "$gene\t$cl2gene{$gene}\n";
}
close O;
#$cmd = "$pigzBin -p $numCor $c90d/compl.incompl.95.fna.clstr.idx\n" unless (-e "$c90d/compl.incompl.95.fna.clstr.idx.gz");
$cmd = "gzip $c90d/compl.incompl.95.fna.clstr.idx\n" unless (-e "$c90d/compl.incompl.95.fna.clstr.idx.gz");
systemW $cmd;


systemW "touch $finalStone" unless ($finalStone eq "");
systemW "rm -r $tmpD\n";


print "Phase II decluter completed\n";






















