#!/usr/bin/env perl
#script to blast a GC against other GCs and calculate the overlap between these
use strict;
use warnings;
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::Subm qw(qsubSystem emptyQsubOpt);
my $lambdaBin = getProgPaths("lambda");
my $lambdaIdxBin = getProgPaths("lambdaIdx"); #lambda ref DB search
my $mini2Bin = getProgPaths("minimap2");#"/g/bork3/home/hildebra/bin/minimap2-master/minimap2";
my $sortSepScr = getProgPaths("sortSepReLen_scr");#"perl $thisDir/secScripts/sepReadLength.pl";
my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
my $smtBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";

my $selfSearch=2;#set to -1 to deactivate

my $DB = "/g/bork3/home/hildebra/data/SNP/GCs/SoilCatv2b1/compl.incompl.95.fna";
$DB = "/g/bork3/home/hildebra/DB/GeneCats/Soil/soil.catalog.95nr.fna";
$DB = "/scratch/bork/hildebra/geneCat/soil.catalog.95nr.fna";
my @fnaRefs = ("/g/bork3/home/hildebra/DB/freeze11/freeze11.genes.representatives.fa","/g/bork3/home/hildebra/DB/GeneCats/Tara/Tara.fna","/g/bork3/home/hildebra/DB/GeneCats/IGC/1000RefGeneCat.fna","/scratch/bork/hildebra/geneCat/soil.catalog.95nr.fna");
my @fnRNms = ("frz11","tara","gut1000","soil");
my $odir = "/g/bork3/home/hildebra/DB/GeneCats/Soil/cmp2otherGC/";
my $tmpDir = "/local/hildebra/cmpGC/";
$tmpDir = "\$TMPDIR/hildebra/cmpGC/";
my @tdirs = ("9305","9350","20250","20205");
my $DBn = "Soil";
if ($selfSearch>=0){
	$DBn = $fnRNms[$selfSearch]; $DB = $fnaRefs[$selfSearch];
	splice @fnaRefs,$selfSearch,1; splice @fnRNms,$selfSearch,1;
}
my $numCor=24;
my $BlastCores = $numCor; my $doBwt2 = 2; my $tag = "mini";
my $bwtCore = $numCor;

my $QSBoptHR = emptyQsubOpt(1,"");
my %QSBopt = %{$QSBoptHR};

system "mkdir -p $odir" unless (-d $odir);
system "mkdir -p $tmpDir" unless (-d $tmpDir);
my $cnt = -1; my $DBdep = "";
foreach my $query (@fnaRefs){
	$cnt++;my $cmd="mkdir -p $tmpDir\n";
	my $tmptaxblastf = "$odir/${DBn}_$fnRNms[$cnt].m8";
	my $tmptaxsamf = "$odir/${DBn}_$fnRNms[$cnt].sam";
	my $bwt35Log = "$odir/${DBn}_$fnRNms[$cnt].bwt2.log";
	next if (-e $tmptaxsamf);
	if ($doBwt2==0){
		if (!-d $DB.".lambda/"  ) {
			print "Building LAMBDA index anew (may take up to an hour)..\n";
			my $cmdIdx = "$lambdaIdxBin -p blastn -t ".int($BlastCores/2)." -d $DB";
			#if (system($cmdIdx)){print("Lamdba ref DB build failed\n$cmdIdx\n",3);}
			my ($dep,$qcmd) = qsubSystem($odir."QsubCmpGeneCatDB.sh",$cmdIdx,$numCor,"10G","cmpGCDB","","",1,[],$QSBoptHR);
			$DBdep = $dep;
			#die;
		}
#		$cmd .= "$lambdaBin $defLopt -q $query -i $DB.lambda -o $taxblastf\n";

		$cmd .= "$lambdaBin -t $BlastCores -id 95 -v 1 -nm 1 -p blastn -e 1e-10 -q $query -i $DB.lambda -o $tmptaxblastf\n";
	}elsif($doBwt2==2){
		$cmd .= "$mini2Bin -x asm5 -I 16G -t $numCor -d $tmpDir/query.mmi $query\n";
		$cmd .= "$mini2Bin -ax asm5 -N 0 -Q -t $numCor  $tmpDir/query.mmi $DB | $smtBin view -F 0x04 > $tmptaxsamf\n" #| $smtBin view -F 0x04
	} else {
		#bowtie2 way
		my $bwtIdx = $DB.".bw2";
		unless (1|| -e $bwtIdx.".rev.2.bt2" && -e $bwtIdx.".4.bt2"&& -e $bwtIdx.".2.bt2"){
			print "Building bwt2 catalog..\n";
			my $cmdIdx = $bwt2Bin."-build --threads $bwtCore -q $DB $bwtIdx\n" ;
			#if (system($cmdIdx)){print("Bowtie2 ref DB build failed\n$cmdIdx\n",3);}
			die "build DB\n";
			my ($dep,$qcmd) = qsubSystem($odir."QsubCmpGeneCatDB.sh",$cmdIdx,$numCor,"12G","cmpGCDB","","",1,[],$QSBoptHR);
			$DBdep = $dep;
			print "Done\n";
			die;
		}
		#take care of long reads
		my $tmpFNA = "$tmpDir/bwt24GC_$cnt.fna";
		$cmd .= "mkdir -p $tmpDir\n";
		$cmd .= "rm -f $tmpFNA;cp $query $tmpFNA\n";
		$cmd .= "$sortSepScr 7000 $tmpFNA\n";
		#and bulk of reads
		$cmd .= $bwt2Bin." --local --norc --no-unal --no-hd --no-sq -p $bwtCore ";
		$cmd .= " -x $bwtIdx -f -U  $tmpFNA >> $tmptaxsamf 2> $bwt35Log\n";

		$cmd .= $bwt2Bin." --local --norc --no-unal --no-hd --no-sq -p 1 ";
		$cmd .= "-x $bwtIdx -f -U  $tmpFNA.long > $tmptaxsamf 2>> $bwt35Log\n";
		
		#fix missing newlines
		$cmd .= "rm -f $tmpFNA\n";
	}
	#die "$cmd\n";
	print "Aligning sample $fnRNms[$cnt]...\n";
	#die;
	#system $cmd;
	$QSBoptHR->{useLongQueue} = 0;
	my ($dep,$qcmd) = qsubSystem($odir."QsubCmpGeneCat_${tag}_${DBn}_$cnt.sh",$cmd,$numCor,"9G","cmpGC$cnt",$DBdep,"",1,[],$QSBoptHR);
	print "$qcmd\n";
	#last if ($cnt == 0);
	#last;
	print "Done with alignment\n";

}






