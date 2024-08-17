#!/usr/bin/env perl
#actually changed to a very specific script to map CAZy onto eggNOG!!!
#script to just do funct annotation of fasta file (NTs)
#perl FuncAnnoFasta.pl /hpc-home/hildebra/DB/Funct/eggNOG10/eggnog4.proteins.all.fa CZy 24 

use warnings;
use strict;
use Mods::FuncTools qw(assignFuncPerGene calc_modules);
use Mods::FuncTools qw(readGene2COG);
use Mods::Subm qw(qsubSystem emptyQsubOpt );
use Mods::GenoMetaAss qw(gzipwrite gzipopen );

die "Not enough args\n" if (@ARGV < 1);

my $query = $ARGV[0];
my $DB = $ARGV[1];
my $ncore = $ARGV[2];
$query =~ m/^(.*\/)[^\/]+$/;
my $inP = $1;
my ($tmpD, $doClean,$fastaSplits) = ($inP,1,0);
my $outD = $inP."/Anno/";
my $QSBoptHR = emptyQsubOpt(1,"");
$QSBoptHR->{qsubDir} = $outD;
	
if (0){
	#my $DB = "NOG";
	#my $ncore = 40; 
	
	
	my $curDB = $DB;#"NOG";#CZy,ABRc,KGM,NOG
	
	my %optsDia = (eval=>1e-20,percID=>90,minPercSbjCov=>0.8,fastaSplits => $fastaSplits,ncore=>$ncore,
			splitPath=>$tmpD,keepSplits=>!$doClean,redo=>$doClean, minAlignLen=>30, minBitScore=>45);
			
			
	my ($allAss,$jdep) = assignFuncPerGene($query,$outD,$tmpD,$curDB,\%optsDia,$QSBoptHR,"diamond");
	#my $tarAnno = "${allAss}geneAss.gz";
	#my $tmpP2 = "$tmpD/CNT_1e-8_25//";
	#create actual COG table
	#$tarAnno =~ s/\.gz$//;
	print "Done\n";
}
#read ref DB
my $DButil = "/hpc-home/hildebra/DB/Funct/eggNOG10/";
my $bl2dbF = "$DButil/NOG.members.tsv";
my ($hr1,$hr2) = readGene2COG($bl2dbF);
my %g2COG = %{$hr1}; my %c2CAT = %{$hr2};

print "Reading annotation\n";

my $tarAnno = "$outD/DIAass_CZy.srt.gzgeneAss.gz";
my ($I,$OK) = gzipopen($tarAnno,"assignments");
my %COG2GZ;
open O,">$outD/COGassigns.txt";
while (my $line = <$I>){
	chomp $line;
	my @spl = split /\t/,$line;
	#die "@spl\n";
	if (exists($g2COG{$spl[0]})){
		if (@spl>2){
			print O "$g2COG{$spl[0]}\t$spl[1]\t$spl[2]\n";
		} else {
			print O "$g2COG{$spl[0]}\t$spl[1]\t\n";
		}
		$COG2GZ{$g2COG{$spl[0]}} = $spl[1];
	} else { print "can't find $spl[0]\n";}
}
close $I;
close O;

print "Done reading annotation\n";

print "Rewriting NOG tables\n";
my $Ndir = "/ei/workarea/users/hildebra/projects/keyTaxa/";
my @subNs = ("NOG.mat.ARC.cnt.1e-9_50.txt.gz","NOG.mat.BAC.cnt.1e-9_50.txt.gz","NOG.mat.EUK.cnt.1e-9_50.txt.gz","NOG.mat.FNG.cnt.1e-9_50.txt.gz");
foreach my $subN (@subNs){
	print "$subN\n";
	my %omat; my $Cfnd=0;
	($I,$OK) = gzipopen($Ndir.$subN,"assignments");
	my @head; my$cnt=0;
	while (my $line = <$I>){
		chomp $line;
		$cnt++;
		my @spl = split /\t/,$line;
		my $rID = shift @spl;
		if ($cnt==1){
			@head = @spl; next;
		}
		#use only genes that are in COG 2 CAZy assingments
		next unless (exists($COG2GZ{$rID}));
		$Cfnd++;
		for (my $i=0; $i< @spl;$i++){
			$spl[$i]=0 if ($spl[$i] eq "-");
		}
		my $newCzy = $COG2GZ{$rID};
		if (exists($omat{$newCzy})){#add
			for (my $i=0; $i< @spl;$i++){
				#next if ($spl[$i] eq "-");
				${$omat{$newCzy}}[$i] += $spl[$i];
			}
		} else {
			$omat{$newCzy} = \@spl;
		}
	}
	close $I;
	print "Found $Cfnd COGs\n";
	my ($OE,$OK2) = gzipwrite($Ndir."CZy.".$subN,"CZy out mat",1);
	print $OE "\t".join("\t",@head)."\n";
	foreach my $CZ (keys %omat){
		print $OE $CZ."\t".join("\t",@{$omat{$CZ}})."\n";
	}
	close $OE;
}


print "Done rewriting NOGs -> CZy\n";


