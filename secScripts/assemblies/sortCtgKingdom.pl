#!/usr/bin/env perl
#takes a set of contigs and assigns them to eukaryote or prokaryote - needed for gene predictions
#./sortCtgKingdom.pl /g/bork3/home/hildebra/data/SNP/SimuB/AssmblGrp_a1/metag/scaffolds.fasta.filt /scratch/bork/hildebra/GNMass3/testSimu/ /scratch/bork/hildebra/GNMass3/ /scratch/bork/hildebra/GNMass3/ 12
use warnings;
use strict;
use Mods::GenoMetaAss qw(readFasta );
use Mods::IO_Tamoc_progs qw( getProgPaths );
sub readNogKingdom($);

#scaffolds, path to save files (tmp), cores to use
die "Not enough input args\n" if (@ARGV < 3);
my ($inScaff, $tmpPath, $krakenDBDirGlobal, $DiaDBdir, $numCore) = @ARGV;

my $krkBin = getProgPaths("kraken");#"/g/scb/bork/hildebra/DB/kraken/./kraken";
my $diaBin = getProgPaths("diamond");#"/g/bork5/hildebra/bin/diamond/./diamond";
my $sizSplitScr = getProgPaths("sizSplit_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/splitFNAbyLength.pl";
my $cmd = "";

#my $metaGD = `cat $refSmpl/assemblies/metag/assembly.txt`; chomp $metaGD;
#my $smplAss = "$metaGD/scaffolds.fasta.filt";
system "mkdir -p $tmpPath" unless (-d $tmpPath);


#----------------- FNAs
my $bactSeqs = "$tmpPath/bact.kraken.fasta";
my $bactDiaSeqs = "$tmpPath/bact.diamond.fasta";
my $krakUncl = "$tmpPath/unclassi.kraken.fasta";
my $EukDiaSeqs = "$tmpPath/euk.kraken.fasta";
if (-s $inScaff <= 500){
	system "touch $bactSeqs $EukDiaSeqs";
	exit(0);
}


#first use normal mini_kraken with confidence to get Bacterial contigs
my $globalKraTaxkDB = "minikraken_2015";
my $curDB = "$krakenDBDirGlobal/$globalKraTaxkDB";
$cmd .= "$krkBin --preload --threads $numCore --fasta-input  --db $curDB  --unclassified-out $krakUncl --classified-out $bactSeqs $inScaff > /dev/null \n";
my $thresh = "0.2";
system "$cmd\n";


#second use diamond mappings to COG to get idea where the rest is hitting best
my $DBpath = "/g/bork3/home/hildebra/DB/FUNCT/eggNOG10/";	my $refDB = "eggnog4.proteins.all.fa"; 
die "diamond DB \"$DiaDBdir$refDB.db.dmnd\" doesnt exist\n" unless (-e "$DiaDBdir$refDB.db.dmnd");
my $outF = "$tmpPath/diaCOG.tab";
#$DBcmd .= "cp $DBpath/NOG.members.tsv $DBpath/NOG.annotations.tsv $DBpath/all_species_data.txt $CLrefDBD\n";
$cmd = "";
$cmd = "$diaBin blastx -f tab  --quiet -t $tmpPath -d $DiaDBdir$refDB.db -k 1 --query-cover 40 -q $krakUncl -k 5 -e 1e-9 -o $outF -p $numCore\n";
system "$cmd";

#read in and determine tax of reads
my $NOGtaxf = "$DiaDBdir/all_species_data.txt";
my %NOGkingd = readNogKingdom($NOGtaxf);

open I,"<$outF" or die "$outF diamond result does not exist\n";
my %refQuerys; my %EukQuerys; my $bacCnt=0; my $EukCnt=0;
while (my $line = <I>){
	chomp $line;
	my ($Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc) = split /\t/,$line;
	next if ($eval > 1e-11 || $id < 70);
	$Subject =~ m/^(\d+)\./; my $taxid = $1;
	if (!exists $NOGkingd{$taxid}){print "Can't find $taxid in ref tax\n";}
	my $kgdm = $NOGkingd{$taxid};
	if ($kgdm==0 || $kgdm==1){
		unless (exists ($refQuerys{$Subject} ) ){
			$refQuerys{$Query} = 1; $bacCnt++;
		}
	} else {
		$EukCnt++;
	}
	#print "XX";
}
close I;

print "Bacs=$bacCnt, Euks=$EukCnt\n";

#write to file
my $hr =readFasta($krakUncl,1); my %fnas = %{$hr};
open O,">$bactDiaSeqs" or die "Can't open $bactDiaSeqs\n";
foreach my $hd (keys %refQuerys){
	if (!exists($fnas{$hd})){die "can't find bac fna for $hd\n";}
	print O ">$hd\n$fnas{$hd}\n";
}
close O;
open O,">$EukDiaSeqs" or die "Can't open $EukDiaSeqs\n";
foreach my $hd (keys %EukQuerys){
	if (!exists($fnas{$hd})){die "can't find euk fna for $hd\n";}
	print O ">$hd\n$fnas{$hd}\n";
}
close O;

system "cat $bactDiaSeqs >> $bactSeqs; rm -f $bactDiaSeqs";
#system "cat $EukDiaSeqs >> $bactSeqs; rm -f $EukDiaSeqs";













sub readNogKingdom($){
	my ($inT) = @_;
	my %ret;
	open I,"<$inT" or die "Can't open NOG tax file\n";
	while (my $line = <I>){
		next if ($line =~ m/^#/);
		my $cls = 0; #bacteria
		$cls = 1 if ($line =~ m/Archaea/);
		if ($line =~ m/Eukaryota/){
			$cls = 2;
			$cls =3 if ($line =~ m/Fungi/);
		}
		$line =~ m/^(\d+)\s/;
		$ret {$1} = $cls;
	}
	close I;
	return %ret;
}



