#!/usr/bin/env perl

#create refDB specific to each metag sample
#script that maps refTarDB against assembly of respective sample, removing contigs that are already in refDB
#./deployMapDB.pl /g/bork3/home/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa /g/bork3/home/hildebra/data/SNP/GNMass3/alien-14-392-3/ /scratch/bork/hildebra/GNMass3/alien-11-376-0/specDB/T2combDB.fna 12

use strict; use warnings;use threads;
use Mods::IO_Tamoc_progs qw(getProgPaths buildMapperIdx);

if (@ARGV < 4){die "not enough args to function \n";}
my ($refTar, $refSmpl, $tmpdb, $ncore,  $basns, $fdirs) = @ARGV;
$tmpdb =~ m/^(.*\/)[^\/]*$/;
my $tmpd = $1;
#die $tmpd."\n";
system "mkdir -p $tmpd";#tmpd is also outD;
my $metaGD = `cat $refSmpl/assemblies/metag/assembly.txt`; chomp $metaGD;
my $smplAss = "$metaGD/scaffolds.fasta.filt";
die "Can't find assembly of sample at\n$smplAss\n" unless (-e "$metaGD/scaffolds.fasta.filt");

my $blatBin = getProgPaths("blat");
my $pigzBin = getProgPaths("pigz");
#$cmd .= "$blatBin -makeOoc=$refTar.11.ooc $refTar\n";
#die $cmd;
my @refs = split /,/,$refTar; my @finalD = split /,/,$fdirs; my @baseN = split /,/,$basns;
die "unequal array lengths: ".@refs." ".@finalD." ".@baseN."\n" if (@refs != @finalD || @refs != @baseN ); 
my @thrs; my @b8s;
for (my $i=0;$i<@refs;$i++){
	my $b8file = "$tmpdb.$i.b8";my $cmd ="";
	$cmd .= "$blatBin -t=dna -q=dna -minIdentity=95 -minScore=100 -out=blast8 $refs[$i] $smplAss $b8file\n";
	push (@b8s, $b8file);
	
	$cmd .= "$pigzBin -p $ncore -c $b8file > $finalD[$i]/$baseN[$i].b8.gz\n";
	#die $cmd;
	$thrs[$i] = threads->create( sub {system $cmd; } );
}
print "@b8s\n\n";
for (my $i=0;$i<@refs;$i++){
	$thrs[$i]->join();
}
my $b8file = "$tmpdb.all.b8";
#concatenate and sort (to get all mappings to single contig)

if (@b8s == 1){
	system "mv ". join(" ",@b8s) . " $b8file";
} else {
	system "cat ". join(" ",@b8s) . " | sort > $b8file";
	#b8s are important info for later steps (inversions etc)
	system "rm ". join(" ",@b8s);
}
#$cmd =  "/g/bork3/home/hildebra/bin/MUMmer3.23/nucmer -o $tmpdb.coords $refTar $smplAss";
#print $cmd;
#system $cmd;

#die;
open I,"<$b8file";
my $curSeq = ""; my $curL = 0;my $curCov=0;my $quAcc=0;
my %excl;
while (<I>){
	chomp;
	my @spl = split /\t/;
	if ($spl[0] ne $curSeq){
		$curSeq = $spl[0]; $curSeq =~ m/_L=(\d+)[;=]/; $curL = $1;$curCov=0;$quAcc=0;
	}
	next if ($quAcc==1);
	if ($spl[2] > 95){$curCov += $spl[3];}
#	if ($curSeq eq "MM33M9__C803_L=58531;"){ print $curCov."\n";}
	if ($curCov >= $curL*0.8){#accept hit (rm from db)
		$quAcc=1; $excl{$curSeq} = 1;
	}
}
close I;
#die;
#system "rm $b8file";
#MM4__C18_L=178675; MM4__C2_L=366023; MM4__C885_L=2840; MM4__C46_L=92476; MM4__C8395_L=403; MM4__C4833_L=601; MM4__C5_L=285576; MM4__C2681_L=944; MM4__C733_L=3783; MM4__C25_L=157425; MM4__C8_L=253037; MM4__C277_L=14651; MM4__C1998_L=1207; MM4__C293_L=13658; MM4__C916_L=2723; MM4__C4663_L=617; MM4__C736_L=3772; MM4__C789_L=3380; MM4__C30_L=122477; MM4__C0_L=694123; MM4__C608_L=4959; MM4__C20_L=167037; MM4__C94_L=47848; MM4__C95_L=47401; MM4__C5363_L=560; MM4__C667_L=4339; MM4__C33_L=117351;

#my @exCtgs = keys %excl;
#die "@exCtgs\n";
#create fasta without set
my $newRefDB = "$tmpdb";
open I,"<$smplAss" or die "Can't open assembly $smplAss\n";;
open O,">$newRefDB";
my $skip =0; my $skipCnt=0; my $skipBpCnt=0;
while (<I>){
	if (m/^>(\S+)$/){
		$skip=0;
		if (exists($excl{$1})){
			#print "SK: $1\n";
			$skip=1; $skipCnt++;
		}
	}
	if ($skip){
		$skipBpCnt += length($_);
		next;
	}
	print O $_;
}
close O; close I;
#and attach tar to new ref
foreach my $refFNA (@refs){
	system "cat $refFNA >> $newRefDB";
}
print "Skipped $skipCnt fasta entries, $skipBpCnt bps, added $basns fasta(s).\n";

#and build bowtie2 DB on $newRefDB
my ($bcmd,$newRef) = buildMapperIdx($newRefDB,$ncore,0,0);

system $bcmd;
exit(0);
#return "$newRef";


