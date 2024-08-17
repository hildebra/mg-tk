#!/usr/bin/env perl
#usage: ./eval_assembly.pl [assembly.fa] [ref.fa] [num cores] [id_of_cmp]
use warnings;
use strict;
use Mods::GenoMetaAss qw(readFasta);
use Mods::IO_Tamoc_progs qw(getProgPaths);

my $blatBin = getProgPaths("blat");
my $pigzBin = getProgPaths("pigz");
my $nucmBin = "";#getProgPaths("nucmer");

if (@ARGV < 3){die "not enough args to function \n";}
my ($assemblFa) = $ARGV[0];
my ( $refFa) = $ARGV[1];
my ( $ncore) = $ARGV[2];
my $smplID="";
$smplID = $ARGV[3] if (@ARGV>3);
#die "$smplID";
my $saveRes = 0;
my $tmpdb = $assemblFa.".tmp";

my $redo=1;

my $b8file = "$tmpdb.b8";my $cmd ="";
unlink($b8file) if ($redo);
$cmd .= "$blatBin -t=dna -q=dna -minIdentity=95 -minScore=100 -out=blast8 $refFa $assemblFa $b8file\n";
system ($cmd) unless (-e $b8file);
my $cfile = $tmpdb.".coords";
$cmd = "$nucmBin -p $tmpdb -o $refFa $assemblFa\n";
#system $cmd unless (-e $cfile);;

#$cmd .= "$pigzBin -p $ncore -c $b8file > $finalD[$i]/$baseN[$i].b8.gz\n" if ($saveRes);

#read ref to get stats on this
my $hr = readFasta($refFa);
my %ref = %{$hr};
my $refL=0;my$refCtg=0;
foreach my $k (keys %ref){
	$refL += length($ref{$k});
	$refCtg++;
}
print "REFERENCE------------\n$refFa\n   Contigs = $refCtg\n   Length  = $refL\n";
my @covered = (0) x $refL;
#read tar to get stats on this
$hr = readFasta($assemblFa);
my %tar = %{$hr};
my $assL=0;my$assCtg=0;
foreach my $k (keys %tar){
	$assL += length($tar{$k});
	$assCtg++;
	#print $k." ".length($tar{$k})."\n";
}
print "ASSEMBLY------------\n$assemblFa\n   Contigs = $assCtg\n   Length  = $assL\n";

my $doCtgRep=0;
if ($smplID ne ""){
	$doCtgRep = 1;
	open O,">$assemblFa.$smplID.hitCtg.txt";
}


#old, not used any longer
print "$b8file\n";
open I,"<$b8file";
my $curSeq = ""; my $curL = 0;my $curCov=0;my $totCov=0; my $covSst=0;my $covSed=0;my $assErrs=0;
my $curPID=0; my $totPID=0; my %hitSb = (); my $curRds=1; my $circCov =0;my $chimeric=0;
my $tooMuchBp=0;
#stores if a hit is convincing enough to reject all other hits to this query
my %qOK;
my %excl;
while (<I>){
	chomp;
	my @spl = split /\t/;
	if ($spl[0] ne $curSeq){
		if (exists($qOK{$curSeq})){
			print O "M	$curSeq\n" if ($doCtgRep);
			#collect erros from prev contig, if successful matched at all
			$assErrs += ($curRds-1);
			foreach my $xk(keys %hitSb){
				if ($hitSb{$xk} == 3){
					$circCov++;
					$assErrs -- ;
				}
			}
			#print $curSeq."\t$assErrs\n";
		} elsif ($curCov > 150 && $curL-$curCov>150){
			print O "C	$curSeq\n" if ($doCtgRep);
			$chimeric++;$totCov += $curCov;
			$tooMuchBp += $curL-$curCov;
			#print "CHIM $curCov $curL\n";
		} else {#total fail
			#$tooMuchBp += $curL;
		}
	
		$curSeq = $spl[0]; $curSeq =~ m/_L=(\d+)[;=]/; $curL = $1;%hitSb = ();	
		$covSst=0;$covSed=0;$curCov=0;$curPID=0;$curRds=0;
	}
	if ($spl[2] < 98 || exists($qOK{$spl[0]}) ){next;}
	
	#die if (scalar(keys %qOK) > 10);
	$curRds++;
	my $qSt=$spl[6];my $qEn=$spl[7]; my $sSt = $spl[8]; my $sEd = $spl[9];
	if ($sSt<=10 || $sEd <= 10){$hitSb{$spl[1]} += 1;}#$covSst=1;} 
	if ($sSt>=($refL-10) || $sEd >= ($refL-10)){$hitSb{$spl[1]} += 2;}#$covSed=1;} 
	$curPID += $spl[3] *$spl[2];
	$curCov += $spl[3];#just how many bases are covered
	if ($curCov >= ($curL*0.995-10)){ #defines a non-chimeric contig
		$qOK{$spl[0]}=1;$totCov += $curCov;	$totPID += $curPID;
	}
	my $start = $sSt; $start = $sEd if ($sEd < $start);
	my $end = $sSt; $end = $sEd if ($sEd > $start);
	for (my $ix = $start; $ix < $end; $ix ++){
		$covered[$ix]++;
	}
	
	#print "$hitSb{$spl[1]}\n";
}
close I;
close O if ($doCtgRep);
print "Total Coverage: $totCov , missed: ". ($refL-$totCov) .", notInRef: ". ($assL -$totCov) . "\n";
my $rounded = 0; 
$rounded = sprintf "%.3f", $totPID/$totCov if ($totCov>0); 
print "Assembly errors: $assErrs\nAvg. PID: ".$rounded." in ".scalar(keys %qOK)." contigs\n";
print "Circular Coverage: $circCov; Chimeric contigs: $chimeric\n";


#eval missing parts
my $miss=0;my %retMiss; my $ab1k=0;my $ab10k=0;
for (my $i=0;$i<@covered;$i++){
	if ($covered[$i]==0){
		$miss ++ ;
	} elsif ($miss>0){
		$retMiss{$miss} =  ($i-$miss)."\t$i";
		$miss=0;
	}
}
my $toosmall=0;my $toosmallCnt=0;my $prCnt=0;
foreach my $k (sort {$b <=> $a} keys %retMiss){
	if ($k>1000 && $prCnt < 8){
		print $k."bp\t".$retMiss{$k}.";\t";
		$prCnt++;
		if (($prCnt)%3==0){print "\n";}
	} else {
		$toosmall += $k; $toosmallCnt ++;
	}
}
print"\n $toosmall bp in $toosmallCnt smaller misses\n";









exit(0);







#nucmer.. doesn't find inversions, don't use
#open I,"<$cfile"; my $lcnt=0;
#while (my $l = <I>){
#	chomp $l;	$lcnt++;next if ($lcnt <= 5);	my @spl = split /\s+\|?\s+/,$l; shift @spl;	die "@spl\n$spl[0],$spl[3]\n";	
#}
#close I;



exit(0);
