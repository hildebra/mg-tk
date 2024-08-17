#!/usr/bin/env perl
#does minimus2 assembly of several assemblies of the same genome from different runs (same patient though)
use warnings;
use strict;
use Mods::GenoMetaAss qw(readMap readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
use Mods::IO_Tamoc_progs qw(getProgPaths jgi_depth_cmd inputFmtSpades createGapFillopt  buildMapperIdx);
sub nucEstDbl;
my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/renameCtgs.pl";
my $nucmerBin = "/g/bork3/home/hildebra/bin/MUMmer3.23/nucmer";
my $mini2Bin = "/g/bork3/home/hildebra/bin/amos-3.1.0_broken/bin/minimus2";


my $runID = "T3"; #name of files
my $preSbst = "";
if (0){
	$preSbst = "MM3,MM12M2,MM35M5";#,MM4,MM6,MM5,MM11,MM19,MM15"; #TEC3, everything with > 100 genes
} else {
	$preSbst = "MM3,MM4,MM1";#,MM36,MM35M5,MM256M2,MM12M2,MM33M9,MM5";#TEC4, everything with > 100 genes
	 $runID = "T4";
}
my @subsets = split (",",$preSbst);
my $inP = "/g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/$runID/R_filt/contigs2";
my $odir = "$inP/mini2/"; #output dir
my $tdir = "$odir/tmp/"; #temp

system "mkdir -p $odir" unless (-d $odir); system "mkdir -p $tdir" unless (-d $tdir);

my $tsmpl = shift @subsets;
my $firstFNA = "$inP/$tsmpl.ctgs.fna";
my $combFNA = "$tdir/comb.seq";
my $cmd = "";
#hash based
#$cmd = "/g/software/bin/minimus-3.1.0 $odir$runID.afg $odir$runID";

#nucmer based
my $cnt = 0; my @iniCtgs; my @iniBP; my @iniDblCov;
foreach my $secFNA (@subsets){
	#my $secFNA = "$inP/MM7.ctgs.fna";
	$secFNA = "$inP/$secFNA.ctgs.fna";
	 my $curBP=0;
	$cmd="";
	$cmd .= "cat $firstFNA $secFNA > $combFNA\n";
	$cmd .= "/g/bork3/home/hildebra/bin/amos-3.1.0_broken/bin/toAmos -s $combFNA -o $tdir$runID.afg\n";
	my $seqCnt = `grep -c '^>' $firstFNA`; chomp $seqCnt;
	my $hr = readFasta($firstFNA); my %fna = %{$hr}; foreach (keys %fna){$curBP += length($fna{$_});}
	push(@iniCtgs,$seqCnt); push (@iniBP,$curBP);
	my $trat = nucEstDbl($firstFNA);
	push (@iniDblCov, $trat);

	if ($cnt ==0){
#		$iniCtgs = $seqCnt;		$iniBP = $curBP;
		#die "$iniBP $iniCtgs\n";
	}
	
	$cmd .= "$mini2Bin $tdir$runID -D REFCOUNT=$seqCnt -D MINID=98 -D OVERLAP=150 > $tdir$runID.mini2.log\n";
	system "$cmd";
	print $cmd."\n";
	$firstFNA = "$tdir/prevRnd.fna";
	system "cat $tdir$runID.fasta  > $firstFNA";
	system "$renameCtgScr $firstFNA $tsmpl";
	$hr = readFasta("$tdir$runID.singletons.seq");	%fna = %{$hr};
	open O,">>$firstFNA";
	foreach (keys %fna){if (m/${tsmpl}__/){print O ">$_\n$fna{$_}\n";}};
	close O;
	
	
	#die "$tlen $tcov $trat\n";
}
#die $cmd ."\n";
#print "first file: $iniBP $iniCtgs\n";
print "BPs\tContigs\n";
my @subsetsX = split (",",$preSbst);

for (my $i=0;$i<@iniCtgs;$i++){
	print "$subsetsX[$i]\t$iniBP[$i]\t$iniCtgs[$i]\t$iniDblCov[$i]\n";
}

my $finCtgs = `grep -c '^>' $firstFNA`; chomp $finCtgs;
my $finBP=0;
my $hr = readFasta($firstFNA); my %fna = %{$hr}; foreach (keys %fna){$finBP += length($fna{$_});}
my $trat = nucEstDbl($firstFNA);

print "Now:        $finBP $finCtgs $trat\n";

system "cp $firstFNA $odir/$runID.mini2.fna";

exit (0);





sub nucEstDbl(){
	#nucmer duplication check
	my ($inFNA) = @_;
	my $cmd =  "$nucmerBin  --coords --prefix $tdir/nuc $inFNA $inFNA   2>&1  > $tdir$runID.nuc.log";
	system $cmd;
	open I, "<$tdir/nuc.coords" or die "can't open $tdir/nuc.coords\n"; my $startCoords=0;
	my %ctgFill;
	while (my $l = <I>){
		if (!$startCoords){
			$startCoords=1 if ($l =~ m/^======/);
			next;
		}
		chomp $l;
		my @spl = split /\|/,$l; for (my $i=0;$i<@spl;$i++){$spl[$i] =~ s/^\s+//;$spl[$i] =~ s/\s+$//;}
		#die "$spl[3]\n$spl[3]+1\n@spl\n";
		if ($spl[3] > 98){
			my $ctID = $spl[4]; $ctID =~ s/\s.*$//;
			my $hitL = $spl[2]; $hitL =~ s/\s.*$//;
			$ctgFill{$ctID} += $hitL;}
	}
	close I;
	my $tcov=0; my $tlen =0;
	foreach my $hd (keys %ctgFill){
		$hd =~ m/L=(\d+)[;=]/; my $len = $1;
		#print "$hd\t$len\t$ctgFill{$hd}\n";
		$tlen += $len; $tcov += $ctgFill{$hd};
	}
	#last;
	my $trat =  $tcov / $tlen;
	return $trat;
}









