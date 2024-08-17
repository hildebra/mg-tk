#!/usr/bin/perl
#script that evaluates MB2 qual within single assembly via checkm, gets assembly stats per Bin and makes assembly-bins ready for collection via compoundBinning.pl
#perl extractMeBat2.pl /g/bork3/home/hildebra/data/SNP/GNMass3_a/alien-11-883-0/assemblies/metag/scaffolds.fasta.filt /g/bork3/home/hildebra/data/SNP/GNMass3_a/alien-11-883-0/assemblies/metag/Binning/MB2/MM20 /tmp/XX
use warnings;
use strict;
use Mods::GenoMetaAss qw(  systemW readFasta);
use Mods::Binning qw(runCheckM runCheckM2);
use Mods::math qw(medianArray);
sub MB2assigns; sub MB2N50;
sub createBinFAA;


my $refFA = $ARGV[0];
my $MB2 = $ARGV[1];
my $tmpD = $ARGV[2];
my $MB2Dir = $MB2; $MB2Dir =~ s/[^\/]+$//;
my $version= 0.1;
my $ncore = 1;
$ncore = $ARGV[3] if (@ARGV > 3);
my $usCheckM2 = 0;
$usCheckM2 = $ARGV[4] if (@ARGV > 4);
my $usCheckM1 = 1;
$usCheckM1 = $ARGV[5] if (@ARGV > 5);
my $BinnerChoice = 0;
$BinnerChoice = $ARGV[6] if (@ARGV > 6);
$tmpD =~ s/\$/\\\$/g;
#die "$ncore\n";
#die "$tmpD\n";
print "Bin postprocessing v$version\n";
print "Extracting Bins, Qual check, reformatting into $MB2\n";

system "mkdir -p $tmpD" unless (-d $tmpD);
my $binD = "$tmpD/bins/";
system "mkdir -p $binD" unless (-d $binD);



#$isSemiBin = 1 if (-d "$MB2Dir/output_recluster_bins"); #!-e $MB2 &&  #likely SemiBin outdir 

my $emptyBin=0; #anything to do here, or no Bins found?
my %MB;

if ($BinnerChoice == 2 && !-s $MB2){
	#needs to create metabat like file..
	print "Detected SemiBat output dir..\n";
	#first prepare to delete all unused files..
	system "mv $MB2Dir/* $tmpD/";
	system "cp $tmpD/*.stone $MB2Dir";
	if (-d "$tmpD/output_recluster_bins"){
		system "mv $tmpD/output_recluster_bins/* $binD" ;
	} elsif (-d "$tmpD/output_bins"){
		system "mv $tmpD/output_bins/* $binD" ;
	}
	
	my $repStr = "";
	opendir(DIR, $binD) or die "Could not open $binD\n";
	while (my $filename = readdir(DIR)) {
		#print "$binD/$filename\n";
		
		my $binN = $filename; $binN =~ s/\.fa$//;
		next unless (-f "$binD/$filename" && -s "$binD/$filename");
		my $FR = readFasta("$binD/$filename");
		my %FAS = %{$FR};
		foreach my $c (keys %FAS){
			$repStr .= "$c\t$binN\n";
		}
	}
	closedir(DIR);
	#write actual output, if something to report..
	if ($repStr ne ""){
		open O,">$MB2" or die $!;
		print O $repStr;
		close O;
	}elsif (!-e $MB2){
		system "touch $MB2";
	}
} else {#MetaBat2 processing..
	print "Detected MetaBat2 output dir..\n";
}

#read in contig to bin assignments
system "touch $MB2" unless (-e $MB2);
my $hr = MB2assigns($MB2);
%MB = %{$hr};

#standardized path to recreate bin groups
createBinFAA($binD);

my $outFile = $MB2.".cm";
my $outFile2 = $MB2.".cm2";

print "found ". int(keys %MB) ." Bins\n";

if ($usCheckM1){
	my $tmpD2 = $tmpD."/CM/";
	if ($emptyBin){#no bins in metag..
		system "touch $outFile";
	} else {
		runCheckM($binD,$outFile,$tmpD2,$ncore,1,"fna") unless (-e $outFile);
	}
}
if ($usCheckM2){
	my $tmpD2 = $tmpD."/CM2/";
	if ($emptyBin){
		#system "touch $outFile2";
		open O,">$outFile2";
		print O "Name\tCompleteness\tContamination\tCompleteness_Model_Used Translation_Table_Used\tAdditional_Notes\n";
		close O;

	} else {
		runCheckM2($binD,$outFile2,$tmpD2,$ncore,1,"fna") unless (-e $outFile2);
	}
}
#and get N50 etc vals for each contig
if (1 || !-e "$MB2.assStat"){
	my $hr = MB2N50(\%MB);
	my %asS = %{$hr};
	open O,">$MB2.assStat" or die $!;
	print O "MB2\ttotalL\tmeanL\tctgN\tN20\tN50\tN80\tG1k\tG10k\tG100k\tG1M\n";
	my @itKeys = qw(tL meanL cN N20 N50 N80 1K 10K 100K 1M);
	foreach my $ak (keys %asS){
		#print "$ak\n"; print "$it\n";
		print O "$ak";
		foreach my $it(@itKeys){print O "\t$asS{$ak}{$it}";} #
		print O "\n";
	}
	close O;
}



#DONE






sub MB2N50($){
	my ($hr) = @_;
	my %ret;
	my %M = %{$hr};
	my %sizes_to_shorthand = (1000     => '1K',
							  10000    => '10K',
							  100000   => '100K',
							  1000000  => '1M',
							  10000000 => '10M');
	foreach my $k(keys(%M)){
		my @mem = @{$M{$k}};
		my $totL=0; my $ctgs=0; my @lengs;
		#die @mem;
		foreach my $x (@mem){
			$x =~ m/_L=(\d+)=/;
			push(@lengs,$1);
			$totL+=$1; $ctgs++;
		}
		my $meanL = $totL/$ctgs;
		
		$ret{$k}{tL} = $totL; $ret{$k}{meanL} = $meanL;$ret{$k}{cN} = $ctgs;
		# find number of sequences above certain sizes
		foreach my $size (1000,10000,100000,1000000){
			$ret{$k}{$sizes_to_shorthand{$size}}=0;
			foreach my $l (@lengs){
				$ret{$k}{$sizes_to_shorthand{$size}} ++ if ($l>=$size);
			}
		}
		#and find N50
		@lengs = sort { $a <=> $b } @lengs;
		my $N20 = int ($totL *0.2); my $N50 = int ($totL *0.5);my $N80 = int ($totL *0.8);
		my $cumL=0;
		foreach my $l (@lengs){
			$cumL += $l;
			if (!exists($ret{$k}{N20}) && $cumL >= $N20){$ret{$k}{N20} = $l;}
			if (!exists($ret{$k}{N50}) && $cumL >= $N50){$ret{$k}{N50} = $l;}
			if (!exists($ret{$k}{N80}) && $cumL >= $N80){$ret{$k}{N80} = $l;}
		}

	}
	return \%ret;
}

sub createBinFAA($){
	my ($binD) = @_;
	$hr = readFasta($refFA);
	system "mkdir -p $binD";
	my %FAS = %{$hr};

	$emptyBin= 1 if (-e $MB2 && -s $MB2 == 0);

	foreach my $bin (keys %MB){
		my @ctgs = @{$MB{$bin}};
		open O,">$binD/$bin.fna" or die $!;
		foreach my $ctg (@ctgs){
			die "can't find contig $ctg\n" unless (exists $FAS{$ctg});
			print O ">$ctg\n$FAS{$ctg}\n";
		}
		close O;
	}
	undef %FAS ;
}

sub MB2assigns($){
	my ($inF) = @_;
	my %ret;
	open I,"<$inF" or die "Can't open maxbin2 output $inF\n";
	while (<I>){
		chomp; my @spl  = split /\t/;
		next if ($spl[1] eq "0");
		next if ($spl[0] eq "Sequence ID");
		push(@{$ret{$spl[1]}}, $spl[0]);
	}
	close I;
	return \%ret;
}
