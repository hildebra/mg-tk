#!/usr/bin/perl
#extracts FAA/FNA for each MGS
#perl extractMGSgenes.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Canopy4_AC/clusters.txt /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Canopy4_AC//extr/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ 200
#perl extractMGSgenes.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/extr/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ 200 /scratch/bork/hildebra/MGStest/



use strict;
use warnings;
use Mods::GenoMetaAss qw(systemW readFasta);
use Mods::Binning qw(runCheckM);



sub readCluster{
	my ($cF) = @_;
	my %canos;
	open I,"<$cF" or die "can;t open $cF\n";
	while (my $line = <I>){
		chomp $line;
		$line =~ s/ //g;
		my @spl = split /\t/,$line;
		push(@{$canos{$spl[0]}}, $spl[1]);
	}
	close I;
	return \%canos;
}
my $ncore = 20;


my $cluF = $ARGV[0];
my $oDir = $ARGV[1];
my $GCd = $ARGV[2];
my $tmpD = "$oDir/tmp/";
my $minGenes=200;
$tmpD = $ARGV[4] if (@ARGV >= 4);
$minGenes = $ARGV[3];
system "rm -r $oDir" if (-e $oDir);
system "mkdir -p $oDir" unless (-d $oDir);
system "mkdir -p $tmpD" unless (-d $tmpD);

my $hr = readCluster($cluF);
my %clust = %{$hr};

print "Reading ref FNA..\n";
$hr = readFasta("$GCd/compl.incompl.95.fna",1);
my %FNA = %{$hr};
#my @test = keys %FNA; print "$test[0] $test[1] $test[123]\n"; print "$FNA{13220655}\n";
foreach my $cl (sort keys %clust){
	my $oF = "$oDir/$cl.fna";
	my @refG = @{$clust{$cl}};
	next if (scalar(@refG) < $minGenes);
	open O,">$oF" or die $!;
	foreach my $rg (@refG){
		chomp $rg;
		die "Can't find gene $rg in gene cat\n" unless(exists($FNA{$rg}));
		my $rn = $rg; $rn =~ s/://;
		print O ">$rn\n$FNA{$rg}\n";
	}
	close O;
}


print "Reading ref FAA..\n";
$hr = readFasta("$GCd/compl.incompl.95.prot.faa",1);
my %FAA = %{$hr};
foreach my $cl (sort keys %clust){
	my $oF = "$oDir/$cl.faa";
	my @refG = @{$clust{$cl}};
	next if (scalar(@refG) < $minGenes);
	open O,">$oF" or die $!;
	foreach my $rg (@refG){
		die "Can't find gene $rg in gene cat\n" unless(exists($FAA{$rg}));
		my $rn = $rg; $rn =~ s/://;
		print O ">$rn\n$FAA{$rg}\n";
	}
	close O;
}

print "Starting checkm...\n";
my $outFile = $cluF.".cm";
runCheckM($oDir,$outFile,$tmpD,$ncore) unless (-e $outFile);



print "Done\n";



