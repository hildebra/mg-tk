#!/usr/bin/perl
#perl kmer_Ngenes.pl /g/scb/bork/hildebra/SNP/Drama1/AssmblGrp_16020939/metag/ContigStats/scaff.pergene.4kmer.gz 5
#this script takes existing kmer stats per gene and calculates over a window of X genes the average kmer stats

use Mods::GenoMetaAss qw(gzipopen);
use Mods::math qw(avgArray roundAr);

use warnings;
use strict;

my $inK = $ARGV[0];
my $numG = $ARGV[1];

my $outF = $inK;
$outF =~ s/4kmer\.gz/4kmer\.pm$numG\.gz/;
system "rm -f $outF";
open O,"| gzip -c > $outF" or die "Can't open $outF\n";
my ($I,$ok) = gzipopen($inK,"K-mer per gene");
die "can't open input kmer\n" if (!$ok);
my $cctg= "";
my @roll=();;my @roGenes=();;
while (my $line = <$I>){
	chomp $line;
	my @spl = split /\t/,$line;
	my $gne = shift(@spl);
	if ($gne eq "Contig"){
		print O $line."\n";
		next;
	}
	unless ($gne =~ m/(^.*)_\d+$/){
		die "can't find contig info $gne\n";
	}
	my $ctg = $1;
	if ($cctg ne $ctg){
		$cctg = $ctg;
		#release old stats
		my $rSize = scalar @roll;
		while ($rSize){ #roll back from the front
			my $aar = avgArray(\@roll); $aar = roundAr($aar,2); my @aa = @{$aar};
			my $idx = $rSize - $numG; $idx=0 if ($idx < 0);
			print O $roGenes[$idx]."\t".join("\t",@aa)."\n";
			shift @roll; shift @roGenes;
			$rSize = scalar @roll;
		}
	}
	#print "$ctg\n";
	push(@roll,\@spl);
	push (@roGenes,$gne);
	my $rSize = scalar @roll;
	if ($rSize >= $numG){
		my $aar = avgArray(\@roll);  $aar = roundAr($aar,2);my @aa = @{$aar};
		print O $roGenes[$rSize - $numG]."\t".join("\t",@aa)."\n";
		#and check if it needs prunning..
		if ($rSize > $numG*2){
			shift @roll; shift @roGenes;
		}
	}
	
}

close $I;
my $rSize = scalar @roll;
while ($rSize){ #roll back from the front
	my $aar = avgArray(\@roll);  $aar = roundAr($aar,2); my @aa = @{$aar};
	my $idx = $rSize - $numG; $idx=0 if ($idx < 0);
	print O $roGenes[$idx]."\t".join("\t",@aa)."\n";
	shift @roll; shift @roGenes;
	$rSize = scalar @roll;
}


close O;

#system "rm-f $outF.gz;gzip $outF;rm -f $outF";