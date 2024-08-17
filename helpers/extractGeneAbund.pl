#!/usr/bin/perl
#./extractGeneAbund.pl /g/scb/bork/hildebra/SNP/GNMass2_singl/GeneCatalog/betaLacta.lst /g/scb/bork/hildebra/SNP/GNMass2_singl/ /g/scb/bork/hildebra/SNP/GNMass2_singl/GeneCatalog/betalact/
use warnings;
use strict;

sub readGeneIdx;

my $lstFile = $ARGV[0];
my $inD = $ARGV[1];
my $outD = $ARGV[2];

system "mkdir -p $outD" unless (-d $outD);

my $GCd = "$inD/GeneCatalog";
my $GCMversion = 0.1; if (-e "$GCd/version.txt"){$GCMversion = `cat $GCd/version.txt`; chomp $GCMversion;}
my $idxF = "$GCd/Matrix.genes2rows.txt";
#my $hr = readGeneIdx($idxF);
#my %gene2row = %{$hr};

#my $listl = `cat $lstFile`;
#my @list = split(/\n/,$listl);
#my @geneIdx;
#my %liHa; $liHa{$_}++ for (@list);
#foreach my $item (keys (%liHa)){
#	unless(exists($gene2row{$item})){die "can;t find gene $item\n";}
#	push(@geneIdx,$gene2row{$item});
#}
#add +2 as 1 is head and 0 is line 1 in mat
#if ($GCMversion <= 0.1){
#	foreach my $ele (@geneIdx){$ele+=2;} 
#}
#my $sedStr = join("p;",@list);
#system "sed -n '1p;$sedStr"."p' $GCd/Matrix.mat > $outD/Abund.mat";
system "/g/bork5/hildebra/dev/C++/rare/rare lineExtr $GCd/Matrix.mat $outD/Abund.mat $lstFile";

print "$outD/Abund.mat\n";





sub readGeneIdx(){
	my ($in) = @_;
	my %ret;
	open I,"<$in";
	while(my $line=<I>){
		chomp $line;
		my @spl = split(/\t/,$line);
		$ret{$spl[2]} = $spl[0];
		#print $spl[2] ." ". $spl[0]."\n";
	}
	close I;
	print "Gene Index read\n";
	return (\%ret);
}