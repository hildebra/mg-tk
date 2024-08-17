#!/usr/bin/env perl
#combines the depth estimates of MATAFILER in a mapping dir
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipopen readMapS);
use Mods::IO_Tamoc_progs qw(getProgPaths);



my $inD = $ARGV[0];#"/g/scb/bork/hildebra/SNP/HMP_moSim/GlbMap/640069301/"; #testDir
$inD .= "/" unless ($inD =~ m/\/$/);
opendir D,$inD or die "Can't open dir $inD\n";
my @files = readdir D;closedir D;
my @tmp; foreach (@files){push @tmp,$_ if (m/.*smd.bam.coverage.gz.window$/);}

print "Found ".@tmp." valid depth wins\n";

#print "@tmp\n";
my %mat;my @smpl; my $genoN = "";
my $ctg = ""; my %ctgL;
foreach my $win (@tmp){
	$win =~ m/([^_]+)_(.*)-\d+-smd.bam.coverage.gz.window/;
	my $csmpl = $2;
	if ($genoN ne ""){die "two different geno names $genoN / $1\n" if ($1 ne $genoN);}
	$genoN = $1;
	push(@smpl,$csmpl);

	#die "$inD$win\n";
	open I,"<$inD$win" or die "Can't open $inD$win\n" or die "Can't find $inD$win\n";
	while (my $line = <I>){
		chomp $line;
		my @spl = split /\t/,$line;
		if ($line =~ m/^\D/){
			$ctg = $spl[0];
			if (exists($ctgL{$ctg})){
				$ctgL{$ctg} = $spl[1] if ($ctgL{$ctg} < $spl[1]);
			} else {
				$ctgL{$ctg} = $spl[1];
			}
			next;
		}
		#print "@spl\n";
		$mat{$ctg}{$spl[0]}{$csmpl} = $spl[1];
	}
}

my $goCnt=0; my $misCnt=0;
my $outF = $inD."$genoN.all.coverage.gz.window";
@smpl = sort(@smpl);
open O,">$outF" or die "Can't open matrix out $outF\n";
print O "#Position";
foreach my $sa (@smpl){
	print O "\t$sa";
}
print O "\n";
foreach my $ctg (keys %mat){
	my @siz = sort {$a <=> $b} (keys %{$mat{$ctg}});
	print "$ctg\t";
	print O "$ctg\t$ctgL{$ctg}\n";
	foreach my $si (@siz){
		print O $si;
		foreach my $sa (@smpl){
			if (exists($mat{$ctg}{$si}{$sa})){
				print O "\t$mat{$ctg}{$si}{$sa}";
				$goCnt++;
			} else {
				print O "\t0";
				$misCnt++;
			}
		}
		print O "\n";
	}
}

print "$goCnt entries and $misCnt missed entried\n$outF\n";