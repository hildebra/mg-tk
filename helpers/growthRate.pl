#!/usr/bin/env perl
#perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/growthRate.pl /g/bork5/hildebra/results/TEC2/v5/Genomes/T6/435591.27.fna /g/bork5/hildebra/results/TEC2/v5/Genomes/T6/435591.27.faa
#ARGS: ./grothRate.pl [high expressed FNA] [all FNA (only complete genes)] 

use warnings;
use strict;
use threads ('yield',
                 'stack_size' => 64*4096,
                 'exit' => 'threads_only',
                 'stringify');
				 
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( systemW readFasta writeFasta);
use Mods::TamocFunc qw( getE100);

my $growthBin = getProgPaths("growthP");
my $tmpDir = "/tmp/MATAFILER/GrPr/";
system "mkdir -p $tmpDir" unless (-d $tmpDir);

die "needs two arguments [fna] [faa]\n" unless (@ARGV == 2);
my ($inFNA,$inFAA) = @ARGV;
my $hr = readFasta($inFNA);my %FNA = %{$hr};
my $fnaChg=0;
foreach my $k (keys %FNA){
	my $tmp = uc ($FNA{$k});
	if ($tmp ne $FNA{$k}){
		$FNA{$k} = $tmp; $fnaChg=1;
	}
}
writeFasta(\%FNA ,$inFNA) if ($fnaChg);
my $inP = $inFNA; $inP =~ s/([^\/]+)$//;
my $baseF = $1; $baseF =~ s/\.[^\.]+$//;
my $oDess = $inP."e100/"; system "mkdir -p $oDess" unless (-d $oDess);

if ( !-e "$oDess/alle100.fna"){
	getE100($oDess,$inFAA,$inFNA);
	systemW("cat $oDess/*.fna >$oDess/alle100.fna") unless (-e "$oDess/alle100.fna");
}
my $cmd = "$growthBin -f $oDess/alle100.fna -o $inP$baseF.grPred -g $inFNA -c 0 -T 37 -s -S --tmp $tmpDir";
if (system ($cmd)){ die "Couldn't execute command\n$cmd\n";}

