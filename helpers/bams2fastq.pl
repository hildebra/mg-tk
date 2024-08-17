#!/usr/bin/perl
#converts in input dirs bams to fastq.. better to do this once in input raw dir
#redo: EGAF00001151487 EGAF00001151489 EGAF00001150355
use warnings; use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( systemW);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
my $smtBin = getProgPaths("samtools");

my $inDir = $ARGV[0];
die "Input dir $inDir is not a dir\n" unless (-d $inDir);
my $logdir = "$inDir/LOGandSUB/";
system "mkdir -p $logdir" unless (-d $logdir);

my @bams = grep { -e } glob "$inDir/*.bam";
@bams = grep { -e } glob "$inDir/*/*.bam" if (@bams == 0);

#die "@bams\n";

if (scalar @bams == 0){die "Could not find any bams!\n";}

my $QSBoptHR = emptyQsubOpt(1);

my $cnt =0;
my $ncore = 3;
#die "@bams";
foreach my $bb (@bams){
	my $o1 = $bb;my $o2 = $bb;  my $os = $bb;
	$o1 =~ s/\.bam$/\.1\.fastq\.gz/;
	$o2 =~ s/\.bam$/\.2\.fastq\.gz/;
	$os =~ s/\.bam$/\.singl\.fastq\.gz/;
	if ( $o1 =~ m/\.bam/ || $o2 =~ m/\.bam/ || $os =~ m/\.bam/ ){die "$o1\n$o2\n$os\n";}
	my $cmd = "$smtBin fastq -c 6 -@ $ncore -1 $o1 -2 $o2 -0 $os $bb\n";
	$cmd .= "rm $bb\n";
	#print "$cmd\n";
	qsubSystem($logdir."/bam2fq_$cnt.sh",$cmd,$ncore,"2G","b2f_$cnt","","",1,[],$QSBoptHR) ;
	$cnt++;
}
