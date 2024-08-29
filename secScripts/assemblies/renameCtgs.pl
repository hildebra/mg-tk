#!/usr/bin/perl
#./renameCtgs.pl /g/scb/bork/hildebra/SNP/GNMass/alien-11-2-0/assemblies/metag/scaffolds.fasta
use warnings;
use strict;

sub newHD;

my $inF = $ARGV[0] ;
my $tag = $ARGV[1];
my $tmpOut = $inF.".tmp";
my $transOut = $inF.".lnk";

if (-z $inF){print "Warning:: $inF is empty!!\n";}

open I,"<$inF" or die "Cant open $inF\n";
open O,">$tmpOut" or die "Cant open $tmpOut\n";
open O2,">$transOut"; my $ohd = ""; my $ntag = ""; my $seq = "";
my $cnt = 0; 
my $line = <I>;
my @circCtgs; 
chomp $line; $ohd = $line;
#$ntag = newHD($line,$seq);
$cnt++;
while ($line = <I>){
	chomp $line;
	if ($line =~ m/>/){
		#next if (length($seq) ==0);
		my $ntag = newHD($seq);
		$seq =~ s/(.{1,80})/$1\n/gs;
		print O "$ntag\n$seq"; 
		print O2 "$ntag\t$ohd\n";
		if ($ohd =~ m/>ctg[\d_x]+c$/){push(@circCtgs,$ntag);} #old metaMDBG
		if ($ohd =~ m/>ctg[\d_x]+ .* circular=yes$/){push(@circCtgs,$ntag);}
		$cnt++;
		$seq = "";
		chomp $line; $ohd = $line;
	} else {
		$seq .= $line;
	}
}
$seq =~ s/(.{1,80})/$1\n/gs;
$ntag = newHD("",$seq);
print O "$ntag\n$seq"; $cnt++;
print O2 "$ntag\t$ohd\n";
close I; close O; close O2;
if (@circCtgs){
	#write out a file listing these contigs..
	open O ,">$inF.circ"; print O join("\n",@circCtgs)."\n"; close O;
}

system("rm $inF; mv $tmpOut $inF");
print "Done renaming contigs\n";


sub newHD($){
	my ( $seq) = @_;
	my $LendTag = "="; #;
	#if ($line =~ m/^>.*_length_(\d+).*/){
	#	$ntag = ">$tag"."__C$cnt"."_L=$1$LendTag";
	#} else {
		$ntag = ">$tag"."__C$cnt"."_L=".length($seq).$LendTag;
	#}

}
