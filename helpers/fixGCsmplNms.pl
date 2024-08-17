#!/usr/bin/env perl
#specific fix for gene cat protein sample names in Soil MG rast samples
use warnings;
use strict;
use Mods::GenoMetaAss qw( readMapS  systemW readFastHD);
use Mods::IO_Tamoc_progs qw(getProgPaths);


my $maps = "/g/scb/bork/hildebra/data2/refData/Soil/PNAS_refs/other_soil.map";#,/g/scb/bork/hildebra/data2/refData/Soil/howe2014/iowa.map,/g/scb/bork/hildebra/data2/Soil_finland/soil_map.txt";
my $GCD = "/g/scb/bork/hildebra/SNP/GCs/SoilCatv2b2/";
my ($hr,$hr2) = readMapS($maps);
my %map = %{$hr};
my @samples = @{$map{opt}{smpl_order}};
my %smplGNms; #big object to save gene names in..
my $gnCnt=0;
foreach (my $i=0;$i<@samples;$i++){
	my $smpl = $samples[$i];
	my $SmplName = $map{$smpl}{SmplID};
	#die "\n\n$SmplName\t$smpl\n";
	print $smpl."\n";
	if ($map{$smpl}{ExcludeAssem} eq "1"){next;}
	my $dir2rd = $map{$smpl}{wrdir};
	my $tarfile = "$dir2rd/assemblies/metag/genePred/proteins.bac.shrtHD.faa";
	open I,"<$tarfile" or die "can t open $tarfile\n";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			chomp $l;
			print "double entry: $l  $samples[$smplGNms{$l}]\n$dir2rd\n" if (exists($smplGNms{$l}) && $smplGNms{$l} != $i );
			$smplGNms{$l} = $i;
		} else {
			$gnCnt++;
		}
	}
	close I;
	#last if ($i>40);
}
#die "test $gnCnt\n";
open I,"<","$GCD/compl.incompl.95.fna.clstr";
open O,">","$GCD/compl.incompl.95.fna.clstr.2";
my $gnCnt=0;
while (my $line = <I>){
	if ($line =~ /^>/){print O $line;next;}
	#chomp $line; 
	unless ($line =~ m/nt, >(.*_\d+)\.\.\..*/){print O $line;next;;}
	#die $line.$1."\n";
	my $x=$1;
	my $y = ">$1";
	unless (exists($smplGNms{$y})){print O $line;next;}
	$line =~ s/^(.*nt, >).*(_\d+\.\.\..*)$/$1$samples[$smplGNms{$y}]__${x}_L=100$2/;
	#die $line."\n$x\n";
	print O $line;
	$gnCnt++;
}
close I; close O;
print "Fixed $gnCnt names\n";