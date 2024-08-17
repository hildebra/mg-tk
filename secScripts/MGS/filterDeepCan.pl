#!/usr/bin/perl
#perl filterDeepCan.pl /hpc-home/hildebra/grp/Links/Chicken_project/data/GC_Chicken2/Binning/MetaBat2/MB2.ext.can.pear.corr /hpc-home/hildebra/grp/Links/Chicken_project/data/GC_Chicken2/Binning/MetaBat2/MB2.ext.can.spear.corr /hpc-home/hildebra/grp/Links/Chicken_project/data/GC_Chicken2/Binning/MetaBat2/MB2.clusters.core 0.35 0.25 /hpc-home/hildebra/grp/Links/Chicken_project/data/GC_Chicken2/Binning/MetaBat2/MB2.ext.can.mrg.corr
#script to filter  deep canopies to more reasonable number of genes/MGS

use strict; use warnings;
sub readCan;
use Mods::math qw(quantileArray);
use Mods::Binning qw(readMGS);

my $maxGenes = 25000;

die "not enough input args" if (@ARGV < 2);
my $pearI = $ARGV[0]; my $speaI = $ARGV[1]; my $MB2cF =  $ARGV[2];
my $pearCut = $ARGV[3];my $speaCut = $ARGV[4];my $outFile = $ARGV[5];
#unless (-e "$pearI.bkup"){
#	system "mv $pearI $pearI.bkup";
#	system "mv $speaI $speaI.bkup";
#}

my $Phref=readCan("$pearI"); 
my $Shref=readCan("$speaI"); 
my $MBhref=readCan("$MB2cF"); 
my %MB2 = %{$MBhref};


print "Read input\n";

my %Pc = %{$Phref};
my %Sc = %{$Shref};
my @uniq = keys %{{%Pc,%Sc}};
my $MB2s=0;my $MB2l=0;
my %Fc; #stores just MGS -> Genes
foreach my $cM (@uniq){
	my %Pc1 ;%Pc1 = %{$Pc{$cM}} if (exists($Pc{$cM}));
	my %Sc1 ;%Sc1 = %{$Sc{$cM}} if (exists($Sc{$cM}));
	my %MB1; %MB1 = %{$MB2{$cM}} if (exists($MB2{$cM}));
	my @cGenes = keys %{{%Pc1,%Sc1}};
	if (@cGenes < $maxGenes){
		foreach (@cGenes){$Fc{$cM}{$_} = 1;} 
		foreach (keys %MB1){if (exists($Fc{$cM}{$_})){next;}$Fc{$cM}{$_} = 1;$MB2s++;} #just add the MB2 core genes.. in case..
		next;
	}
	#get all out that fit both Sc and Pc
	my @Pv = values %Pc1;	my @Sv = values %Sc1;
	my $nP = scalar(@Pv);	my $nS = scalar(@Sv);
	$nP = 1 if ($nP==0); $nS = 1 if ($nS == 0);
	#get median to get a grip of deviation..
	#my $medP = quantileArray(0.5,@Pv);	my $medS = quantileArray(0.5,@Sv);
	#print "$cM\nmedianS: $medS , medianP = $medP\n";
	my $Pperc = scalar(@Pv)/(scalar(@Pv) + scalar(@Sv));
	my $Sperc = 1-$Pperc;
	$Pperc = 10000/$maxGenes if ($Pperc*$maxGenes < 10000); 
	$Sperc = 10000/$maxGenes if ($Sperc*$maxGenes < 10000);
	
	#my $Pbr = ($maxGenes/scalar(@Pv)); my $Sbr = ($maxGenes/scalar(@Sv)); #breakpoint for quantile..
	#my $Pv8 = $pearCut*0.8;#
	#my $Sv8 = $speaCut*0.8;
	#$Pv8 = quantileArray($Pbr,@Pv) if ($Pbr<1);$Sv8 = quantileArray($Sbr,@Sv) if ($Sbr<1);
	
	#calc more strict cutoffs..
	my $Sv5 = $speaCut*0.8;	my $Pv5 = $pearCut*0.8;	
	my $Pbr = ($maxGenes*$Pperc/$nP); my $Sbr = ($maxGenes*$Sperc/$nS); #strict breakpoint for quantile..
	$Sv5 = quantileArray($Sbr,@Sv) if ($Sbr<1);$Pv5 = quantileArray($Pbr,@Pv) if ($Pbr<1);
	#print "valAccpets: $Sv8, $Pv8\n";
	foreach (@cGenes){
		my $eP = exists($Pc1{$_}); my $eS = exists($Sc1{$_});
		if ( ($eP && $eS ) || ( $eP && $Pc1{$_} <= $Pv5 || $eS && $Sc1{$_} <= $Sv5)	){ 
			$Fc{$cM}{$_} = 1;#condition hard: either good corr, and present in both sets
		}
	} 
	
	#ok in case still too extreme, rework again..
	if (scalar(keys(%{$Fc{$cM}}))*1.5 > $maxGenes){
		$Pv5 = $Pv5*0.8; $Sv5 = $Sv5*0.8;
		foreach (@cGenes){
			my $eP = exists($Pc1{$_}); my $eS = exists($Sc1{$_});
			if ( ($eP && $eS ) || ( $eP && $Pc1{$_} <= $Pv5 || $eS && $Sc1{$_} <= $Sv5)	){ 
				$Fc{$cM}{$_} = 1;#condition hard: either good corr, and present in both sets
			}
		} 
	}
	#just add the MB2 core genes.. in case..
	my $MB2lL = 0;my $MB2eL = 0;
	foreach (keys %MB1){
		if (exists($Fc{$cM}{$_})){$MB2eL++;next;}
		$Fc{$cM}{$_} = 1;$MB2lL++;} 
	$MB2l += $MB2lL;
	print "$cM: Reduced from ".scalar @cGenes." to ".scalar(keys(%{$Fc{$cM}})).", added $MB2lL/". ($MB2lL+$MB2eL) . " MB2 genes\n";
	#die;
}

my $nMGS = scalar(keys(%Fc));
my $nGenes = 0;

open O,">$outFile" or die $!;
	foreach my $mgs (keys %Fc){
		my %locH = %{$Fc{$mgs}};
		foreach my $ge (keys %locH){
			print O "$mgs\t$ge\n";
		}
		$nGenes += scalar (keys %locH);
	}
close O;

print "All Done\nWrote $nMGS MGS, with total of $nGenes genes, added $MB2l MB2 genes\n";
#system "gzip $pearI";system "gzip $speaI";







sub readCan{
	my ($inF) = @_;
	my %ret;
	open my $I1, "<$inF" or die "Can't open $inF\n";

	while (my $lin1 = <$I1>){
		chomp $lin1; my @sp1 = split /\t/,$lin1;
		$ret{$sp1[0]}{$sp1[1]} = $sp1[2];
	}
	close $I1;
	return \%ret;
}
