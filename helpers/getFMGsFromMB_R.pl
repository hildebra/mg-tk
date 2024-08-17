#!/usr/bin/perl
#redundant with getMarkersMultSmpl.pl
#the R routines can give a set of genes that correspond to markers; these are extracted from faa/fna here
#if RfileD emtpy, do this for all TECs; if number, extract this number
#./getFMGsFromMB_R.pl /g/scb/bork/hildebra/SNP/GNMass2_singl2/alien-11-374-0/Binning/MaxBin/FMG/ /g/scb/bork/hildebra/SNP/GNMass2_singl2/alien-11-374-0 /g/scb/bork/hildebra/SNP/GNMass2_singl2/alien-11-374-0/Binning/MaxBin/FMG/genes/
#./getFMGsFromMB_R.pl "" /g/scb/bork/hildebra/SNP/GNMass2_singl2/alien-11-374-0 
use warnings;
use strict;
# another implementation is in getFMGfnaAllSmpls.pl
sub readFasta;

my $rfilesD = $ARGV[0];
my $inD = $ARGV[1];
#my $outD = $ARGV[2]; #/g/scb/bork/hildebra/SNP/GNMass2_singl2/alien-11-374-0/Binning/MaxBin/FMG/genes/
my $mode = "bat";
my $outD = $inD."/Binning/MetaBat/FMG/genes/";
system "mkdir -p $outD";
my $cntStart = 0;
my $metaGD = `cat $inD/assemblies/metag/assembly.txt`; chomp $metaGD;
my $FMGfile = $FMGfile."/ContigStats/FMG/FMGids.txt";
my $primBinFile = "$inD/Binning/MaxBin/..";
if ($mode eq "bat"){
	my $nms = chomp `cat $inD/Binning/MetaBat/MeBa.sto`;
	$primBinFile = "$inD/Binning/MetaBat/$nms.fasta.fna";
	$cntStart=1;
}

system "rm -rf $outD; mkdir -p $outD";

my $metaGD = `cat $inD/assemblies/metag/assembly.txt`; chomp $metaGD;
#die("$metaGD/genePred/proteins.shrtHD.faa");

if (-d $rfilesD){
	opendir(DIR, $rfilesD) || die "can't opendir $rfilesD: $!";
	my @TECs = grep { /^TEC/ && -f "$rfilesD/$_" } readdir(DIR);
	closedir DIR;
{ else { #create these files on the fly
	system "mkdir -p $rfilesD";
	#get list of TECs by Binning method
	my %TECs;
	open I,"<$primBinFile";
	while (<I>){
		chomp; my @spl = split /\t/; next if ($spl[1] <$cntStart);
		$TECs{$spl[0]} = $spl[1];
	}
	close I;
	open I,"<$FMGfile" or die "cant open $FMGfile";
	while (<I>){
		chomp; my @spl = split /\s/; my @s2 = split /;_/,$spl[0];
		if (exists($TECs{$s2[0]}){
			die $s2[0]." XX ".$spl[0]."\n";
		}
	}
	close I;
}

#read full gene set
my $hr = readFasta("$metaGD/genePred/proteins.shrtHD.faa");
my %prots = %{$hr};
$hr = readFasta("$metaGD/genePred/genes.shrtHD.fna");
my %genes = %{$hr};


my %seen;
foreach my $T (@TECs){
	open I,"<$rfilesD/$T" or die "doesnt exist: $rfilesD/$T\n";
	my $Thd = $T; $Thd =~ s/\.txt//;
	while (my $line = <I>){
		$line =~ s/\r\n?/\n/; chomp($line); 
		my @spl = split("\t",$line);
		if (exists($seen{$Thd}{$spl[0]})){print "Double COG: $Thd\t$spl[0]\n";next;}
		if (!exists($prots{$spl[0]})){ die "Can't find protein for $spl[0]\n";}
		if (!exists($genes{$spl[0]})){ die "Can't find protein for $spl[0]\n";}
		
		open O,">>$outD/".$spl[1].".faa";
		print O ">$Thd\n$prots{$spl[0]}\n";
		close O;
		open O,">>$outD/g".$spl[1].".fna";
		print O ">$Thd\n$genes{$spl[0]}\n";
		close O;
	}
	close I;
}



sub readFasta($){
  my ($fil) = @_;
  my %Hseq;
  if (-z $fil){ return \%Hseq;}
  open(FAS,"<","$fil") || die("Couldn't open FASTA file $fil.");
    
     my $temp; 
     my $line; my $hea=<FAS>; chomp ($hea);
      my $trHe = ($hea);
      my @tmp = split(" ",$trHe);
      $trHe = substr($tmp[0],1);
      # get sequence
    while($line = <FAS>)
    {
      #next if($line =~ m/^;/);
      if ($line =~ m/^>/){
        chomp($line);
        $Hseq{$trHe} = $temp;
        $trHe = ($line);
        @tmp = split(" ",$trHe);
	$trHe = substr($tmp[0],1);
	$trHe =~ s/\|//g;
	
#        print $trHe."\n";
        $temp = "";
        next;
      }
    chomp($line);
    $line =~ s/\s//g;
    $temp .= ($line);
    }
    $Hseq{$trHe} = $temp;
  close (FAS);
    return \%Hseq;
}
