#!/usr/bin/env perl
#combines ribosomal predictions from several TamoC runs
use warnings;
use strict;
use File::Basename;
use List::MoreUtils 'first_index'; 

sub readFasta;
sub FQ_nameFmt;
use Mods::GenoMetaAss qw(readMap);


#my $baseOut = "/g/scb/bork/hildebra/SNP/GNMass/combinations/";my $mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM.txt";
#Soil samples
my $mapF = "/g/scb/bork/hildebra/data2/Soil_finland/soil_map2.txt";
my $baseOut = "/g/scb/bork/hildebra/FinSoilAll/";


my $protOutD = $baseOut."FMD";
my $riboOutD = $baseOut."Ribo";
system ("rm -r $protOutD\nmkdir -p $riboOutD");
my ($hr,$hr2) = readMap($mapF);
my %map = %{$hr};
my @samples = @{$map{opt}{smpl_order}};

print "Combining markers in ".@samples" samples..\n";

my @baseRiF = ("reads_SSU","reads_LSU","reads_ITS");
my $rlogs="Sample\tSSU_RiboReads\t18S\t16Sbac\t16Sarc\tLSU_RiboReads\t28S\t23S_bac\t23S_arc\tITS\n";
#cleanup
foreach my $bRF (@baseRiF){open O,">$riboOutD/$bRF.c.1.fq"or die "Can't open Ribo\n"; close O;open O,">$riboOutD/$bRF.c.2.fq"; close O;}

foreach my $smpl(@samples){
	my $dir2rd = $baseOut.$map{$smpl}{dir};
	my $SmplName = $map{$smpl}{SmplID};
	print $SmplName."\n";
	#FMG genes
	if (0){
		my $inFMGd = "$dir2rd/assemblies/metag/ContigStats/FMG/";
		next unless (-e "$dir2rd/assemblies/metag/scaffolds.fasta.filt" && -e "$dir2rd/assemblies/metag/genePred/genes.fna.shrtHD.fna" && -d $inFMGd);
		print "==== ".$dir2rd." ====\n";
		opendir(DIR, $inFMGd) or die "$inFMGd:".$!;	my @protsFiles = sort ( grep { /^COG\d+\.faa$/ && -f "$inFMGd/$_" } readdir(DIR) );	rewinddir DIR;
		my @geneFiles = sort ( grep { /^COG\d+\.fna$/  && -f "$inFMGd/$_"} readdir(DIR) );	close(DIR);
		foreach my $gf (@geneFiles){
			next unless (-e $inFMGd.$gf);
			#die $inFMGd.$gf."\n";
			$hr = readFasta($inFMGd.$gf);
			my %genes = %{$hr};
			open O,">>$protOutD/$gf";
			foreach (keys %genes){print O ">$_\n$genes{$_}\n";}
			close O;
			$gf =~ s/\.fna/\.faa/;
			#die $inFMGd.$gf."\n";
			my $hr = readFasta($inFMGd.$gf);
			my %prots = %{$hr};
			open O,">>$protOutD/$gf";
			foreach (keys %prots){print O ">$_\n$prots{$_}\n";}
			close O;
	#		die "$protOutD/$gf"."\n";
		}
	}
	if (1){
		#ribosomal genes
		my $inFMGd = "$dir2rd/ribos/";
		foreach my $bRF (@baseRiF){ #"reads_SSU","reads_LSU","reads_ITS"
			#change FQ headName
			FQ_nameFmt("$dir2rd/ribos/$bRF.r1.fq","$riboOutD/$bRF.c.1.fq",$SmplName."__RIB_",1);
			FQ_nameFmt("$dir2rd/ribos/$bRF.r2.fq","$riboOutD/$bRF.c.2.fq",$SmplName."__RIB_",2);
			#system "cat $dir2rd/ribos/$bRF.r1.t.fq >> $riboOutD/$bRF.c.1.fq";
			#system "cat $dir2rd/ribos/$bRF.r2.t.fq >> $riboOutD/$bRF.c.2.fq";
			my $log = `cat $dir2rd/ribos/$bRF.log`;
			$log =~ m/Total reads passing E-value threshold = (\d+) \(/;
			my $totRds=$1;
			$log =~ m/silva-euk-.*\.fasta		(\d\.\d+)%\n.*silva-bac-.*\.fasta		(\d\.\d+)%\n.*silva-arc-.*\.fasta		(\d\.\d+)%/;
			#euk,bac,arc
			if ($bRF eq "reads_SSU"){
				$rlogs .= "$SmplName";
			}
			$rlogs .= "\t$totRds";
			$rlogs .= "\t".$totRds*$1/($1+$2+$3)."\t".$totRds*$2/($1+$2+$3)."\t".$totRds*$3/($1+$2+$3) unless ($bRF eq "reads_ITS");
		#die "$riboOutD/$bRF.c.1.fq\n";
		}
		#check in log file for % different DBs
	}
	$rlogs.="\n";
}
print "FINISH\n $protOutD\n";
open O,">$riboOutD/riboAccum.log";
print O $rlogs;
close O;
sub FQ_nameFmt(){
	my ($inF,$oF,$NN,$read_pair)=@_;
	my $cntR=0;my $cnt=0;
	open I,"<$inF" or die "Can't open $inF\n";
	open O,">>$oF" or die "Can't open $oF\n";
	
	while (my $line = <I>){
		
		if ($cnt % 4 ==0 ){
			if ($line !~m/^@/){
				die "FQ formatter out of read order: $cnt\n$line\n";
			}
			print O "\@$NN$cntR/$read_pair\n"; $cntR++;
		} else {
			print O $line;
		}
		$cnt++;
	}
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

sub readMapX{
	my ($inF) = @_;
	my $cnt =-1; my %ret; my @order = ();
	my $dirCol = 1; my $smplCol = 0; my $rLenCol = -1; my $SeqTech = -1;
	open I,"<$inF" or die "Can't open $inF\n";
	while (<I>){
		chomp;	$cnt++;my @spl = split(/\t/);
		if ($cnt == 0){
			$dirCol = first_index { /Path/ } @spl;
			$smplCol = first_index { /#SmplID/ } @spl;
			$SeqTech = first_index { /SeqTech/ } @spl;
			$rLenCol = first_index { /ReadLength/ } @spl;
			next;} #maybe later check for col labels etc
		
		#die $spl[0]." ".$spl[1]."\n";
		$ret{$spl[0]}{dir} = $spl[$dirCol];
		$ret{$spl[0]}{SmplID} = $spl[$smplCol];
		if ($SeqTech >= 0) {$ret{$spl[0]}{SeqTech} = $spl[$SeqTech]; } else {$ret{$spl[0]}{SeqTech} = "";}
		if ($rLenCol >= 0){$ret{$spl[0]}{readLength} = $spl[$rLenCol];} else {$ret{$spl[0]}{readLength} = 0;}
		push(@order,$spl[0]);
	}
	close I;
	$ret{opt}{smpl_order} = \@order;
	return \%ret;
}
