#!/usr/bin/env perl
#perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFNAs.fna /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFAAs.faa /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//categories4ete.txt /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2/ 12 0 0.8
use warnings;
use strict;
use threads ('yield',
                 'stack_size' => 64*4096,
                 'exit' => 'threads_only',
                 'stringify');
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( readFasta readGFF);
sub synPosOnly;



die "STOP: did this in R in the end, this code is not functional at the moment!!\n";




my @refFAs = ("/g/bork5/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa","/g/bork5/hildebra/results/TEC2/v5/T3/T3.mini2.3smpl.fna",
"/g/bork5/hildebra/results/TEC2/v5/T4/T4.mini2.3smpl.fna",
#"/g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T5/R_filt/contigs/MM3.ctgs.fna",
"/g/bork5/hildebra/results/TEC2/v5/T6/TEC6.ctgs.rn.fna");
my @refGFFs = ("/g/bork5/hildebra/results/TEC2/v5/Genomes/T2/6666666.214148.gff","","","","");
my @name = ("T2d","T3d","T4d","T6d");
my $inDir = "/g/scb/bork/hildebra/SNP/GNMass3/GlbMap/$name[$idx]/";
my $mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v5_T2subset.txt";
my $MAP = "30"; my $BQUAL = "25"; my $LOC = "ali";

my $refVCF = "$odir2/FB_reGT_${LOC}_$MAP"."_$BQUAL.vcf"
my $inGenome = $refFAs[0];
my $refGFF = $refGFFs[0];

#read reference fasta
my $hr = readFasta($inGenome); my $gFNA = %{$inGenome};

#read gff
my %gff;
open I,"<$refGFF";
while (my $line = <I>){
	chomp $line; my @spl = split /\t/,$line;
	$gff{$spl[0]}{$spl[3]}{sta} = $spl[3];
	$gff{$spl[0]}{$spl[3]}{sto} = $spl[4];
	$gff{$spl[0]}{$spl[3]}{str} = $spl[6];
	$gff{$spl[0]}{$spl[3]}{typ} = $spl[2];
}
close I;

#read vcf line by line
open I,"<$refVCF";
close I;




sub synPosOnly($ $ $ $){#not finished, I used the AA version instead
	my ($inMSA,$inAAMSA,$outMSA, $ffold) = @_;
	#print "Syn NT";
	my %convertor = (
    'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',    # Serine
    'TTC' => 'F', 'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L', 'TTG' => 'L',    # Leucine
    'TAC' => 'Y',  'TAT' => 'Y',    # Tyrosine
    'TAA' => '*', 'TAG' => '*', 'TGA' => '*',    # Stop
    'TGC' => 'C', 'TGT' => 'C',    # Cysteine   
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',    # Leucine
    'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',    # Proline
    'CAC' => 'H', 'CAT' => 'H',    # Histidine
    'CAA' => 'Q', 'CAG' => 'Q',    # Glutamine
    'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',    # Arginine
    'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',    # Threonine
    'AAC' => 'N','AAT' => 'N',    # Asparagine
    'AAA' => 'K', 'AAG' => 'K',    # Lysine
    'AGC' => 'S', 'AGT' => 'S',    # Serine
    'AGA' => 'R','AGG' => 'R',    # Arginine
    'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',    # Valine
    'GCA' => 'A','GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',    # Alanine
    'GAC' => 'D', 'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E', 'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G','GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',    # Glycine
    );
	my %ffd;
	if ($ffold){ #calc 4fold deg codons in advance to real data
		foreach my $k (keys %convertor){
			my $subk = $k; my $iniAA = $convertor{$subk} ;
			my $cnt=0;
			foreach my $sNT ( ("A","T","G","C") ){
				
				substr ($subk,2,1) = $sNT;
				#print $subk ." " ;
				$cnt++ if ($convertor{$subk} eq $iniAA);
				
			}
#			if( $cnt ==4){ $ffd{$k} = 4;
#			} else {$ffd{$k} = 1;}
			if( $cnt ==4){ $ffd{$iniAA} = 4;
			} else {$ffd{$iniAA} = 1;}
		}
	}

	#assumes correct 3 frame for all sequences in inMSA
	my $hr = readFasta($inMSA); my %FNA = %{$hr};
	my %FAA;
	if (1  || !$ffold){
		$hr = readFasta($inAAMSA); %FAA = %{$hr};
	}
	#print "$inMSA\n$inAAMSA\n$outMSA\n";
	my @aSeq = keys %FAA; my %outFNA;
	for (my $j=0;$j<@aSeq;$j++){$outFNA{$aSeq[$j]}="";}
	my $len = length ($FAA{$aSeq[0]});
	my $nsyn=0;my $syn=0;
	for (my $i=0; $i< $len; $i+=1){
		my $iniAA = substr $FAA{$aSeq[0]},$i,1; my $isSame = 1;
		next unless (!$ffold || $ffd{$iniAA} == 4);
	#print $i." $iniAA ";
		for (my $j=1;$j<@aSeq;$j++){
			my $newAA = substr $FAA{$aSeq[$j]},$i,1;
			if ($iniAA ne $newAA && $newAA ne "-"){
				$isSame =0; last;
			}
		}
		if ($isSame){#add nts to file
			for (my $j=0;$j<@aSeq;$j++){
				if ($ffold){
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},($i*3)+2,1;
				} else {
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},$i*3,3;
				}
				#print substr $FNA{$aSeq[$j]},$i*3,3 . " ";
			}
			$syn++;
		} else {$nsyn++;}
	}
	open O ,">$outMSA" or die "Can't open outMSA $outMSA\n";
	for (my $j=0;$j<@aSeq;$j++){
		print O ">$aSeq[$j]\n$outFNA{$aSeq[$j]}\n";
	}
	close O;
	$aSeq[0] =~ m/^.*_(.*)$/;
	#die "$outMSA\n";
	print "$1 ($syn / $nsyn) ".@aSeq." seqs \n";
	#print " only\n";
	#print "\n";
}

