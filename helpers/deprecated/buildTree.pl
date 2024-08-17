#!/usr/bin/env perl
#perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl -fna /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFNAs.fna -aa /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//allFAAs.faa -cats /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2//renameTEC2//categories4ete.txt -outD /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5//T2/tesssst/ -cores 12 -useEte 0 -NTfilt 0.8 -NonSynTree 1 -SynTree 1
#ARGS: ./buildTree.pl [FNA] [FAA] [categoryFile] [outDir] [CPUs] [useEte? [1=ETE,0=this script]] [filter]

use warnings;
use strict;
use threads ('yield',
                 'stack_size' => 64*4096,
                 'exit' => 'threads_only',
                 'stringify');
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( systemW readFasta);
use Getopt::Long qw( GetOptions );


sub convertMultAli2NT;
sub mergeMSAs;
sub synPosOnly;
sub runRaxML;

my $doPhym= 0;

my $pal2nal = getProgPaths("pal2nal");#"perl /g/bork3/home/hildebra/bin/pal2nal.v14/pal2nal.pl";
#die $pal2nal;
my $clustaloBin = getProgPaths("clustalo");#= "/g/bork3/home/hildebra/bin/clustalo/clustalo-1.2.0-Ubuntu-x86_64";
my $fastq2phylip = getProgPaths("fastq2phylip_scr");
my $phymlBin = getProgPaths("phyml");
my $raxmlBin = getProgPaths("raxml");
my $msapBin = getProgPaths("msaprobs");
my $trimalBin = getProgPaths("trimal");
my $pigzBin  = getProgPaths("pigz");

#die "TODO $trimalBin\n";
#trimal -in /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/COG0185.faa -out /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/tst.fna -backtrans /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/inMSA0.fna -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon
#some runtim options...
#my $ncore = 20;#RAXML cores
 my $ntFrac =2; 
 my $clustalUse = 1; #do MSA with clustal (1) or msaprobs (0) 
 if ($clustalUse == 0){print "Warning:  MSAprobs with trimal gives warnings (ignore them)\n";}

my $ntCnt =0; my $bootStrap=0;
my ($fnFna, $aaFna,$cogCats,$outD,$ncore,$Ete, $filt,$smplDef,$smplSep,$calcSyn,$calcNonSyn) = ("","","","",12,0,0.8,1,"_",1,0);
my $outgroup="";

die "no input args!\n" if (@ARGV == 0 );


GetOptions(
	"fna=s" => \$fnFna,
	"aa=s"      => \$aaFna,
	"cats=s"      => \$cogCats,
	"outD=s"      => \$outD,
	"cores=i" => \$ncore,
	"useEte=i"      => \$Ete,
	"NTfilt=f"      => \$filt,
	"smplDef=i"	=> \$smplDef, #is the genome somehow quantified with a delimiter (_) ?
	"smplSep=s" => \$smplSep, #set the delimiter
	"outgroup=s"	=> \$outgroup,
	"NonSynTree=i"	=> \$calcNonSyn,
	"SynTree=i"	=> \$calcSyn,
	"bootstrap=i" => \$bootStrap,
) or die("Error in command line arguments\n");

if ($filt <1){$ntFrac=$filt; print "Using filter with $ntFrac fraction of nts\n";}
if ($outgroup ne ""){print "Using outgroup $outgroup\n";}
if ($bootStrap>0){print "Using bootstrapping in tree building\n";}
else {$ntCnt = $filt;}
my $tmpD = $outD;
system "mkdir -p $outD" unless (-d $outD);
my $cmd =""; my %usedGeneNms;

#------------------------------------------
#sorting by COG, MSA & syn position extraction
if (!$Ete){
	my $treeD = "$outD/TMCtree/";
	system "rm -fr $treeD" if (-d $treeD);
	system "mkdir -p  $outD/MSA/" unless(-d "$outD/MSA/");
	my $multAli = "$outD/MSA/MSAli.fna";
	my $multAliSyn = $multAli.".syn.fna";
	my $multAliNonSyn = $multAli.".nonsyn.fna";

	my @MSAs; my @MSAsSyn; my @MSAsNonSyn;#full MSAs and MSAs with syn / nonsyn pos only
	my @MSrm; 
	my %FAA ; my %FNA ;
	#my @xx = keys %FAA; die "$xx[0] $xx[1]\n$FAA{HM29_COG0185}\n";
	if (1 && $cogCats ne ""){
		my $hr = readFasta($aaFna,1); my %FAA = %{$hr};
		#die $FAA{"642204217"}."\n";
		$hr = readFasta($fnFna,1); my %FNA = %{$hr};
		print "ReadFasta\n";
		open I,"<$cogCats" or die "Can't open cogcats $cogCats\n";
		my $cnt = 0; 
		my %samples; my $ogrpCnt=0;#my %genCats; 
		while (<I>){
			chomp; my @spl = split /\t/;
			@spl = grep !/NA/, @spl;#remove NAs
			next if (@spl ==0);
			$spl[0] =~ m/^(.*)$smplSep(.*)$/;
			my @spl2 = ($1,$2);#split /$smplSep/,$spl[0] ;
			my $ogrGenes = "";
			if ($outgroup ne ""){
				foreach my $seq (@spl){
					#my @spl3 = split /$smplSep/,$seq ; 
					$seq =~ m/^(.*)$smplSep(.*)$/;#my @spl3 =($1,$2);
					if ($1 eq $outgroup){$ogrGenes = $seq; $ogrpCnt ++ ; last;}
				}
			}
			
			#die "$ogrGenes $spl2[0]   $spl2[1] \n";
			
			die "Double gene name in tree build pre-concat: $spl2[1]\n" if (exists($usedGeneNms{$spl2[1]}));
			$usedGeneNms{$spl2[1]} = 1;
			
			my $tmpInMSA = "$tmpD/inMSA$cnt.faa";
			my $tmpInMSAnt = "$tmpD/inMSA$cnt.fna";
			my $tmpOutMSA2 = "$outD/MSA/$spl2[1].$cnt.faa";
			my $tmpOutMSA = "$outD/MSA/$spl2[1].$cnt.fna";
			my $tmpOutMSAsyn = "$tmpD/MSA/$spl2[1].$cnt.syn.fna";
			my $tmpOutMSAnonsyn = "$tmpD/MSA/$spl2[1].$cnt.nonsyn.fna";
			open O,">$tmpInMSA" or die "Can;t open tmp faa file for MSA: $tmpInMSA\n";
			open O2,">$tmpInMSAnt" or die "Can;t open tmp fna file for MSA: $tmpInMSAnt\n";
			foreach my $seq (@spl){
				#print "$seq\n";
				#my @spl2 = split /$smplSep/,$seq; 
				my $seq2 = $seq;
				if (!$smplDef){#create artificial head tag
					#TODO.. don't need it now for tec2, since no good NCBI taxid currently...
				}
				$seq2 =~ m/^(.*)$smplSep(.*)$/;#my @spl2 =($1,$2);
				$samples{$1} = 1; #$genCats{$spl2[1]} = 1; 
				die "can't find AA seq $seq\n" unless (exists ($FAA{$seq}));
				die "can't find fna seq $seq\n" unless (exists ($FNA{$seq}));
				$FAA{$seq} =~ s/\*//g if (!$clustalUse);
				print O ">$seq2\n$FAA{$seq}\n";
				print O2 ">$seq2\n$FNA{$seq}\n";
			}
			close O;close O2;
			if ($clustalUse){
				$cmd = $clustaloBin." -i $tmpInMSA -o $tmpOutMSA2 --outfmt=fasta --threads=$ncore --force\n";
			} else {
				$cmd = "$msapBin -num_threads $ncore $tmpInMSA > $tmpOutMSA2";
			}
			#die $cmd;
			system $cmd unless (-s $tmpOutMSA2);
			#die;
			convertMultAli2NT($tmpOutMSA2,$tmpInMSAnt,$tmpOutMSA);
			
			#die("XX\n") if ($spl2[1] eq "COG0081");
			($tmpOutMSAsyn,$tmpOutMSAnonsyn) = synPosOnly($tmpOutMSA,$tmpOutMSA2,$tmpOutMSAsyn,$tmpOutMSAnonsyn,0,$ogrGenes,$calcSyn,$calcNonSyn);
			system "rm -f $tmpInMSA $tmpInMSAnt";# $tmpOutMSA2";
			push (@MSAs,$tmpOutMSA);
			push (@MSAsSyn,$tmpOutMSAsyn) if ($tmpOutMSAsyn ne "");
			push (@MSAsNonSyn,$tmpOutMSAnonsyn) if ($tmpOutMSAnonsyn ne "");
			
			push (@MSrm,$tmpOutMSA2,$tmpOutMSA);
			$cnt ++;
			print "$cnt ";

		}
		close I;
		if ($outgroup ne ""){print "Found $ogrpCnt of $cnt outgroup sequences\n";}
		#die;
		#merge cogcats - can go to tree from here
		mergeMSAs(\@MSAs,\%samples,$multAli,0);
		mergeMSAs(\@MSAsSyn,\%samples,$multAliSyn,1) if ($calcSyn);
		mergeMSAs(\@MSAsNonSyn,\%samples,$multAliNonSyn,1) if ($calcNonSyn);
		
		#die();
		system "rm -f ".join(" ",@MSrm)." ".join(" ",@MSAsSyn); 
	} elsif (1) {#no marker way, single gene
		my $tmpInMSA = $aaFna;
		my $tmpInMSAnt = $fnFna;
		my $tmpOutMSA2 = "$tmpD/outMSA.faa";
		my $tmpOutMSA = $multAli;#"$tmpD/outMSA.fna";
		my $tmpOutMSAsyn = $multAliSyn;#"$tmpD/outMSA.syn.fna";
		my $tmpOutMSAnonsyn = $multAliNonSyn;
		my $numFas = `grep -c '^>' $tmpInMSA`;
		chomp $numFas;
		if ($numFas <= 1){print "Not enough Sequences\n"; exit(0);}
		if ($clustalUse){
			$cmd = $clustaloBin." -i $tmpInMSA -o $tmpOutMSA2 --outfmt=fasta --threads=$ncore --force\n";
		} else {
			$cmd = "sed -i 's/\\*//g' $tmpInMSA\n";
			$cmd .= "$msapBin -num_threads $ncore $tmpInMSA > $tmpOutMSA2\n";
		}
		#die $cmd;
		system $cmd; print "finished MSA\n";
		convertMultAli2NT($tmpOutMSA2,$tmpInMSAnt,$tmpOutMSA);
		synPosOnly($tmpOutMSA,$tmpOutMSA2,$tmpOutMSAsyn,$tmpOutMSAnonsyn,0,"",$calcSyn,$calcNonSyn);
		#system "rm $tmpInMSA $tmpInMSAnt $tmpOutMSA2";
		system "rm $tmpOutMSA2";
		push (@MSAs,$tmpOutMSA);
		#$multAli = $tmpOutMSA; $multAliSyn = $tmpOutMSAsyn;
	}

	
	#-------------------------------------------
	#Tree building part with RaxML
	#die $multAli."\n";
	#convert fasta again
	
	my @thrs;
	my $tcmd = "$fastq2phylip -c 50 $multAli > $multAli.ph\n";
	$tcmd .= "$fastq2phylip -c 50 $multAliSyn >$multAliSyn.ph\n" if ($calcSyn);
	$tcmd .= "$fastq2phylip -c 50 $multAliNonSyn >$multAliNonSyn.ph\n"if ($calcNonSyn);
	if (system $tcmd) {die "fasta2phylim failed:\n$tcmd\n";}
	my $BStag = ""; if ($bootStrap>0){$BStag="_BS$bootStrap";}
	runRaxML("$multAli.ph",$bootStrap,$outgroup,"$treeD/RXML_allsites$BStag.nwk");
	runRaxML("$multAliSyn.ph",$bootStrap,$outgroup,"$treeD/RXML_syn$BStag.nwk") if ($calcSyn);
	runRaxML("$multAliNonSyn.ph",$bootStrap,$outgroup,"$treeD/RXML_nonsyn$BStag.nwk") if ($calcNonSyn);

	#phyml
	if ($doPhym){
		my $nwkFile = "$treeD/tree_phyml_all.nwk";
		my $nwkFile2 = "$treeD/tree_phyml_syn.nwk";
		$tcmd = "$phymlBin --quiet -m GTR --no_memory_check -d nt -f m -v e -o tlr --nclasses 4 -b 2 -a e -i $multAli.ph > $nwkFile\n";
		push(@thrs, threads->create(sub{system $tcmd;}));
		$tcmd = "$phymlBin --quiet -m GTR --no_memory_check -d nt -f m -v e -o tlr --nclasses 4 -b 2 -a e -i $multAliSyn.ph > $nwkFile2\n";
		push(@thrs, threads->create(sub{system $tcmd;}));
	}
	for (my $t=0;$t<@thrs;$t++){
		my $state = $thrs[$t]->join();
		if ($state){die "Thread $t exited with state $state\nSomething went wrong with RaxML\n";}
	}
	system "rm  $multAli.ph $multAliSyn.ph ";
	
	#die "$distTree_scr -d -a --dist-output $raxD/distance.syn.txt $raxD/RXML_sym.nwk\n";

	
	###################### ETE ######################3
} else {
	$cmd = "ete3 build -n $fnFna -a $aaFna -w clustalo_default-none-none-none  -m sptree_raxml_all --cpu $ncore -o $outD/tree --clearall --nt-switch 0.0 --noimg  --tools-dir /g/bork3/home/hildebra/bin/ete/ext_apps-latest"; #--no-seq-checks
	$cmd .= " --cogs $cogCats" unless ($cogCats eq "");
	print "Running tree analysis ..";
	system $cmd . "> $outD/tree/ETE.log";
	print " Done.\n$outD/tree\n";
}

print "All done: $outD \n\n";
exit(0);













##########################################################################################
##########################################################################################

sub runRaxML{
	my ($mAli,$bootStrap,$outgroup,$outTree) = @_;
	my ($raTmpF,$raTmpF2,$raTmpF3) = ("RXMtmp","R2XM","RX3M"); #my $raxFile = "RXMall";
#	my $raxFile2 = "RXMsyn";
	my $raxD = $outTree; $raxD =~ s/[^\/]+$/tmp\//;
	my $raxLogD= $outTree; $raxLogD =~ s/[^\/]+$/Logs\//;
	system "mkdir -p  $raxLogD" unless (-d $raxLogD);
	system "rm -rf $raxD;mkdir -p $raxD";
	my $outGrpRXML = "";
	if ($outgroup ne ""){$outGrpRXML = "-o $outgroup";}
	my $raxDef = " --silent -m GTRGAMMA -p 312413 ";
	#raxml - on all sites
	my $tcmd="";
	my $bootCmd = "";
	if ($bootStrap > 0){$bootCmd.= "-N $bootStrap -b ". int(rand(10000));}
	$tcmd =  "$raxmlBin -T$ncore -f d -s $mAli $raxDef -n $raTmpF -w $raxD $outGrpRXML > $raxLogD/ini.log\n";
	#die "$tcmd\n";
	if (system $tcmd){#failed.. prob optimization problem, use other tree instead
		$raxDef = " --silent -m GTRGAMMAI -p 312413 ";	system "rm $raxD/*$raTmpF";
		systemW "$raxmlBin -T$ncore -f d -p 31416 -s $mAli $raxDef -n $raTmpF -w $raxD $outGrpRXML > $raxLogD/optimized.log\n";
	}
	#decide which support vals to calc
	if ($bootStrap==0){
		systemW "$raxmlBin -T$ncore -f J -s $mAli $raxDef -n $raTmpF3 -w $raxD -t $raxD/RAxML_bestTree.$raTmpF $outGrpRXML > $raxLogD/optimized.log\n";
		system "mv $raxD/RAxML_fastTreeSH_Support.$raTmpF3 $outTree";
	} else {
		systemW  "$raxmlBin -T$ncore -f d  -s $mAli $raxDef -n $raTmpF2 -w $raxD $bootCmd $outGrpRXML > $raxLogD/ini.log\n";
		systemW  "$raxmlBin -T$ncore -f b  -s $mAli $raxDef -n $raTmpF3 -w $raxD -t $raxD/RAxML_bestTree.$raTmpF -z $raxD/RAxML_bootstrap.$raTmpF2 $outGrpRXML > $raxLogD/ini.log\n";
		system "mv $raxD/RAxML_bipartitions.$raTmpF3 $outTree";
	
	}
	#and create distance matrix
	systemW  "$raxmlBin -T$ncore -f x -s $mAli $raxDef -n all -w $raxD -t $raxD/RAxML_bestTree.$raTmpF $outGrpRXML > $raxLogD/.log\n";
	my $outDist = $outTree; $outDist =~ s/\.[^\.]+$/\.dist/;
	#die "$outDist\n";
	system "mv $raxD/RAxML_distances.all $outDist";

	
	#clean up
	system "rm $raxD/*_info* $raxD/*_log* $raxD/*_parsimony* $raxD/*_result*";
}




sub mergeMSAs($ $ $ $){
	my ($MSAsAr,$samplesHr,$multAliF,$del) = @_;
	my @MSAs = @{$MSAsAr}; my %samples = %{$samplesHr};
	my %bigMSAFAA;foreach my $sm (keys %samples){$bigMSAFAA{$sm} ="";}
	foreach my $MSAf (@MSAs){
		#print $MSAf."\n"; 
		my $hit =0; my $miss =0;
		my $hr = readFasta($MSAf,1); my %MFAA = %{$hr};
		system "rm $MSAf" if ($del);
		my @Mkeys = keys %MFAA;
		#die "$Mkeys[0]\n";
		
#		my @spl2 = split /$smplSep/,$Mkeys[0];
		$Mkeys[0] =~ m/^(.*)$smplSep(.*)$/;my @spl2 =($1,$2);
		my $gcat = $spl2[1];
		my $len = length( $MFAA{$Mkeys[0]} );
		#die $len;
		foreach my $sm (keys %samples){
			my $curK = $sm.$smplSep.$gcat; #print $curK. " ";
			#print "$MFAA{$curK}\n";
			if (exists $MFAA{$curK}){
				$bigMSAFAA{$sm} .= $MFAA{$curK}; $hit++;
			} else {
				$bigMSAFAA{$sm} .= "-"x$len; $miss++;#print "nooooooo ";
			}
		}
		
		#die "$hit - $miss\n";
	}
	#filter part - count "-" in each seq
	
	my @ksMSAFAA = keys %bigMSAFAA;
	my $iniSeqNum = @ksMSAFAA; my $remSeqNum = 0;
	my %ntCnts; my $maxNtCnt=1;
	foreach my $kk (@ksMSAFAA){
		my $strCpy = $bigMSAFAA{$kk};
		my $num1 = $strCpy =~ tr/[\-N]//;
		#$num1++ while ($bigMSAFAA{$kk} !~ m/-/g);
		#now get actual num of bps
		$num1 = length($strCpy) - $num1;
		if ($num1 > $maxNtCnt){$maxNtCnt = $num1;}
		$ntCnts{$kk} = $num1;
	}
	foreach my $kk (@ksMSAFAA){
		my $num1 = $ntCnts{$kk};
		if ( ( ($num1 / $maxNtCnt ) < $ntFrac) || ($num1 < $ntCnt ) ){
			delete $bigMSAFAA{$kk}; $remSeqNum++; 
			print "$kk $num1  ".int($num1*1000 / ($maxNtCnt) )/1000 ." $ntFrac \n";
		}
		#print "$num1  $kk \n";#$bigMSAFAA{$kk}\n\n"; last;
	}
	open O,">$multAliF" or die "Can't open MSA outfile $multAliF\n";
	foreach my $kk (keys %bigMSAFAA){
		print O ">$kk\n$bigMSAFAA{$kk}\n";
	}
	close O;
	#die "$multAliF\n";
	print "Removed $remSeqNum of $iniSeqNum sequences\n";
}
sub convertMultAli2NT($ $ $){
	my ($inMSA,$NTs,$outMSA) = @_;
	my $tmpMSA=0;
	if ($inMSA eq $outMSA){$outMSA .= ".tmp"; $tmpMSA=1;}
	my $cmd = "";
	#"$pal2nal $inMSA $NTs -output fasta -nostderr -codontable 11 > $outMSA\n";
	
	#$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon\n";
	$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -ignorestopcodon  -gt 0.1 -cons 60\n";
	#die "$inMSA,$NTs,$outMSA\n";
	#my $hr1= readFasta($inMSA);
	#my %MSA = %{$hr1};
	#$hr1= readFasta($NTs);
	#my %NTs = %{$hr1};
	if ($tmpMSA){$cmd .= "rm $inMSA;mv $outMSA $inMSA;\n";}
	#print $cmd;
	die "Can't execute $cmd\n" if (system $cmd);
}

sub synPosOnlyAA($ $){#only leaves "constant" AA positions in MSA file.. 
#stupid, don't know if pal2nal can handle this.. prob not
	my ($inMSA,$outMSA) = @_;
	print "Syn";
	my $hr = readFasta($inMSA,1); my %FNA = %{$hr};
	my @aSeq = keys %FNA;
	my $len = length ($FNA{$aSeq[0]});
	for (my $i=0; $i< $len; $i+=3){
		my $cod = substr $FNA{$aSeq[0]},$i,3;
		my $iniAA = "A";
		for (my $j=1;$j<@aSeq;$j++){
		}
	}
	print " only\n";

}

sub synPosOnly{#now finished, version is cleaner
	my ($inMSA,$inAAMSA,$outMSA, $outMSAns, $ffold, $outgroup, $doSyn, $doNSyn) = @_;
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
	my $hr = readFasta($inMSA,1); my %FNA = %{$hr};
	#my %FAA;
	#if (0  || !$ffold){
	#	$hr = readFasta($inAAMSA); %FAA = %{$hr};
	#}
	#print "$inMSA\n$inAAMSA\n$outMSA\n";
	my @aSeq = keys %FNA; 
	my %outFNA;#syn
	my %outFNAns;#non syn
	for (my $j=0;$j<@aSeq;$j++){$outFNA{$aSeq[$j]}="";}
	my $len = length ($FNA{$aSeq[0]});
	my $nsyn=0;my $syn=0;
	for (my $i=0; $i< $len; $i+=3){ #goes over every position
		my $j =0;
		my $iniAA = "-";
		while (1){ #check for first informative position
			my $iniCodon = substr $FNA{$aSeq[$j]},$i,3;
			if ($iniCodon =~ m/---/ || $iniCodon =~ m/N/){$j++; last if ($j >= @aSeq); next;}
			die "error: $iniCodon\n" if ($iniCodon =~ m/-/); #should not happen
			$iniAA = $convertor{$iniCodon};#substr $FAA{$aSeq[0]},$i,1; 
			last;
		}
		#die "$iniAA\n";
		my $isSame = 1;
		next unless (!$ffold || $ffd{$iniAA} == 4);
	#print $i." $iniAA ";
		for (;$j<@aSeq;$j++){
			if ($aSeq[$j] eq $outgroup){next;}#print "HIT   $aSeq[$j] eq $outgroup";
			my $newCodon = substr $FNA{$aSeq[$j]},$i,3;
			my $newAA = "-";
			if ($newCodon !~ m/-/ && $newCodon =~ m/[ACTG]{3}/i){
				die "Unkown AA $newCodon\n" unless (exists $convertor{$newCodon} );
				$newAA = $convertor{$newCodon} ; # substr $FAA{$aSeq[$j]},$i,1;
			} elsif ($newCodon =~ m/---/ || $newCodon =~ m/[NWYRSKMDVHB]/){
			} else {
				die "newCodon wrong $newCodon\n" ;
			}
			if ($iniAA ne $newAA && $newAA ne "-"){
				$isSame =0; last;
			}
		}
		if ($isSame){#add nts to file
			for (my $j=0;$j<@aSeq;$j++){
				if ($ffold){
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},($i)+2,1;
				} else {
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},$i,3;
				}
				#print substr $FNA{$aSeq[$j]},$i*3,3 . " ";
			}
			$syn++;
		} else {
			for (my $j=0;$j<@aSeq;$j++){
				$outFNAns{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},$i,3;
			}
			$nsyn++;
		}
	}
	#die $inMSA."\n";
	if ($doSyn){
		if ($syn ==0){
			$outMSA = "";
		} else {
			open O ,">$outMSA" or die "Can't open outMSA $outMSA\n";
			for (my $j=0;$j<@aSeq;$j++){
				print O ">$aSeq[$j]\n$outFNA{$aSeq[$j]}\n";
			}
			close O;
		}
	}
	if ($doNSyn){
		if ($nsyn ==0){
			$outMSAns = "";
		} else {
			open O ,">$outMSAns" or die "Can't open outMSA $outMSAns\n";
			for (my $j=0;$j<@aSeq;$j++){
				print O ">$aSeq[$j]\n$outFNAns{$aSeq[$j]}\n";
			}
			close O;
		}
	}
	$aSeq[0] =~ m/^.*_(.*)$/;
	#die "$outMSA\n";
	print "$1 ($syn / $nsyn) ".@aSeq." seqs \n";
	#print " only\n";
	#print "\n";
	return ($outMSA,$outMSAns);
}






