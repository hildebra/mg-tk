package Mods::Binning;
use Exporter qw(import);
our @EXPORT_OK = qw(runMetaBat runSemiBin  runMetaDecoder runCheckM runCheckM2 
				
				createBin2 createBinFAA createBinCtgs
				readMGS readMGSrev deNovo16S readMGSrevRed minQualFilter 
				filterMGS_CM MB2assigns calcLCAcompl readCMquals);

use warnings;
use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw (systemW readFasta gzipopen getAssemblPath);
use Mods::TamocFunc qw (cram2bsam);




sub readCMquals{
	my ($IQ) = @_;
	
	my %rQ;
	my $CM2mode=0; $CM2mode = 1 if ($IQ =~m/\.cm2/);
	open I,"<$IQ" or die "Can't open maxbin2 quality $IQ\n";
	while (<I>){
		chomp; my @spl  = split /\t/;
		next if ($spl[0] eq "Bin Id" || $spl[0] eq "Name");
		#die "can't find Bin \"$spl[0]\"\n" unless (exists($ret{$spl[0]}));
		my $Bin = shift @spl;
		if ($CM2mode){
			$rQ{$Bin}{compl} = $spl[0];
			$rQ{$Bin}{conta} = $spl[1];
			$rQ{$Bin}{hetero} = 0;
			#$rQ{$Bin}{line} = join("\t",@spl);
		} else {
			$rQ{$Bin}{compl} = $spl[10];
			$rQ{$Bin}{conta} = $spl[11];
			$rQ{$Bin}{hetero} = $spl[12];
			#$rQ{$Bin}{line} = join("\t",@spl);
		}
	}
	close I;
	return (\%rQ);
}

sub minQualFilter($ $ $ $ $){
	my ($hr1,$hr2,$Compl, $Conta, $LCAcompl) = @_;
	#%MB = %{$hr1}; %MBQ = %{$hr2};
	#foreach my $bin (keys %MB){
	my $prevS = scalar(keys %{$hr1});
	#print "$prevS , " . scalar(keys %{$hr2}) . "size\n";
	foreach my $bin (keys %{$hr1}){
		if ($hr2->{$bin}{compl}< $Compl || $hr2->{$bin}{conta}> $Conta || scalar(@{$hr1->{$bin}}) == 0){
			delete ($hr1->{$bin});
			delete ($hr2->{$bin});
		}
	}
	my $delEntries = $prevS - scalar(keys %{$hr1});
	print "Deleted $delEntries MAGs due to not meeting min qual criteria (Compl:$Compl, Conta:$Conta, LCA:$LCAcompl)\n" if ($delEntries);
	return ($hr1, $hr2);
}

sub MB2assigns($ $){
	my ($inF,$IQ) = @_;
	my %ret;
	#print "$inF\n";
	open I,"<$inF" or die "Can't open Binner output $inF\n";
	while (<I>){
		chomp; my @spl  = split /\t/;
		next if ($spl[1] eq "0");
		push(@{$ret{$spl[1]}}, $spl[0]);
	}
	close I;
	
	my $rQHR = readCMquals($IQ);
	#print "$inF, $IQ ". scalar(keys %{$rQHR}) ."\n";

	return (\%ret,$rQHR);
}


#return how congruent LCA assignments of marker genes in MAG are
sub calcLCAcompl{
	my ($genesAR,$LCAHR)  = @_;
	my @genes = @{$genesAR}; my %LCA  = %{$LCAHR};
	#print "@genes\n";
	my %LCAcnts; my $hits =0 ;
	foreach my $gen (@genes){
		next unless (exists($LCA{$gen}));
		$hits++;
		my @locLCA = @{$LCA{$gen}};
		for (my $i=0;$i<@locLCA;$i++){
			$LCAcnts{$i}{$locLCA[$i]}++;
		}
	}
	#calc completeness per LCA level
	my @lvls = sort { $a <=> $b } keys %LCAcnts;
	my @complPerLvl = (); my @HitsPerLvl = ();
	foreach my $LV (@lvls){
		my %entr = %{$LCAcnts{$LV}};
		#my $ignCnt = $entr{"?"} if (exists($entr{"?"}));
		my $maxCnt = 0; my $sumCnt=0;
		foreach my $k (keys %entr){
			#print "$k $entr{$k} ; ";
			if ($k eq "?"){ next;}
			$sumCnt += $entr{$k};
			$maxCnt = $entr{$k} if ($maxCnt < $entr{$k});
		}
		if ($sumCnt > 0){
		$complPerLvl[$LV] = $maxCnt/$sumCnt;
		} else {$complPerLvl[$LV] =  0;}
		$HitsPerLvl[$LV] = $sumCnt;
		#print "  .. $maxCnt $sumCnt $complPerLvl[$LV] \n";
	}
	#scale error on higher tax levels for compound scoring
	my $maxL = @lvls;
	if ($maxL ==0) {return 0;}
	my $maxScore = 0;my $LCAcompl = 0;
	foreach my $LV (@lvls){
		$LCAcompl += ($maxL-$LV) * $HitsPerLvl[$LV]/$hits * $complPerLvl[$LV];
		$maxScore+=($maxL-$LV) * $HitsPerLvl[$LV]/$hits;
		#print "L$LCAcompl $maxScore :: ";
	}
	if ($maxScore == 0) {return 0;}
	#print "	$LCAcompl /= $maxScore;";
	$LCAcompl /= $maxScore;
	#print "Found $hits genes with LCA\n";
	#print $LCAcompl."\n";
	#die;
	return $LCAcompl;
}



#get list of MGS, based on checkM filtering
sub filterMGS_CM{
	my $CMfile=$_[0];
	my $complT=50; my $contT=5; my $retBetter=1;
	$complT=$_[1] if (@_ > 1); 
	$contT=$_[2] if (@_ > 2); 
	$retBetter=$_[3] if (@_ > 3); 
	my %ret; my $cnt =0;
	my $CMtag = "CM";
	my $contIdx = 12; my $complIdx = 11;
	if ($CMfile =~ m/\.cm2/){
		$CMtag = "CM2";
		$contIdx = 2; $complIdx = 1;
	}
	open I,"<$CMfile" or die "Binning.pm::filterMGS_CM: can't open $CMtag file $CMfile\n";
	while (<I>){
		$cnt++; next if ($cnt == 1);
		chomp;my @spl=split/\t/;
		#next if ($spl[0] eq "Bin Id");
		if (($spl[$complIdx] <= $complT || $spl[$contIdx] > $contT) ){
			$ret{$spl[0]} = [$spl[$complIdx],$spl[$contIdx]] if (!$retBetter);
			#die "$spl[11],$spl[12]\n" if ($spl[0] eq "MB2bin12");
		} elsif ($retBetter) {
			$ret{$spl[0]} = [$spl[$complIdx],$spl[$contIdx]];
		}
	}
	close I;
	return(\%ret);
}

 sub deNovo16S{
	#returns 1)hash with rDNAs and 2) name of best seq
	my ($fasFile,$outfile) = @_;
	die "Fasta $fasFile doesn't exist\n" unless (-e $fasFile);
	my %ret; my $refSize = 1580;
	print "Detecting de novo 16S\n";
	
	my $newGff = "$outfile.ribo.gff";
	system "rm -f $newGff" if (-e $newGff);
	my $rnaBin = getProgPaths("rnammer");
	my $cmd = "$rnaBin -S bac -m ssu -gff $newGff < $fasFile";
	system $cmd."\n";
	#FP929038        RNAmmer-1.2     rRNA    507757  508732  1202.9  -       .       16s_rRNA
	my $genoR = readFasta($fasFile,1); my %geno = %{$genoR};
	#my @k = keys(%geno); die @k.join(" ",@k)."\n";
	open II,"<",$newGff; my $gffcnt=0;
	while(<II>){next if (/^#/ || length($_) < 1);my @spl = split(/\s+/);
		#print $_;
		my $newS = substr($geno{$spl[0]},$spl[3],$spl[4]-$spl[3]);
		if ($spl[6] eq "-"){$newS = reverse_complement_IUPAC($newS);}
		#die $newS;
		$ret{$spl[0]."_rrn_".$gffcnt} = $newS;
		$gffcnt++;
	}
	close II;
	#check for new better 16S
	my $bestV=100000;
	my @head = keys(%ret);	my $best = 0; 	my $cnt =0;
	if (@head == 0){
		#die "Empty array\n$newGff\n\n$cmd\n";
		print "Empty 16S array\n$newGff\n";
		return (\%ret,"");
	} else {
		foreach my $hd (@head){
			my $tmp = abs($refSize - length($ret{$hd}));
			if ($tmp < $bestV){$bestV = $tmp; $best = $cnt;}
			$cnt++;
			#print $tmp."\n";
		}
	}
	#print "16S ".$bestV."\n";
	my %ret3 = ($head[$best],$ret{$head[$best]});
	return (\%ret,$head[$best]);
}

sub readBinSB($){
	my ($inF) = @_;
	my %can2gene;
	print "reading SemiBin file: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		push(@{$can2gene{$spl[1]}},$spl[0]);
	}
	close I;
	return (\%can2gene);

}

sub readMGS($){
	my ($inF) = @_;
	my %can2gene;
	print "reading MGS file: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		push(@{$can2gene{$spl[0]}},$spl[1]);
	}
	close I;
	return (\%can2gene);
}
sub readMGSrev($){
	my ($inF) = @_;
	my %g2c;
	my $dbl=0;
	print "reading MGS reverse: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		#next unless (defined $spl[0] && defined $spl[1]);
		$dbl++ if (exists($g2c{$spl[1]}));
		$g2c{$spl[1]} = $spl[0];
	}
	close I;
	print "$dbl double gene assignments\n";
	return (\%g2c);
}

sub readMGSrevRed($){
	my ($inF) = @_;
	my %g2c;
	my $dbl=0;
	print "reading MGS reverse: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		#next unless (defined $spl[0] && defined $spl[1]);
		$dbl++ if (exists($g2c{$spl[1]}));
		$g2c{$spl[1]}{ $spl[0]} = 1;
	}
	close I;
	print "$dbl double gene assignments\n";
	return (\%g2c);
}


sub getRepresentBins{
	my ($guide) = @_;
	print "Reading guidMGS file $guide\n";
	my ($I,$OK) = gzipopen($guide,"MGSvsGC",1);
	my @hd = (); my $lastMGS = "";
	my $repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;
	my %ret;
	while (my $line = <$I>){
		chomp $line; my @spl =  split /\t/,$line;
		if (!@hd){
			@hd = @spl;
			while ($hd[$repIdx] ne "Representative4MGS" && ($repIdx+2) < @hd){$repIdx++;}
			die "Couldn't find \"Representative4MGS\" string in @hd\n!\n" if ($repIdx >= @hd) ;
			while ($hd[$ComplIdx] ne "Completeness" && ($ComplIdx+2) < @hd){$ComplIdx++;}
			while ($hd[$ContaIdx] ne "Contamination" && ($ContaIdx+2) < @hd){$ContaIdx++;}
			while ($hd[$LCAidx] ne "LCAcompleteness" && ($LCAidx+2) < scalar(@hd)){$LCAidx++;}
			while ($hd[$N50idx] ne "N50" && ($N50idx+2) < @hd){$N50idx++;}
			#die "$repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;\n@hd \n";
			next;
		}
		
		my $curMGS = $spl[1];
		if ($curMGS ne $lastMGS){ #new cycle
			$lastMGS = $curMGS;
			$ret{$curMGS} = $spl[0] if ($spl[0] !~ m/^Cano__/); #set default to last one..
		}
		if ($spl[0] !~ m/^Cano__/ && ($spl[$repIdx] eq "*" || !defined($ret{$curMGS}) ) ){
			$ret{$curMGS} = $spl[0] ;
			#print "MGS $curMGS repr: $spl[0]\n";
		}
	}
	close $I;
	
	print "Found representative for ". scalar(keys %ret). " MGS\n";
	return \%ret;
}




sub getRepresentBinsPerFamily{ #needs some work
	my ($guide,$hrMap) = @_;
	print "Reading guidMGS file $guide\n";
	my ($I,$OK) = gzipopen($guide,"MGSvsGC",1);
	my @hd = (); my $lastMGS = "";
	my $repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;
	my $rejQuali = 0; my $famFnd=0; my $assmGrpFnd=0; my $noGrpFnd=0;
	my %ret;
	my %map = %{$hrMap};
	
	
	my %conta; my %compl; my %fam; my %score; my %MGS2bin;
	
	
	#read in family info for each sample
	my %famSmpl;
	my @smpls = @{$map{opt}{smpl_order}};
	foreach my $smpl (@smpls){
		#$famSmpl{$smpl} = $map{$smpl}{AssGroup};
		if ($map{$smpl}{FamGroup} ne ""){
			$famSmpl{$smpl} = $map{$smpl}{FamGroup};$famFnd++;
		} elsif ($map{$smpl}{AssGroup} ne ""){#fallback assembly group
			$famSmpl{$smpl} = $map{$smpl}{AssGroup}; $assmGrpFnd++;
		} else {
			$famSmpl{$smpl} = $smpl; $noGrpFnd++;
		}
	}

	while (my $line = <$I>){
		chomp $line; my @spl =  split /\t/,$line;
		if (!@hd){
			@hd = @spl;
			while ($hd[$repIdx] ne "Representative4MGS" && ($repIdx+2) < @hd){$repIdx++;}
			die "Couldn't find \"Representative4MGS\" string in @hd\n!\n" if ($repIdx >= @hd) ;
			while ($hd[$ComplIdx] ne "Completeness" && ($ComplIdx+2) < @hd){$ComplIdx++;}
			while ($hd[$ContaIdx] ne "Contamination" && ($ContaIdx+2) < @hd){$ContaIdx++;}
			while ($hd[$LCAidx] ne "LCAcompleteness" && ($LCAidx+2) < scalar(@hd)){$LCAidx++;}
			while ($hd[$N50idx] ne "N50" && ($N50idx+2) < @hd){$N50idx++;}
			#die "$repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;\n@hd \n";
			next;
		}
		
		
		
		my $curMGS = $spl[1]; my $curBin = $spl[0];
		if ($curMGS ne $lastMGS){
			#eval last MGS for family reps
			foreach my $ff (keys %fam){
				my @curF = @{$fam{$ff}};
				#select best MGS in each family
				my $bestSc = 0; my $bestBin = "";
				foreach my $lBin (@curF){
					if ($score{$lBin} > $bestSc){
						if ($compl{$lBin} > 60 && $conta{$lBin} < 20){
							$bestSc = $score{$lBin}; $bestBin = $lBin;
						} else {
							$rejQuali++;
						}
					}
				}
				my $MGSfam = $ff . "." .$MGS2bin{$bestBin};
				#print "$MGSfam     $bestBin  $bestSc $compl{$bestBin} $conta{$bestBin} \n";
				$ret{$MGSfam} = $bestBin;
			}
			
			#reset params
			$lastMGS = $curMGS;
			%conta = (); %compl = (); %fam=(); %score = (); %MGS2bin = ();
		}
		
		if ($spl[0] !~ m/^Cano__/){
			#$ret{$curMGS} = $spl[0] ;
			$conta{$curBin} = $spl[$ContaIdx];$compl{$curBin} = $spl[$ComplIdx];
			$score{$curBin} = $spl[$ComplIdx] - (2. * $spl[$ContaIdx]);
			$spl[0] =~ m/^(.+)__(.+)$/;my $cFam = $famSmpl{$1};
			push(@{$fam{$cFam}}, $curBin);
			$MGS2bin{$spl[0]} = $curMGS;
		}	
	}
	close $I;
	
	print "Found representative for ". scalar(keys %ret). " MGS, $rejQuali rejected on Quali\n";
	print "Identified $famFnd families, $assmGrpFnd assembly groups, $noGrpFnd fallbacks\n";
	return \%ret;
}

#extract 1 reference genome per MGS, choosing the default rep and extracting its contigs from original assembly file
sub createBinCtgs{
	#$binDctg,$hrM,"$logDir/MAGvsGC.txt.gz
	my ($outD,$hrMap,$guideF,$perFam) = @_;
	
	my $hr;
	if ($perFam){
		print "Getting per family ref genomes\n";
		$hr = getRepresentBinsPerFamily($guideF,$hrMap);
	} else {
		$hr = getRepresentBins($guideF);
	}
	my %repBins = %{$hr};
	
	my %map = %{$hrMap};
	my @allReps =  sort { $repBins{$a} cmp $repBins{$b} } keys %repBins; #@allReps = sort @allReps;
	#die "@allReps\n";
	my $hr2;my$hr1; my $lastSmpl = "";
	foreach my $MGS (@allReps){
		my $MAG = $repBins{$MGS};
		if ($MAG =~ m/^Cano__/){
			print "Could not retrieve MAG for $MGS : $MAG , because is Canopy\n";
			next;
		}
		my $smpl = ""; my $bin= "";
		if ($MAG =~ m/^(.+)__(.+)$/){  $smpl = $1;  $bin=$2;
		} else { #skip this MAG completely.. but not good
			print "Could not match \"$MAG\" to sample and contig!\n";
			next;
		}
			
		
		$smpl = $map{altNms}{$smpl} if ( defined($map{altNms}{$smpl}) );
		if ($lastSmpl ne $smpl){
			$lastSmpl = $smpl;
			print "Reading $smpl\n";
			my $dirIn = $map{$smpl}{wrdir}; 
			my $assDir = getAssemblPath($dirIn);
			my $BinDir = "$assDir/Binning/SB/"; my $BinFile = "$BinDir/$smpl";
			$hr1 = readBinSB($BinFile);
			$hr2 = readFasta("$assDir/scaffolds.fasta.filt");
		}
		
		
		my @ctgs = @{${$hr1}{$bin}};
		#die "@ctgs\n". @ctgs . "\n";
		
		my $outF = "$outD/$MGS.ctgs.$MAG.fna";
		#print "writing  $MGS.ctgs.$MAG.fna\n";
		open O,">$outF" or die "Couldn't open $outF\n";
		foreach my $ctg (@ctgs){
			print O ">$ctg\n${$hr2}{$ctg}\n";
		}
		close O;
		#die "$MAG :: $smpl $bin\n$dirIn\n$assDir\n$BinFile\n";
	}

	print "----------------------\nDone\nWrote representative MAGs (contigs) to $outD\n----------------------\n";
}

sub createBin2{
	my ($binD,$cnopyF,$refFA) = @_;
	my $hr = readMGSrevRed($cnopyF);
	my %G2MGS = %{$hr};
	my ($I,$OK) = gzipopen($refFA,"reference gene cat",1);
	my $seq=""; my $hd=""; my %MGSfxa;
	my $MGScnt=0; my $geneCnt=0;
	my $fileEnd = ".fna"; $fileEnd = ".faa" if ($refFA =~ m/\.faa$/);
	while (my $line = <$I>){
		chomp $line;
		if ($line =~ m/^>/){
			#take care of old fna..
			if ($hd =~ m/^>(\d+)/ && exists($G2MGS{$1})){
				$geneCnt++;
				my @tars = keys %{$G2MGS{$1}};
				foreach my $MGS (@tars){
					$MGSfxa{$MGS}{$hd} = $seq;
				}
			}
			$hd = $line; $seq = "";  
			next;
		}
		$seq .= $line;
	}
	close $I;
	print "Found $geneCnt genes in " . scalar(keys%MGSfxa). " MGS (avg " . int($geneCnt/scalar(keys%MGSfxa)*100)/100  . " per MGS). Writing to $binD\n";
	system "mkdir -p $binD" unless (-d $binD);
	foreach my $MGS (keys %MGSfxa){
		open O,">$binD/$MGS$fileEnd";
		foreach my $gen (keys %{$MGSfxa{$MGS}}){
			print O "$gen\n$MGSfxa{$MGS}{$gen}\n";
		}
		close O;
	}
	
}

sub createBinFAA{
	my ($binD,$cnopyF,$refFA) = @_;
	my $suffix = "faa";
	$suffix = $_[3] if (@_ > 3);
	print "Reading reference MGS $cnopyF\n";
	my $hr = readMGS($cnopyF);my %clust = %{$hr};
	print "Reading ref FAA $refFA\n";
	$hr = readFasta($refFA,1);
	my %FAA = %{$hr};
	
	system "mkdir -p $binD" unless (-d $binD);
	foreach my $cl (sort keys %clust){
		my $oF = "$binD/$cl.$suffix";
		my @refG = @{$clust{$cl}};
		my %alreadySeen;
		my $ostr = ""; my $gcnt=0;
		foreach my $rg (@refG){
			unless(exists($FAA{$rg})){
				print "Can't find gene $rg in gene cat\n" ;
				next;
			}
			my $rn = $rg; $rn =~ s/://;
			next if (exists($alreadySeen{$rn}));
			$ostr .= ">$rn\n$FAA{$rg}\n";
			$alreadySeen{$rn} = 1;
			$gcnt ++;
		}
		if ($ostr ne ""){
			open O,">$oF" or die $!;
			print O $ostr;
			close O;
		}
	}
}



sub runCheckM{#runs checkM on *.faa files (each file one Bin)
	my ($binD,$outFile,$tmpD,$ncore) = ($_[0],$_[1],$_[2],$_[3]);
	my $runNow = 1; $runNow = $_[4] if (@_ > 4);
	my $ext = "faa"; $ext = $_[5] if (@_ > 5);
	my $gtag = "--genes"; $gtag = "" if ($ext eq "fna");
	#system "rm -rf $tmpD/CM/";
	#system "mkdir -p $tmpD/tmp/" unless(-d "$tmpD/tmp/");
	my $cmC = "";
	$cmC .= "rm -rf $tmpD/CM/;mkdir -p $tmpD/tmp/\n";
	#my $p2a = getProgPaths("py2activate");
	#my $pd = getProgPaths("pydeacti");
	my $checkMBin = getProgPaths("checkm");
	#$cmC .= "$p2a\n";
	$cmC .= "$checkMBin lineage_wf $gtag -x $ext -t $ncore --tab_table -f $outFile -q --pplacer_threads 3 --tmpdir $tmpD/tmp/ $binD $tmpD/CM/\n";
	#$cmC .= "$pd\n";
	$cmC .= "rm -rf $tmpD/CM/ $tmpD/tmp/\n";
	
	#first check there's something in the bin file..
	my $outFile2 = $outFile; $outFile2 =~ s/\.cm$//;
	if (-e $outFile2){
		open I,"<$outFile2" or die "No input for runCheckM function:$outFile\n";
		my %binsFnd;
		while (my $l =<I>){
			chomp $l; my @spl = split /\t/,$l;
			$binsFnd{$spl[1]} = 1;
		}
		close I;
		my @bins = keys %binsFnd;
		print "Found ".(scalar( @bins) - 1)  . " metag Bins\n";
		if (@bins <= 1){$cmC = "\ntouch $outFile\n";}
	}
	print "$cmC\n";
	
	systemW "$cmC" if ($runNow > 0);
	return $cmC;
}

sub runCheckM2{#runs checkM2 on *.faa files (each file one Bin)
	my ($binD,$outFile,$tmpD,$ncore) = ($_[0],$_[1],$_[2],$_[3]);
	my $runNow = 1; $runNow = $_[4] if (@_ > 4);
	my $ext = "faa"; $ext = $_[5] if (@_ > 5);
	$tmpD .= "/CM2";
	my $gtag = "--genes"; $gtag = "" if ($ext eq "fna");
	#system "rm -rf $tmpD/CM/";
	#system "mkdir -p $tmpD/tmp/" unless(-d "$tmpD/tmp/");
	my $cmC = "";
	$cmC .= "rm -rf $tmpD/;  mkdir -p $tmpD/\n";
	my $outD = $outFile; $outD =~ s/\/[^\/]+$/\//; $outD .= "/CHM2/";
	#my $pd = getProgPaths("pydeacti");
	#my $condaA = getProgPaths("CONDA");
	my $checkM2Bin = getProgPaths("checkm2");
	my $chm2DB = getProgPaths("checkm2DB");
	
	#$cmC .= "$condaA\nconda activate checkm2\n";
	if ($checkM2Bin =~ m/ activate /){
		$checkM2Bin =~ s/activate (\S+)/activate $1\nexport CHECKM2DB=$chm2DB\n/;
	}
	
	
	$cmC .= "$checkM2Bin predict --input $binD $gtag -x $ext --force -t $ncore  --output-directory $tmpD \n";#--tmpdir $tmpD/tmp/ $binD $tmpD/CM/\n";
	#$cmC .= "$pd\n";
	#debugging only	
#	$cmC .= "mkdir -p $outD;cp -r $tmpD $outD\n"; #replace with more targeted function later
	$cmC .= "cp $tmpD/quality_report.tsv $outFile\n"; #replace with more targeted function later
	$cmC .= "touch $outFile\n";
	$cmC .= "rm -rf $tmpD\n";
	
	if ($runNow > 0){
		print "running checkM2 local..\n";
		systemW "$cmC" ;
	}
	return $cmC;
}



sub runSemiBin{
	my ($jgO,$outDir, $tmpDir, $nm, $fna, $cores, $dirsAR, $seqTec, $giveSBenv ) = @_;
	#human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global
	my $SBbin = getProgPaths("SemiBin2");
	#my $semibinGTDB = getProgPaths("semibinGTDB");
	#get list of bams/crams..
	my @dirSS = @{$dirsAR}; #= split(',',$dirs);
	#go through each dir and find sample name
	my $comBAM = "";
	my $isCram=0; my $numBams=0;
	my @BAMS;  
	my $uncramCmd = "";
	foreach my $DDI (@dirSS){
		$numBams++; my %iBAMS;
		my $iBAM = $DDI;
		if ( $DDI =~ m/\/$/ ||  $DDI !~ m/bam$/ ){
			unless (-e "$DDI/mapping/done.sto"){print "Can't find $DDI/mapping/done.sto!! \n Aborting SemiBin\n"; return "" ;}
			my $SmplNm = `cat $DDI/mapping/done.sto`;#$SmplNm =~ s/-smd.bam\n?//;
			chomp $SmplNm;	my $tbam = "$DDI/mapping/$SmplNm";
			if (!-e $tbam){$isCram=1;$tbam =~ s/\.bam/\.cram/;}
			unless (-e $tbam){
				print "runSemiBin:::Can't find either bam nor cram at $DDI\nAborting Binning prep for current sample group $nm\n";
				return "";
			}
			$iBAM = $tbam;
		} 
		next if (-s $iBAM < 15000000 ); #cram/bam too small.. prob no good   (15000000)
		my $oBAM = "$tmpDir/$nm.$numBams.bam";
		$uncramCmd .= cram2bsam($iBAM,$fna,$oBAM,1,$cores) if ($isCram && ! -e $oBAM);
		push @BAMS, $oBAM;
		#check explicitly for suppl mapping being present..
		if ($iBAM =~ m/\.sup/){
			$iBAM =~ s/\.sup//;
		} else {
			$iBAM =~ s/-smd/\.sup-smd/;
		}
		if (-e $iBAM){
			$oBAM = "$tmpDir/$nm.$numBams.sup.bam";
			$uncramCmd .= cram2bsam($iBAM,$fna,$oBAM,1,$cores) if ($isCram && ! -e $oBAM);
			push @BAMS, $oBAM;
		}
	}
	if (@BAMS == 0){
		#fake empty entry...
		print "runSemiBin::No bams found, creating fake output\n";
		system "mkdir -p $outDir;touch $outDir/$nm;touch $outDir/$nm .assStat";
		open O,">$outDir/$nm.cm2";
		print O "Name\tCompleteness\tContamination\tCompleteness_Model_Used Translation_Table_Used\tAdditional_Notes\n";
		close O;
		return "" ;
	}
#die;
	
	
	# --environment human_gut, dog_gut, ocean, soil, cat_gut, human_oral, mouse_gut, pig_gut, built_environment, wastewater, chicken_caecum, global
	my $smode = "single_easy_bin ";
	my $senv = "--environment human_gut";$senv =  "--environment $giveSBenv" if ($giveSBenv ne "");
	my $dflags = " --random-seed 555 --tmpdir $tmpDir -p $cores";
	my $seqType = "--sequencing-type=short_read ";
	$seqType = "--sequencing-type=long_read " if ($seqTec eq "PB" || $seqTec eq "ONT" || $seqTec eq "hybrid");#PAcBIo/ONT
	my $cmd1 = "###preparing BAMs..\n$uncramCmd\n\n";
	$cmd1 .= "echo \"CRAM->BAM finished\"\n";
#	my $cmd = "";
	#my $output = "$outD/$nm.semibin";
	my $cmd .= "\n\n###Running SemiBin...\n";
	if ($giveSBenv ne ""){
		$cmd .= "$SBbin $smode -i $fna -b  ". join(" ",@BAMS). " $seqType --output $outDir $dflags\n";

	}elsif ($numBams == 1 && $jgO ne ""){
		$cmd .= "$SBbin $smode --depth-metabat2 $jgO -i $fna $senv $seqType --output $outDir $dflags\n";
	} else {
		#--reference-db-data-dir $semibinGTDB --training-type self 
		$cmd .= "$SBbin $smode -i $fna -b  ". join(" ",@BAMS). " $seqType --output $outDir $dflags\n";
	}
	
	#die $cmd;
	
	return $cmd1.$cmd;

}


sub runMetaDecoder{
	my ($jgO,$outDir, $tmpDir, $nm, $fna, $cores, $dirsAR ) = @_;
	#human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global
	my $MDbin = getProgPaths("MetaDecoder");
	my $baseN = "$tmpDir/$nm";
	#get list of bams/crams..
	my @dirSS = @{$dirsAR}; #= split(',',$dirs);
	#go through each dir and find sample name
	my $comBAM = "";
	my $isCram=0; my $numBams=0;
	my @SAMS;  
	my $uncramCmd = "";
	foreach my $DDI (@dirSS){
		$numBams++;
		my $oSAM = "$tmpDir/$nm.$numBams.sam";
		my $iBAM = $DDI;
		if ( $DDI =~ m/\/$/ ||  $DDI !~ m/bam$/ ){
			my $SmplNm = `cat $DDI/mapping/done.sto`;#$SmplNm =~ s/-smd.bam\n?//;
			chomp $SmplNm;	my $tbam = "$DDI/mapping/$SmplNm";
			if (!-e $tbam){$isCram=1;$tbam =~ s/\.bam/\.cram/;}
			die "runMetaDecoder:::Can't find either bam nor cram at $DDI\n" unless (-e $tbam);
			$iBAM = $tbam;
		} 
		next if (-s $iBAM < 3000);
		$uncramCmd .= cram2bsam($iBAM,$fna,$oSAM,2,$cores) unless (-e $oSAM);
		push @SAMS, $oSAM;
	}
	if (@SAMS == 0){
		#fake empty entry...
		system "mkdir -p $outDir;touch $outDir/$nm;touch $outDir/$nm .assStat";
		open O,">$outDir/$nm.cm2";
		print O "Name\tCompleteness\tContamination\tCompleteness_Model_Used Translation_Table_Used\tAdditional_Notes\n";
		close O;
		return "" ;
	}
	
	# --environment human_gut, dog_gut, ocean, soil, cat_gut, human_oral, mouse_gut, pig_gut, built_environment, wastewater, chicken_caecum, global
	my $cmd = "###preparing SAMs..\n$uncramCmd\n\n";
	$cmd .= "$MDbin coverage -s ". join(" ",@SAMS). " -o $baseN.coverage --threads $cores --mapq 10\n";
	$cmd .= "$MDbin seed  --threads $cores -f $fna -o $baseN.SEED\n";
	$cmd .= "$MDbin cluster -f $fna -c $baseN.coverage -s $baseN.SEED -o $baseN --random_number 511 --no_clusters\n"; 
	$cmd .= "mkdir -p $outDir;\n \ncp $baseN.cluster $outDir/$nm;\n";
	
	#still need to transfer files from here..
	#die $cmd;
	
	return $cmd;

}

sub runMetaBat{
	my ($jgO,$outD,$nm, $fna ) = @_;
	if (-e "$outD/MeBa.sto"){
		#print "metabat2 done: $outD\n";
		return "";
	}
	my $numC = 0;
	$numC = $_[4] if (@_ > 4);
	#print "Running MetaBat..\n";
	my $metab2Bin = getProgPaths("metabat2");
	#my $mbBin = "/g/bork3/home/hildebra/bin/metabat/./metabat";
	my $outCtgNms = "$outD/$nm.ctgs.txt";
	my $outFna = "$outD/$nm.fasta.fna";
	my $outD2 = "$outD/$nm";
	my $outMat = "$outD/$nm.mat.txt";
	#my $jgO = "$tmp/depth.jgi";
	if (!-e $fna){die "Can't find requried scaffold file in metabat routine: $fna\n"; }
	#start metabat
	my $cmd = "";#$before."\n";
	#$cmd .= "$mbBin -i $fna -a $jgO.depth.txt -o $outFna -p $jgO.pairs.sparse -l $outCtgNms --minCVSum  10 -t $numCore --saveCls $outMat -v\n";
	$cmd .= "mkdir -p $outD\n";
	$cmd .= "$metab2Bin -i $fna -a $jgO -o $outD2 -l --saveCls -t $numC --noBinOut -m 2500 --minCVSum 1\n"; #$outMat
	$cmd .= "echo \"$nm\" > $outD/MeBa.sto\n";
	return $cmd;
	#die $cmd."\n";
	#my $jobName = "MBat";
	#my ($jobNameX, $tmpCmd) = qsubSystem($logDir."metaBat.sh",$cmd,$numCore,"60G",$jobName,"","",1,[],$QSBoptHR);

}

