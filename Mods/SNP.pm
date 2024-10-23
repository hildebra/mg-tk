package Mods::SNP;
use warnings;
use strict;

use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( gzipopen systemW readFasta readGFF writeFasta reverse_complement_IUPAC );
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystem2);

use List::Util qw/shuffle/;


use Exporter qw(import);
our @EXPORT_OK = qw(SNPconsensus_vcf SNPconsensus_fasta SNPconsensus_vcf2);

sub regionsFromFAI($){
	my ($inF ) =@_;
	my @ret;
	open I,"<$inF" or die "can't open fai $inF\n";
	while ( my $line = <I>){
		chomp $line;
		my @fields = split /\t/,$line;
		push(@ret,$fields[0] . ":0-" . $fields[1]);
	}
	close I;
	return (@ret);
}

sub getRegionsBamDepth{
	my ($depthPC,$totalSpl,$maxSNPcores) = @_;
	my %LperC; my %DperC; my $tDep=0; my %contigNum; my $tLen=0;
	my $cnt=0;
	#print "$depthPC\n";
	#open I ,"<$depthPC" or die "can;t oopen depth file $depthPC\n";
	my ($IN ,$status) = gzipopen($depthPC,"contig depth file",0);
	while (<$IN>){
		chomp; my @spl = split /\t/;
		
		if ($spl[0] =~ m/__C(\d+)_L=(\d+)=/){
			$contigNum{$1} = $spl[0];
			$LperC{$spl[0]} = $2;
			$tLen += $2;
			$tDep += $2 * $spl[1] ;
		} else {
			$contigNum{$cnt} = $spl[0];
			$LperC{$spl[0]} = 1000;
			$tLen += 1000;
			$tDep += 1000 * $spl[1] ;
		}
		$DperC{$spl[0]} = $spl[1];
		$cnt++;
	}
	close $IN;
	#some hardcoded rules to make small samples into smaller job numbers..
	if ($tDep ==0 || $tLen ==0){
		return [],[];
	}
	
	if ($tDep <1e6){
		$totalSpl = 1;
	} elsif ($tDep <5e6){
		$totalSpl = int($maxSNPcores/3);
	} elsif ($tDep <20e6){
		$totalSpl = $maxSNPcores/2;
	} else {$totalSpl = $maxSNPcores/2;}
	#print "totalSpl $totalSpl $maxSNPcores $tDep\n";
	#expected depth per bin
	my $exD = $tDep/$totalSpl;
	#print "$tDep/$tLen\n";
	my $avgD = $tDep/$tLen;
	#now count up per bin to get to this number..
	my $curD = 0; my @regions; $cnt =0;
	my @regOrd;
	my @srtK = sort {$a <=> $b} keys %contigNum;
	my $startCtg =0;
	my $idx=0;
	foreach my $id (@srtK){
		push (@regOrd,$contigNum{$id});
	}
	while ( $idx<@srtK ){
		#print "$idx ";
		my $k1 = $srtK[$idx];
		my $k = $contigNum{$k1};
		my $thisD = $DperC{$k} * ($LperC{$k} - $startCtg);
		#print "$idx ".@srtK." $curD+$thisD < $exD $startCtg\n";
		
		if ($DperC{$k} > $avgD*0.85){ #deep sample, uses more processing power..
			$thisD *=  1.2;
		} 
		if ($DperC{$k} < $avgD*0.1){ #deep sample, uses more processing power..
			$thisD *=  0.8;
		}
		#if ($k =~ m/MM2__C152_L=33826=/){die "$curD + $thisD) < $exD \n";}
		if ( ($curD + $thisD) > $exD ){
			#die "$curD + $thisD) < $exD \n";
			if ($LperC{$k} < 1000){
				$regions[$cnt] .= "$k\t0\t$LperC{$k}\n";
				$curD = 0; $cnt ++; $startCtg = 0; $idx++;
			} else {
				#how much into the contig?
				my $stopCtg = int($LperC{$k} * (($exD - $curD) / $thisD) + $startCtg);
				if ($stopCtg < 150){#just don't take this, reset counter..
				#print "A";
					#$curD = $thisD; $cnt ++;$regions[$cnt] .= "$k\t0\t$LperC{$k}\n";
					$startCtg = 0;$curD =0; $idx ++; 
				} elsif ( ($LperC{$k} - $stopCtg) < 150){ #just take whole...
				#print "B";
					$regions[$cnt] .= "$k\t$startCtg\t$LperC{$k}\n"; $curD = 0; $cnt ++; $startCtg = 0; $idx++;
				} else { #take part of contig, reset counter..
				#print "C";
					#die "$stopCtg : $LperC{$k} * ($exD - $curD) / $thisD;\n$k\t0\t$stopCtg\n$k\t$stopCtg\t$LperC{$k}\n";
					$regions[$cnt] .= "$k\t$startCtg\t$stopCtg\n";$startCtg = $stopCtg;
					$curD = 0; $cnt ++; 
				}
			}
			#reset
		} else {
			$regions[$cnt] .= "$k\t$startCtg\t$LperC{$k}\n"; $idx++;
			$curD += $thisD; $startCtg = 0;
		}
	}
	#die @regions." regions\n";
	return (\@regions,\@regOrd);
}

sub getRegionsBam{
	my ($splitFAsizeL,$refFA,$tmpD) = @_;
	
	my $smtBin = getProgPaths("samtools");
	my $py_genRegs = getProgPaths("genRegions_scr"); #python $frDir/fasta_generate_regions.py
	my $regionFile = "$refFA.reg";
	systemW "$smtBin faidx $refFA" unless (-e "$refFA.fai");
	my @curReg = regionsFromFAI("$refFA.fai");
	#die "@curReg\n";
	#if (!-e $regionFile ){		systemW "python $py_genRegs $refFA.fai $splitFAsizeL > $regionFile\n"	}
	#open my $handle, '<', $regionFile;	chomp(@curReg = <$handle>);	close $handle;
	
	#now comes the real work: translate to bed format
	#@curReg = shuffle @curReg;
	my @regions; my @regOrd;
	my $curCnt = 0; my $cnt=0;
	foreach my $reg (@curReg){
		my @spl = split(/:/,$reg);
		push(@regOrd,$spl[0]);
		my @spl2 = split(/-/,$spl[1]);
		my $sum =  $spl2[1]-$spl2[0];
		if ($curCnt + $sum > $splitFAsizeL){
			#how much can still be added?
			my $max = $curCnt + $sum - $splitFAsizeL;
			my $start = $spl2[0];
			if ( (($spl2[0]+$max-1)) + 100 > $spl2[1]){
				#print $spl2[0]+$max-1 .":".$spl2[1]."\n";
				$regions[$cnt] .= "$spl[0]\t$spl2[0]\t$spl2[1]\n";
				$cnt++;
				$curCnt=0; 
			} else {
				$start = ($spl2[0]+$max);
				$regions[$cnt] .= "$spl[0]\t$spl2[0]\t". ($start) ."\n";
				$cnt++;
				$regions[$cnt] .= "$spl[0]\t".($start)."\t".$spl2[1]."\n";
				$curCnt = $spl2[1] - $start;
			}
			#print "$regions{$cnt}\n";
			#die "new\n$regions{$cnt}\n";
		} else {
			$regions[$cnt] .= "$spl[0]\t$spl2[0]\t$spl2[1]\n";
			$curCnt += $sum;
			
		}
	}
#	die "@regions\n";
	return (\@regions,\@regOrd);
}


sub pileupcall{
	
	my ($tarR,$tag,$SNPIHR,$QSBoptHR,$scrDir,$tmpOut,$myParL,$crAR) = @_;
	
	my @curReg = @{$crAR};
	
	#die "@curReg\n";
	my $frbBin = getProgPaths("freebayes");
	my $bcftBin = getProgPaths("bcftools");
	my $refFA = $SNPIHR->{assembly};
	my $tmpdir = $SNPIHR->{nodeTmpD};
	my $smplNm = $SNPIHR->{smpl};
	my $qsubDirE = $SNPIHR->{qsubDir};
	my $runLocalTmp = $SNPIHR->{runLocal};
	my $x = $SNPIHR->{JNUM};
	my $threads = 0;#$SNPIHR->{threads};
	my $run2ctg = $SNPIHR->{run2ctg};
	my $rdep = $SNPIHR->{rdep} ;


	my $useFB = 1;	$useFB = 0 if (uc($SNPIHR->{SNPcaller}) eq "MPI");
	my $overwrite = $SNPIHR->{overwrite};

	#basic caller options..
	my $minBQ=30; my $minMQ=30;
	#freebayes std options
	my $frAllOpts= "-u -i -m $minMQ -q $minBQ -C 1 -F 0.1 -k -X --pooled-continuous --report-monomorphic  --min-repeat-entropy 1 --use-best-n-alleles 2 -G 1 ";
	#bcftools options #-q = map qual -Q = base qual
	my $bcfAllOpts = "--count-orphans --min-BQ $minBQ -d 12000 --threads $threads --skip-indels --min-MQ $minMQ -a DP,AD,ADF,ADR,SP"; #Pernille
	
	if ($tag eq "sup-"){
		if ($SNPIHR->{SeqTechSuppl} eq "ONT"){$bcfAllOpts.=" -X ont ";
		} elsif ($SNPIHR->{SeqTechSuppl} eq "PB"){$bcfAllOpts.=" -X pacbio-ccs ";
		} else{$bcfAllOpts.=" -X illumina ";}
	} else {
		if ($SNPIHR->{SeqTech} eq "ONT"){$bcfAllOpts.=" -X ont ";
		} elsif ($SNPIHR->{SeqTech} eq "PB"){$bcfAllOpts.=" -X pacbio-ccs ";
		} else{$bcfAllOpts.=" -X illumina ";}
	}
	
	my $cmd = ""; 
	my $locXtrCmd = ""; $locXtrCmd = " &" if ($runLocalTmp);
	#my $tag = "primary";
	#
	#$cmd .= "mkdir -p $tmpdir\n";
	my $cmdAll2 = "";
	$cmdAll2 .= "echo \"Processing bams - mpileup $tag\"\n";
	if ($useFB){
		$cmd = "ulimit -s unlimited\n$frbBin -f $refFA  $frAllOpts ";
	} else {
		$cmd = "$bcftBin mpileup --fasta-ref $refFA $bcfAllOpts ";
	}
	my @allDeps2; my @checkF;
	#implement in parallel as too slow in single core mode :/
	system "rm -f $tmpOut.$tag*" if ($overwrite);
	system "mkdir -p $qsubDirE" unless (-d $qsubDirE);
	my $bedJobs = 0;
	for (my $i=0;$i<@curReg;$i++){ #go over regions in bed file, submit a job for each "region"
		#$tar[0] = bam file;  $bedF = bedfile with regions
		next if (!$run2ctg);
		system "rm -f $qsubDirE/$smplNm.$tag*.bed" if ($i==0) ;
		my $cmd2 = $cmd ;
		my $bedF = $qsubDirE."$smplNm.${tag}$i.bed";
		push @checkF, $bedF;
		next if (!-e $bedF && -e "$tmpOut.$tag$i" && !$overwrite);
		if ($myParL){
			open O,">",$bedF or die $!;print O $curReg[$i];close O;
			if ($useFB){
				$cmd2 .= " -t $bedF $tarR->[0] > $tmpOut.$tag$i && rm $bedF $locXtrCmd\n"; #--region '$curReg[$i]'
			} else {
				$cmd2 .= " -R $bedF $tarR->[0] | $bcftBin call --output-type v --ploidy 1 --multiallelic-caller -M --output-type v | bgzip -c > $tmpOut.$tag$i.gz  && rm $bedF $locXtrCmd\n"; 
			}
		} else {
			die "incomplete control structure SNP.pm\n";
		}
		$cmdAll2 .= $cmd2."\n" if ($run2ctg);
		if (!$runLocalTmp){
			my ($dep,$qcmd) = qsubSystem($qsubDirE. $SNPIHR->{cmdFileTag} . ".ac$tag.$smplNm.$i.sh",$cmd2,1,"15G","FBC$x.$i",$rdep,"",1,[],$QSBoptHR);
			push (@allDeps2,$dep);
		}
		$bedJobs++;
		#last if ($i == 1000);
	}
	#$bedJobs =1 if ($bedJobs<1);
	$cmdAll2 .= "wait \$(jobs -p);\n" if ($run2ctg);
	$cmdAll2 .= "rm -f $tarR->[0].crai $tarR->[0].bai;\n";
	$cmdAll2 .= "\necho \"Finished mpileup $tag\"\n\n";	
	#$cmdAll .= "rm -f $tarR->[0];\n" if ($run2ctg && $tarR->[0] ne "");
	if ($bedJobs ==0){$cmdAll2="";}
	return (\@allDeps2, \@checkF, $cmdAll2);
}


sub SNPconsensus_vcf{
	my ($SNPIHR)  = @_;
	#my $kpathBin = getProgPaths("kpath");
	my $smtBin = getProgPaths("samtools");
	my $vcfcnsScr = getProgPaths("vcfCons_scr");
	#my $tabixBin = getProgPaths("tabix");
	#my $vcfLD = getProgPaths("vcfLib_dir");#"/g/bork3/home/hildebra/bin/vcflib/bin/";
	my $ctg2fas = getProgPaths("contig2fast_scr");
	my $pigzBin  = getProgPaths("pigz");
	#my $py3 = getProgPaths("py3activate",0);my $py3d = getProgPaths("pydeacti",0);
	my $bcftBin = getProgPaths("bcftools");

	#get parameteres
	my $samcores = 12;
	my %SNPinfo = %{$SNPIHR};
	my $QSBoptHR = $SNPIHR->{QSHR};
	my $x = $SNPIHR->{JNUM};
	my $jdep = ""; $jdep = $SNPIHR->{jdeps} if (exists($SNPIHR->{jdeps}));
	my $tmpdir = $SNPIHR->{nodeTmpD};
	my $smplNm = $SNPIHR->{smpl};
	my $refFA = $SNPIHR->{assembly};
	my $qsubDirE = $SNPIHR->{qsubDir};
	my $scrDir = $SNPIHR->{scratch};
	my $bamcram = $SNPIHR->{bamcram};
	my $splitFAsize = $SNPIHR->{bpSplit};
	my $overwrite = $SNPIHR->{overwrite};
	my $runLocalTmp = $SNPIHR->{runLocal};
	my $maxSNPcores= $SNPIHR->{maxCores};
	my $SNPstone = $SNPIHR->{STOconSNP}; my $SNPsuppStone = $SNPIHR->{STOconSNPsupp};
	my $minDepth = 0;
	$minDepth =$SNPIHR->{minDepth}  if (exists($SNPIHR->{minDepth} ));
	my $minCallQual = 20;
	#my $SNPstone = $ofasConsDir."SNP.cons.stone";
	#my $memReq = "20G";
	my $memReq = $SNPIHR->{memReq};
	my $vcfFile = ""; $vcfFile = $SNPIHR->{vcfFile} if (exists ($SNPIHR->{vcfFile}));
	my $vcfFileS = ""; $vcfFileS = $SNPIHR->{vcfFileSupp} if (exists ($SNPIHR->{vcfFileSupp}));
	my $cmdFTag = $SNPIHR->{cmdFileTag};
	my $firstInSample = 0;$firstInSample = $SNPIHR->{firstInSample} if (exists($SNPIHR->{firstInSample}));
	my $useFB = 1;	$useFB = 0 if (uc($SNPIHR->{SNPcaller}) eq "MPI");
	$vcfcnsScr = getProgPaths("vcfCons_FB_scr") if ($useFB);
	if ($runLocalTmp){
		$scrDir = $tmpdir;
		$samcores = $maxSNPcores;#$SNPIHR->{split_jobs};
	}


	my $ofasCons = $SNPIHR->{ofas};
	$ofasCons =~ m/(^.*)\/[^\/]+$/;
	my $ofasConsDir = $1."/";
	#$ofasCons .= ".gz" unless ($ofasCons =~ m/\.gz$/); #don't change, needed without..
	my $run2ctg=1; #flag to determine if I run the cram to bam, mpileup, consensus contig steps..
	system "rm -f $ofasConsDir/*" if ($overwrite);
	#die "$ofasCons\n";
	if (-e "$ofasCons.gz" && -s "$ofasCons.gz" > 200){
		$run2ctg =0 ;
	} else {
		system "rm -f $ofasConsDir/*"; #better safe than sorry..
	}
	#first all important regions on finalDir
	my @curReg = ("1"); my @regOrd;
	my $myParL=0;
	if ($splitFAsize>0){$myParL=1;}
	if ($myParL && $run2ctg){ #no, don't redo freebayes part
		my ($refAR,$refAR2);
		if (exists($SNPIHR->{depthF}) && $SNPIHR->{depthF} ne ""){
			 ($refAR,$refAR2) = getRegionsBamDepth($SNPIHR->{depthF},$SNPIHR->{split_jobs},$maxSNPcores);
		} else {
			 ($refAR,$refAR2) = getRegionsBam($splitFAsize,$refFA,$tmpdir);
		}
		@curReg = @{$refAR}; 
		@regOrd = @{$refAR2};  #use .fai instead for this..
		if (@curReg == 0){return($SNPIHR,"");}
		
		#open O,">$refFA.reg" or die "can't open region file $refFA.reg\n";		print O join("\n",@regOrd);		close O;
		#die "$refFA.reg\n@curReg\n";
	}
	

	
	my $rdep="";
	#prepare files..
	my $cleanCmd = ""; 
	my $xtra = "";
	$xtra .= "echo \"Preparing data\"\n";
	$xtra .= "mkdir -p $scrDir;\n";
	$xtra .= "$smtBin faidx $refFA;\n" unless (-e "$refFA.fai");
	#$xtra .= "cp $refFA $refFA.fai $scrDir;\n";$refFA =~ m/\/([^\/]+$)/;$refFA = "$scrDir/$1";
	#my $preTar = 
	
	my @tar = ("");
	$xtra .= "echo \"Creating c/bams indexes primary reads\"\n";
	$tar[0] = ${$SNPIHR->{MAR}}[0]; #$preTar;
	if ($bamcram eq "cram"){ #create index for bam/cram
		$xtra .= "if [ ! -e $tar[0].crai ] || [ ! -s $tar[0].crai ]; then rm -f $tar[0].crai; $smtBin index -@ $samcores  $tar[0]; fi\n";
	} else {
		$xtra .= "if [ ! -e $tar[0].bai ] || [ ! -s $tar[0].bai ]; then rm -f $tar[0].bai; $smtBin index -@ $samcores  $tar[0]; fi\n";
	}
	
	if (!$runLocalTmp && $run2ctg && (!-e $tar[0] || !-e $refFA) ){
		my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.CramToBam$x.sh",$xtra,2,"17G","CtB$x",$jdep,"",$samcores,[],$QSBoptHR);
		$cleanCmd .= "rm -r $scrDir\n";
		$rdep = $dep;
		$xtra = "";
	}
	$SNPIHR->{run2ctg} = $run2ctg;
	$SNPIHR->{rdep} = $rdep;

	#$SNPIHR->{assembly} = $refFA;
	my $cmdAll = ""; $cmdAll .= $xtra if ($run2ctg);
	my $tmpOut = "$scrDir/$smplNm.cons.vcf";
	my ($dAR,$cAR,$pilecmd) =  pileupcall(\@tar,"",$SNPIHR,$QSBoptHR,$scrDir,$tmpOut,$myParL,\@curReg);
	my @allDeps2 = @{$dAR}; my @checkF = @{$cAR};
	$cmdAll .= $pilecmd;
	
	
	
		#supplementary mappings?
	my @tarS = ("");
	my $tmpOut2 = "$scrDir/$smplNm.X.cons.vcf";
	if ($SNPsuppStone ne "" && exists ($SNPIHR->{MARsupp} ) ){
		my $xtra2 .= "echo \"Creating c/bams indexes supplemental reads\"\n";
		$tarS[0] = ${$SNPIHR->{MARsupp}}[0];
		if ($bamcram eq "cram"){ #create index for bam/cram
			$xtra2 .= "if [ ! -e $tarS[0].crai ] || [ ! -s $tarS[0].crai ]; then rm -f $tarS[0].crai; $smtBin index -@ $samcores  $tarS[0]; fi\n";
		} else {
			$xtra2 .= "if [ ! -e $tarS[0].bai ] || [ ! -s $tarS[0].bai ]; then rm -f $tarS[0].bai; $smtBin index -@ $samcores  $tarS[0]; fi\n";
		}
		($dAR,$cAR,$pilecmd) =  pileupcall(\@tarS,"sup-",$SNPIHR,$QSBoptHR,$scrDir,$tmpOut2,$myParL,\@curReg);
		$cmdAll .= $xtra2.$pilecmd;
		push(@allDeps2, @{$dAR}); push(@checkF, @{$cAR});
	}


	

	
	my $postcmd ="";
	
	#from here on: merge XX vcf's into one
	if ($myParL ){
		#this string simply sorts all output files in correct numerical order.. doesn't touch file contents!
		my $sortedFileList = " | awk -F '.' '{print \$(NF-1),\$0}'  | sort -n -k1 | cut -f2 -d' '";
		#DEBUG
		#$postcmd .= "#DEBUG:\ncat `ls $tmpOut.*.lz4 $sortedFileList` > $ofasConsDir/Dbg.all.lz4\ncp $tmpOut.0.lz4 $ofasConsDir\n";
		$postcmd .= "mkdir -p $ofasConsDir;\n";
		$postcmd .= "if ls $qsubDirE/$smplNm.*.bed 1> /dev/null 2>&1 ;then echo \"Bed files still present, probably incorrect run\"; exit 33; else echo \"bed files deleted, looks good\"; fi\n\n";
		#old way to save file.. too much data for production environment
		my $saveVCF=1;
		if ($vcfFile eq ""){
			$saveVCF=0;
			$vcfFile = "$scrDir/$smplNm.fin.vcf.gz";
			$vcfFileS = "$scrDir/$smplNm.fin-sup.vcf.gz";
		}
		$vcfFile .= ".gz" unless ($vcfFile =~ m/\.gz$/);
		$postcmd .= "cat `ls $tmpOut.*.gz $sortedFileList` >$vcfFile ;\nrm -f $tmpOut.*.gz;\n";
		my $vcfSuff="";
		if ($SNPsuppStone ne "" ){
			$vcfFileS .= ".gz" unless ($vcfFileS =~ m/\.gz$/);
			$postcmd .= "cat `ls $tmpOut2.*.gz $sortedFileList` >$vcfFileS ;\nrm -f $tmpOut2.*.gz;\n";
			$postcmd .= "$bcftBin index $vcfFileS; $bcftBin index $vcfFile;\n";
			$vcfSuff = ".mrg.gz";
			$postcmd .= "$bcftBin concat -a --threads $samcores -O z -o $vcfFile$vcfSuff $vcfFile $vcfFileS;\n\n";				
		}
		$postcmd .= "zcat $vcfFile$vcfSuff | $vcfcnsScr $ofasCons.depStat $minDepth $minCallQual | $pigzBin -p $samcores -c >$ofasCons.gz ;\n\n"; #$refFA.fai
		$postcmd .= "rm -f $vcfFile$vcfSuff $vcfFileS.csi $vcfFile.csi;\n" if ($vcfSuff ne "");
		$postcmd .= "rm  -f $vcfFileS $vcfFile;\n" if (!$saveVCF );
		#} else {
		#	if ($SNPsuppStone ne "" ){die "support reads activated. combined SNP calling only works current with use of the \"-SNPsaveVCF 1\" MG-TK option. Aborting\n";}
			#$postcmd .= "#DEBUG\ncp $tmpOut.lz4 $ofasConsDir\n\n";
		#	$postcmd .= "zcat `ls $tmpOut.*.gz $sortedFileList`  |   $vcfcnsScr $ofasCons.depStat $minDepth  $minCallQual | $pigzBin -p $samcores -c >$ofasCons.gz ;\n\n"; #$refFA.fai 
		#}
		
		$postcmd .= "\necho \"Finished depthStat\"\n\n";
		$postcmd .= "rm -f $tmpOut*;\n";
		#$postcmd .= "$pigzBin -p $samcores $ofasCons;\n";

		$cmdAll .= "\n$postcmd\n" if ($run2ctg != 0);
	}
	if (exists($SNPIHR->{gffFile}) && !-e $SNPIHR->{genefna}){
		#requires python3 environment
		$cmdAll .= "\n$ctg2fas --gff $SNPIHR->{gffFile} --contig $ofasCons.gz --outFNA $SNPIHR->{genefna} --outFAA $SNPIHR->{genefaa};\n";
		$cmdAll .= "touch $SNPstone\n";
		$cmdAll .= "touch $SNPsuppStone\n" if ($SNPsuppStone ne "");
		#die $cmdAll."\n";
	}
	$cmdAll .= "\necho \"Finished contig to fasta\"\n\n";
	
	#die "$run2ctg\n$cmdAll\n";
	if ($myParL && !$runLocalTmp && $cmdAll ne ""){
		if ( ($overwrite || !-e "$ofasCons")){
			my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.cacSNP.sh",$postcmd,1,$memReq."G","Cons$x",join(";",@allDeps2),"",1,[],$QSBoptHR);
			$rdep =$dep;
		}
	}
	#die "$cmdAll\n $SNPIHR->{genefna}\n";
	if ($runLocalTmp && $cmdAll ne ""){#qsub all together now
		#this is the new way of doing this
		my $tmpS = $QSBoptHR->{tmpMinG};
		$QSBoptHR->{tmpMinG} = 70; #in GB
		my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.oSNPc.sh",$cmdAll,1,"5G","Cons$x",join(";",@allDeps2),"",1,[],$QSBoptHR);
		$rdep =$dep;
		$QSBoptHR->{tmpMinG} = $tmpS;
	}
	#$SNPIHR->{intermedVCF} = $oVcfCons;
	$SNPIHR->{cleanCmd} = $cleanCmd;
	#die;
	return ($SNPIHR,$rdep);
}


sub SNPconsensus_vcf2{
	die"not used after all..\n";
	#second version, that takes advantage of newer bcftools versions
	my ($SNPIHR)  = @_;
	my $frbBin = getProgPaths("freebayes");
	#my $kpathBin = getProgPaths("kpath");
	my $smtBin = getProgPaths("samtools");
	my $bcftBin = getProgPaths("bcftools");
	my $vcfcnsScr = getProgPaths("vcfCons_scr");
	#my $tabixBin = getProgPaths("tabix");
	#my $vcfLD = getProgPaths("vcfLib_dir");#"/g/bork3/home/hildebra/bin/vcflib/bin/";
	my $ctg2fas = getProgPaths("contig2fast_scr");
	my $pigzBin  = getProgPaths("pigz");
	#my $py3 = getProgPaths("py3activate",0);my $py3d = getProgPaths("pydeacti",0);

	#get parameteres
	my $samcores = 12;
	my %SNPinfo = %{$SNPIHR};
	my $QSBoptHR = $SNPinfo{QSHR};
	my $x = $SNPinfo{JNUM};
	my $jdep = ""; $jdep = $SNPinfo{jdeps} if (exists($SNPinfo{jdeps}));
	my $tmpdir = $SNPinfo{nodeTmpD};
	my $smplNm = $SNPinfo{smpl};
	my $refFA = $SNPinfo{assembly};
	my $qsubDirE = $SNPinfo{qsubDir};
	my $scrDir = $SNPinfo{scratch};
	my $bamcram = $SNPinfo{bamcram};
	my $splitFAsize = $SNPinfo{bpSplit};
	my $overwrite = $SNPinfo{overwrite};
	my $runLocalTmp = $SNPinfo{runLocal};
	my $maxSNPcores= $SNPinfo{maxCores};
	my $jobDepen= ""; $jobDepen=$SNPinfo{dependency} if (defined($SNPinfo{dependency}));
	#my $memReq = "20G";
	my $memReq = $SNPinfo{memReq};
	my $vcfFile = ""; $vcfFile = $SNPinfo{vcfFile} if (exists ($SNPinfo{vcfFile}));
	my $cmdFTag = $SNPinfo{cmdFileTag};
	my $firstInSample = 0;$firstInSample = $SNPinfo{firstInSample} if (exists($SNPinfo{firstInSample}));
	my $useFB = 1;
	$useFB = 0 if (uc($SNPinfo{SNPcaller}) eq "MPI");
	$vcfcnsScr = getProgPaths("vcfCons_FB_scr") if ($useFB);
	if ($runLocalTmp){
		$scrDir = $tmpdir;
		$samcores = $maxSNPcores;#$SNPinfo{split_jobs};
	}
	#basic caller options..
	my $minBQ=30; my $minMQ=30;
	

	#freebayes std options
	my $frAllOpts= "-u -i -m $minMQ -q $minBQ -C 1 -F 0.1 -k -X --pooled-continuous --report-monomorphic  --min-repeat-entropy 1 --use-best-n-alleles 2 -G 1 ";
	#bcftools options #-q = map qual -Q = base qual
	my $bcfAllOpts = "--count-orphans --min-BQ $minBQ -d 12000 --skip-indels --min-MQ $minMQ -a DP,AD,ADF,ADR,SP"; #
	if ($SNPinfo{SeqTech} eq "ONT"){
		$bcfAllOpts.=" -X ont ";
	} elsif ($SNPinfo{SeqTech} eq "PB"){
		$bcfAllOpts.=" -X pacbio-ccs ";
	} else{
		$bcfAllOpts.=" -X illumina ";
	}

	my $ofasCons = $SNPinfo{ofas};
	my @tar = ("");
	my $preTar = ${$SNPinfo{MAR}}[0];
	#die "$preTar\n";
	$ofasCons =~ m/(^.*)\/[^\/]+$/;
	my $ofasConsDir = $1."/";
	my $SNPstone = $ofasConsDir."SNP.cons.stone";
	system "rm -f $ofasCons*" if ($overwrite);
	my $rdep="";
	#prepare files..
	my $xtra = "";
	#prep ref FA
	$xtra .= "echo \"Preparing data\"\n";
	$xtra .= "mkdir -p $scrDir;\n";
	$xtra .= "$smtBin faidx $refFA;\n" unless (-e "$refFA.fai");
	$xtra .= "cp $refFA $refFA.fai $scrDir;\n";
	$refFA =~ m/\/([^\/]+$)/;
	$refFA = "$scrDir/$1";
	
	if (1){#new way: just use cram/bam that is available..
		if ($bamcram eq "cram"){ $tar[0] = "$scrDir/$smplNm.tmp.cram";} else {$tar[0] = "$scrDir/$smplNm.tmp.bam";}
		$xtra .= "ln -s $preTar $tar[0];\n";
	}elsif ($bamcram eq "cram"){
		my $tarFile1 = "$scrDir/$smplNm.tmp.cram";my $tarFile = "$scrDir/$smplNm.tmp.bam";
		$xtra .= "mkdir -p $scrDir\n";$xtra .= "$smtBin view -T $refFA -@ $samcores -b $preTar > $tarFile;\n";$tar[0] = $tarFile;
	} else {
		my $tarFile = "$scrDir/$smplNm.tmp.bam";$xtra .= "mkdir -p $scrDir;\ncp $preTar $tarFile;\n";$tar[0] = $tarFile;
	}
	my $indexBam = "$preTar.bai"; $indexBam =~ s/\.cram/\.bam/;
	if ( -e $indexBam){
		$xtra .= "cp $indexBam $tar[0].bai;\n"; 
	} else {
		$xtra .= "$smtBin index -@ $samcores -T $refFA  $tar[0];\n"; 
	}
		
	$xtra .= "echo \"Processing bams - mpileup\"\n";
	$SNPinfo{assembly} = $refFA;

	my $cmdAll = "";
	$cmdAll .= $xtra ;

	my $tmpOut = "$scrDir/$smplNm.cons.vcf";
	my $cmd2 = ""; 
	if ($useFB){
		$cmd2 .= "ulimit -s unlimited\n$frbBin -f $refFA  $frAllOpts  $tar[0] "; #--region '$curReg[$i]'
	} else {
		$cmd2 .= "$bcftBin mpileup --fasta-ref $refFA $bcfAllOpts --threads $samcores ";
		$cmd2 .= " $tar[0] | $bcftBin call --threads $samcores --ploidy 1 --multiallelic-caller -M --output-type v "; #this gets piped now.. #-o $tmpOut  \n"; 
	}
	
	my $sortedFileList = " | awk -F '.' '{print \$NF,\$0}'  | sort -n -k1 | cut -f2 -d' '";
	$cmd2 .= $sortedFileList . " | $vcfcnsScr $refFA.fai >$ofasCons 2> $ofasCons.depStat;\n\n ";
	
	$cmdAll .= $cmd2;
	$cmdAll .= "\necho \"Finished mpileup\"\n\n";
	
	$cmdAll .= "$pigzBin -p $samcores $ofasCons;\n";
	
	if (exists($SNPinfo{gffFile}) && !-e $SNPinfo{genefna}){
		#requires python3 environment
		$cmdAll .= "\n$ctg2fas --gff $SNPinfo{gffFile} --contig $ofasCons.gz --outFNA $SNPinfo{genefna} --outFAA $SNPinfo{genefaa};\n";
		$cmdAll .= "touch $SNPstone\n";
		#die $cmdAll."\n";
	}
	$cmdAll .= "\necho \"Finished contig to fasta\"\n\n";
	
	#die "\n$cmdAll\n";
	#die "$cmdAll\n $SNPinfo{genefna}\n";
	if ($runLocalTmp && $cmdAll ne ""){#qsub all together now
		my $tmpS = $QSBoptHR->{tmpMinG};
		$QSBoptHR->{tmpMinG} = 70; #in GB
		my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.oSNPc.sh",$cmdAll,$samcores,"5G","Cons$x",$jobDepen,"",1,[],$QSBoptHR);
		$QSBoptHR->{tmpMinG} = $tmpS;
	}
	#$SNPinfo{intermedVCF} = $oVcfCons;
	return (\%SNPinfo,$rdep);
}


#get fata consensus (scaffold) and extracts protein seqs from this..
sub SNPconsensus_fasta{
	my ($SNPIHR,$jdep) = @_;
	die "use SNPconsensus_vcf instead!\n";
	my %SNPinfo = %{$SNPIHR};
	my $qsubDirE = $SNPinfo{qsubDir};
	my $QSBoptHR = $SNPinfo{QSHR};
	my $x = $SNPinfo{JNUM};
	my $gff = $SNPinfo{gff};
	my $oVcfCons = $SNPinfo{intermedVCF};
	my $ofasCons = $SNPinfo{ofas};


	my $vcfcnsScr = getProgPaths("vcfCons_scr");
	my $postcmd = "zcat $oVcfCons.gz | $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat;\n\n";
	
	#if (system "$fix\n"){system "rm $oVcfCons*"; print "RM\n";}
	my ($dep,$qcmd) = qsubSystem($qsubDirE."CallFasta.sh",$postcmd,1,"40G","Cons$x",$jdep,"",1,[],$QSBoptHR);
}























