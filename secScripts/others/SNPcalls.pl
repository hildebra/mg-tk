#!/usr/bin/env perl
#builds from a contig file several stats to seperate contigs into single species

#deprecated, no longer used in overall MATAFILER run
#still contains a lot of cool tricks that could be later re-integrated with MATAFILER..

#./SNPcalls.pl /g/bork3/home/hildebra/data/AnnaPry/161013_M00758_0551_000000000-ATWLT/Bsubtillis168.fasta /g/scb/bork/hildebra/Tamoc/ANNA/GlbMap/BS168 BS168
#./SNPcalls.pl -nJobs 400 -reportAllSites 1 -refGenome /g/bork3/home/hildebra/data/SNP/Drama1/GlbMap/BacVulgatus/BacVulgatus.fasta -inD /g/bork3/home/hildebra/data/SNP/Drama1/GlbMap/BacVulgatus -name BacVulgatus -map /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/drama_HMP.map -outDir /g/bork3/home/hildebra/data/SNP/Drama1/GlbMap/BacVulgatus/SNP/
#./SNPcalls.pl -nJobs 400 -reportAllSites 1 -refGenome /g/bork3/home/hildebra/data/SNP/Drama1/GlbMap/PrevCopri/PrevCopri.fasta -inD /g/bork3/home/hildebra/data/SNP/Drama1/GlbMap/PrevCopri/ -name PrevCopri -map /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/drama_HMP.map -outDir /g/bork3/home/hildebra/data/SNP/Drama1/GlbMap/PrevCopri/SNP/
use warnings;
use strict;
use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw(systemW readMap readFasta  median);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Statistics::Descriptive;
use threads ('yield',
                 'stack_size' => 64*4096,
                 'exit' => 'threads_only',
                 'stringify');
sub countRdHits;
sub freebayes_multi;
sub lowfreqSNP;
sub freebayes_genotype;
sub correctBamHD;
sub stat_vcf;
sub sinvictCall;
sub createConsensus;
sub pilonCons;
sub help;



#SNP discovery stage
#Freebayes -> VCF (does indel realignments) alt: Varscan
#Genotyping
#Platypus (assumes stable SNP freq), free bayes (assumes stable SNP freq)
#--debug
#binaries
my $frDir = "/g/bork3/home/hildebra/bin/freebayes/bin/";
my $frbBin = "$frDir/freebayes";
my $vcfTools = "/g/bork3/home/hildebra/bin/vcflib/bin/";
my $vtBin = "/g/bork3/home/hildebra/bin/vt/vt";
my $smtBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $consCntupVCFscr = getProgPaths("consCntupVCF_scr");

my $varscanBin = "java -jar /g/bork3/home/hildebra/bin/varscan2/VarScan.v2.3.8.jar ";
my $vcfFilBin = "$vcfTools/vcffilter";
my $bcBin = "/g/bork3/home/hildebra/bin/bcftools-1.3.1/./bcftools"; 
my $vcftls = "/g/bork3/home/hildebra/bin/vcftools/vcft/bin/./vcf-merge";
#my $lfBin = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/SebSNP/tools/./lofreq";
my $lfBin = "/g/bork3/home/hildebra/bin/lofreq_star-2.1.2/bin/./lofreq";
my $concatScr = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/SNP/./concatVCF.pl";
my $vcfcnsScr = "perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/SNP/vcf2cons.pl ";
my $comDepWinScr = "perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/SNP/comDepWins.pl";
my $colStatsScr = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/SNP/collectStatsCons.pl";

my $frAllOpts= "-u -i -C 1 -F 0.1 -k -X --pooled-continuous --report-monomorphic  --min-repeat-entropy 1 --use-best-n-alleles 2 -G 1 ";


my $seqtkBin = getProgPaths("seqtk");
my $bam2cnsBin = getProgPaths("bam2cns");
my $sinvictBin = getProgPaths("sinvict");
my $br2dCntBin = getProgPaths("br2dCnt");

my $QSBoptHR = emptyQsubOpt(1,"");



#CONFIG------------------
my $idx = 0;
#my $LOC = "ali"; #ALL, ali, DK, US, HD, SWE
my @LOCs = ("ALL");#,"ALL", "ali", "US", "HD");
my @MAPs = (20); my @BQUALs = (25);  #lowfreq
#choice of calling algo
my $DoFreeBayes=1;
my $DoConvert2fasta=1; my $addConsSNPs = 1;

my $getStats=0; #get statistics on callable positions and average coverage / bam
#my @MAPs = (0); my @BQUALs = (0); #freebayes
#@BQUALs = (6,10,20,25,30,32,33,34); @MAPs = (0,5,10,15,20,25,30); @LOCs = ("ali");
my $redo=0; 
my $overwrite = 0; my $redoFreeBayes = 0; #controls freebays and perl script
#CONFIG------------------


#hard coded paths for internal testing
my @refFAs = ("/g/bork5/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa","/g/bork5/hildebra/results/TEC2/v5/T3/T3.mini2.3smpl.fna",
"/g/bork5/hildebra/results/TEC2/v5/T4/T4.mini2.3smpl.fna",
#"/g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T5/R_filt/contigs/MM3.ctgs.fna",
"/g/bork5/hildebra/results/TEC2/v5/T6/TEC6.ctgs.rn.fna");
my @name = ("T2d","T3d","T4d","T6d");
my $inDir = "";#/g/scb/bork/hildebra/SNP/GNMass3/GlbMap/$name[$idx]/";
my $mapF = "";#g/bork5/hildebra/data/metaGgutEMBL/MM_at_v5_T2subset.txt";
my $odir = "";#/g/bork3/home/hildebra/data/TEC2_related/v5/SNPcalls/$name[$idx]/";
my $tmpdir = "/scratch/bork/hildebra/SNP/";

#internal flow control
my $doPar = 0; #really stupid script, use myPar instead
my $myPar=1; my $doQsub = 1;
my $ncore = 24;my $splitFAsize = 50000;
my $reportAllSites=0; my $jobs2split2=-1;
my @thrs; my @allDeps;

if (1){#report all with freebayes.. extremely expensive job
	$splitFAsize = 20000;$reportAllSites=1;
}



if (@ARGV > 0){#non-hardcoded paths
	my $refFApre = ""; my $nmPre = "";my $mapFpre="";my$inDpre="";
	GetOptions(
		"help|?" => \&help,
		"map=s" => \$mapFpre,
		"inD=s" => \$inDpre,
		"name=s" => \$nmPre,
		"refGenome=s" => \$refFApre,
		"outDir=s" => \$odir,
		"reportAllSites=i" => \$reportAllSites,
		"nJobs=i" => \$jobs2split2,
		"ntPerRun=i" => \$splitFAsize,
	);
	if ($odir eq ""){$odir = $inDpre."/SNP/";}
	if ($nmPre eq "" || $refFApre eq "" || $inDpre eq ""){help(); die;} #$mapFpre eq "" || 
	$mapF = $mapFpre;$inDir = $inDpre;
	$idx=0;@LOCs=("ALL");
	@refFAs = split/,/,$refFApre; @name=split /,/,$nmPre; 
	#$inDir =~ m/^(.*\/)GlbMap\/.*/;
	#$mapF = $1 . "/LOGandSUB/inmap.txt";
	die "Can;t find mapping file at $mapF!\n" if (!-e $mapF);
} else { #old tec2 calls
	$inDir = "/g/scb/bork/hildebra/SNP/GNMass3/GlbMap/$name[$idx]/";
	$mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v5_T2subset.txt";
	$odir = "/g/bork3/home/hildebra/data/TEC2_related/v5/SNPcalls/$name[$idx]/";

}
die "doQsub is invalid without myPar\n" if ($doQsub && !$myPar);

$inDir .= "/" unless ($inDir =~ m/\/$/);
my $refFA = $refFAs[$idx];

#file location checks
die "Can't find indir $inDir\n" unless (-d $inDir);
die "Can't find refgenome $refFA\n" unless (-e $refFA);
system "$comDepWinScr $inDir\n" unless (-e "$inDir/$name[$idx].all.coverage.gz.window");

my $ldir = $odir."/LOGs/";
system "mkdir -p $odir" unless (-d $odir);system "mkdir -p $ldir" unless (-d $ldir);
my $qsubDir = "$odir/qsubs/";
system "mkdir -p $qsubDir" unless (-d $qsubDir);
system "$smtBin faidx $refFA\n" unless (-e "$refFA.fai");

#remove old cons dir (moved to new location)
#system "rm -r $inDir/consFasta/";

my $conFastaDir = "$odir/consFasta_0.5/"; #_0.5
system "mkdir -p $conFastaDir" unless (-d $conFastaDir);
system "mkdir -p $tmpdir" unless (-d $tmpdir);


#my $mapF = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v4_T2subset.txt";
#my $refDir = "/g/bork1/hildebra/SNP/GNMassSimu//simulated_metaG2SI_3/";
my $refDir = $ARGV[0];
my $global=0;
my %map; my %AsGrps;
my ($hr,$hr2) = readMap($mapF,0,\%map,\%AsGrps,1); %map = %{$hr};
my @samples = @{$map{opt}{smpl_order}}; my @bams; my @bamSiz; my @bamRdCnts;
my @consFastas; my @vcfConsDep; my @vcfCons;
#die "@samples\n";

print "Looking at indir $inDir\nlocation @LOCs\nquals @BQUALs and MAQ @MAPs\n";
print "Writing to $odir\n";

my $bamListFiles = "$odir/".$LOCs[0]."_$name[$idx]_inbams.txt";
#my $tmps = `grep -c '>' $refFA`; $tmps=~m/^(\d+)/;my $tarCtgNum=$1;

$hr = readFasta($refFA); my @tarCtgs = keys(%{$hr});  my $tarCtgNum = scalar(@tarCtgs); 
my $tarLength=0; foreach my $ctg (@tarCtgs){$tarLength+= length(${$hr}{$ctg});}


if ($jobs2split2 > 1){
	$splitFAsize = int(($tarLength-1000)/$jobs2split2);
}

my $statsStr = ""; #get some stats on run
foreach my $LOC (@LOCs){

	my $smpCnt=0;
	$bamListFiles = "$odir/${LOC}_inbams.txt";
	open O,">$bamListFiles";  my %alreadyIncl; open L1,">$ldir/${LOC}_mapdRds.txt";
	my $atB=0;
	foreach my $cs (@samples){
		
		#subsets of bams
		if ($LOC eq "ali"){
			#only alien
			#print "$cs $map{$cs}{wrdir}";
			next unless ($map{$cs}{wrdir} =~ m/alien-/);
		} elsif ($LOC eq "HD"){
			#only HD
			next unless ($map{$cs}{wrdir} =~ m/alien/ || $map{$cs}{wrdir} =~ m/donald/ || $map{$cs}{wrdir} =~ m/halbarad/ || $map{$cs}{wrdir} =~ m/pluto/|| $map{$cs}{wrdir} =~ m/thistle/);
		} elsif ($LOC eq "SWE"){
			#only US
			next unless($map{$cs}{SmplID} =~ m/^SWE/);
		} elsif ($LOC eq "DK"){
			#only DK
			die;
		} elsif ($LOC eq "US"){
			#only US
			next unless($map{$cs}{SmplID} =~ m/^HM/);
		} elsif ($LOC eq "ALL") {;}
		
		$atB++;
		next if ($atB == 1);
		#print "$LOC $cs $map{$cs}{wrdir}\n";
		#debug
		#next if ($map{$cs}{mapFinSmpl} =~ m/MM3/ || $map{$cs}{mapFinSmpl} =~ m/MM4/);
		my $tar = $inDir."/$name[$idx]_".$map{$cs}{mapFinSmpl}."-0-smd.bam";
		print "$name[$idx]_".$map{$cs}{mapFinSmpl}."\n";
		unless (-f $tar){print "file $tar does not exist\n";  next;}
		next if (exists $alreadyIncl{$tar});
		print $tar unless (-e $tar);
		$alreadyIncl{$tar} = 1;
		my $mbSiz = (-s $tar) / (1024 * 1024);
		my $ofasCons = "$conFastaDir/$name[$idx]_".$map{$cs}{mapFinSmpl}.".cons.fna";
		my $oVcfCons = "$conFastaDir/$name[$idx]_".$map{$cs}{mapFinSmpl}.".pre.vcf";
		push (@consFastas,$ofasCons);
		my $totHits = countRdHits($tar,30);
		if ($getStats){
			my ($callable,$consCnt) = getConsStats($ofasCons);
			my $cov=0;my$bp;
			if (-e "$tar.coverage.gz.percontig"){
				my $ctmp = `cat $tar.coverage.gz.percontig`;
				foreach (split /\n/,$ctmp){
					die "pattern not found: $_\n" unless(m/_L=(\d+)=\t([\d\.]+)$/);
					$bp += $1; $cov += $1*$2;
				}
				$cov /=$bp;
			}
			$statsStr .= "$cs\t$totHits\t$callable\t$consCnt\t$cov\n";
		}

		#next if ($mbSiz < 0.7);
		#print $map{$cs}{mapFinSmpl}."\t".$totHits."\n";
		print L1 $map{$cs}{mapFinSmpl}."\t" .$mbSiz."\t".$totHits."\n";
		if (0){
			print $map{$cs}{mapFinSmpl}."\t" .$mbSiz."\t".$totHits."\n";
			next;
		}
		push(@vcfCons,$oVcfCons);
		if ($DoConvert2fasta){ #creates a genome consensus sequence
			#if ($totHits < 1e6){$myPar=0;
			#} else {$myPar=1;}
			my $splitFAsizeLoc=-1;#$splitFAsize;
			if ($totHits > 1e6){#dynamic attempt
				$splitFAsizeLoc = int(($tarLength-1000)/int($totHits / 2e5));
			}
			print "$splitFAsizeLoc XXX\n";
			#next unless ($tar =~ m/MH3-/);
			push(@vcfConsDep,createConsensus($tar,$oVcfCons,$ofasCons,$atB,$overwrite,$tarCtgNum,$splitFAsizeLoc));
			#die;
			next;
		}
		if ($totHits < 1000){print "skip $map{$cs}{mapFinSmpl}  $totHits\n";next;}
		print O $tar."\n";
		#$tar = correctBamHD($tar,$contigMarker);
		system "$smtBin index $tar" unless (-e "$tar.bai");
		push(@bams,$tar);
		push(@bamSiz,$mbSiz);
		$smpCnt++;
		push (@bamRdCnts,$totHits );
		print $map{$cs}{mapFinSmpl}."\t" .$mbSiz."\t".$totHits."\n";
	}
	close O; close L1;
	if ($getStats){
		print "writing stats\n";
		open O,">$odir/runStats.txt";
		print O $statsStr;
		close O;
	}
	if ($DoConvert2fasta){
		my $outVCF = "$conFastaDir/$name[$idx].$LOC.mrg.vcf";
		#die "@vcfCons\n";
		print "Merge VCF job..\n";
		merge_vcf($conFastaDir,\@vcfConsDep,$refFA,$outVCF);
		
	}


	#die "not needed\n";
	if ($smpCnt ==0 ){print "No bams given\n" ; next;}

	#die;
	#/g/bork3/home/hildebra/bin/bamtools-master/bin/bamtools coverage -in aln.bam | coverage_to_regions.py ref.fa 500 >ref.fa.500.regions

	#old way: just call with freebayes everything
	my $custID=0;
	foreach my $MAP (@MAPs){ #loops over MAPs and BQs to benchmark algo's
	#last;
		foreach my $BQUAL (@BQUALs){
			
			my $odir2 = $odir."${LOC}_lowf2.1_$MAP"."_$BQUAL/";
			if ($redo){system "rm -r $odir2";}
			
			
			if ($DoFreeBayes){
				#$odir2 = $odir."FB_${LOC}_${MAP}_$BQUAL/";
				freebayes_multi($bamListFiles,$odir,$LOC,$MAP,$BQUAL ,$reportAllSites);
			}
			if (0){#sinvict
				sinvictCall(\@bams,$odir,$LOC,$MAP,$BQUAL,$name[$idx],$custID);
				die;
			}
			if ($DoConvert2fasta){#consensus
				my $mergF = "$odir/$name[$idx].cons.vcf";
				system "$consCntupVCFscr $refFA ".join (" ",@consFastas) . " > $mergF";
				system "cat $mergF | /g/bork3/home/hildebra/bin/vcflib/bin/./vcfstreamsort -a | /g/bork3/home/hildebra/bin/vcflib/bin/./vcfuniq > $mergF.1\n";
				system "rm $mergF;mv $mergF.1 $mergF";
				system "\nbgzip -f $mergF; tabix -f -p vcf $mergF.gz;\n";
				$mergF.=".gz";
				freebayes_genotype($odir2,$mergF,$refFA,$name[$idx],"",$MAP,$BQUAL,$LOC,"cons");
			} 
			if (0){#lowfreq
				$qsubDir = "$odir2/qsubs/";
				my $ltag = "lowfr";
				if ($addConsSNPs){
					#die "$consCntupVCFscr $refFA ".join " ",@consFastas;
					#$ltag .="_cons";
				}
				my ($mrgF,$dep) = lowfreqSNP($odir2,\@bams,$refFA,$name[$idx],$MAP,$BQUAL,$custID);
				freebayes_genotype($odir2,$mrgF,$refFA,$name[$idx],$dep,$MAP,$BQUAL,$LOC,$ltag);
			} 
			$custID++;
			#last;
		}
	}
}
my $colCmd .= "$colStatsScr $conFastaDir\n";

print "collect stats:\n$colCmd\n";
print "All SNPs Done!\n$odir\n"; 


exit(0);



die;
my $LOC = $LOCs[0]; $bamListFiles = "$odir/".$LOC."_inbams.txt";



#-------------------------------------------------------------
#stats on num alleles:
#-------------------------------------------------------------

my %varmat; my %avgMat; my %quantMap;
#print "\n\n".join("\t",(sort {$a <=> $b} @BQUALs));
foreach my $MAP (sort {$a <=> $b} @MAPs){
	#print "\n$MAP";
	foreach my $BQUAL (sort {$a <=> $b} @BQUALs){
		
		my $odir2 = $odir."${LOC}_lowf2.1_$MAP"."_$BQUAL/";
		my ($numSNPs,$quntQual,$avgQual) = stat_vcf("$odir2/FB_reGT_${LOC}_$MAP"."_$BQUAL.vcf");
		chomp ($numSNPs);
		$varmat{$MAP}{$BQUAL} = $numSNPs;
		$avgMat{$MAP}{$BQUAL} = $avgQual;
		$quantMap{$MAP}{$BQUAL} = $quntQual;
		#print "\t$numSNPs";
	}
}

print "\nNumHits\n".join("\t",(sort {$a <=> $b} @BQUALs));
foreach my $MAP (sort {$a <=> $b} @MAPs){
	print "\n$MAP";
	foreach my $BQUAL (sort {$a <=> $b} @BQUALs){
		print "\t".$varmat{$MAP}{$BQUAL}
	}
}
print "\n";

print "\nAvgQual\n".join("\t",(sort {$a <=> $b} @BQUALs));
foreach my $MAP (sort {$a <=> $b} @MAPs){
	print "\n$MAP";
	foreach my $BQUAL (sort {$a <=> $b} @BQUALs){
		print "\t".$avgMat{$MAP}{$BQUAL}
	}
}
print "\n";
print "\n95th percentile qual\n".join("\t",(sort {$a <=> $b} @BQUALs));
foreach my $MAP (sort {$a <=> $b} @MAPs){
	print "\n$MAP";
	foreach my $BQUAL (sort {$a <=> $b} @BQUALs){
		print "\t".$quantMap{$MAP}{$BQUAL}
	}
}
print "\n";
print "\n\nAll done analysis\n";exit(0);




sub merge_vcf{
	my ($vcfDir,$arDeps,$ref,$finalVCF) = @_;
	#die "@$arVCFs\n";
	my $mrgCmd = "$bcBin merge $vcfDir/*.pre.vcf.gz > $finalVCF.1 \n";
	$mrgCmd .= "$vtBin normalize -r $ref $finalVCF.1 >  $finalVCF";
	$mrgCmd .= "\nbgzip $finalVCF; tabix -p vcf $finalVCF.gz;\n";
	$mrgCmd .= "rm $finalVCF.1\n";
	#$mrgCmd .= "rm $vcfDir/*.pre.vcf.gz";

	#die "$mrgCmd\n";
	my ($dep,$qcmd) = qsubSystem($qsubDir."mergeVCF.sh",$mrgCmd,1,"40G","ConsMrg",join(";",@$arDeps),"",1,[],$QSBoptHR);
	
}
sub stat_vcf(){
	my ($inVcf) = @_;
	#my $numSNPs = `zcat $inVcf | grep -v -c '^#'`;
	my $numSNPs = `grep -v  '^#' $inVcf | grep -c -v 'MM4TEC2__C1_L=46754='`;
	my $mpt = `grep -v '^#' $inVcf | grep -v 'MM4TEC2__C1_L=46754=' | cut -f6`;
	my @quals = sort {$a <=> $b} split /\n/,$mpt;
#get AO and RO	
	my @AO; my @RO; my @AORO;
	$mpt = `grep -v '^#' $inVcf | grep -v 'MM4TEC2__C1_L=46754=' | cut -f8`;
	my @tmp =  split /\n/,$mpt;
	foreach (@tmp){m/;RO=(\d+)/;my $RR=$1;push(@RO,$RR);m/;AO=(\d+)/;my $AR=$1;push(@AO,$AR);push(@AORO,$AR/($AR+$RR));}
	@AORO = sort {$a <=> $b} @AORO;
	#die "@AORO\n";
	#die "\n$mpt\n";
	my $qunQuals = $AORO[int ($#AORO*0.10)];
	my $avgQual = $AORO[int ($#AORO*0.5)];
	return ($numSNPs,$qunQuals,$avgQual);
}

sub downSmpl(){
	my ($fAR,$tar) = @_;
}

sub getRegionsBam($){
	my ($splitFAsizeL) = @_;
	my @curReg;
	my $regionFile = "$odir/regions_par_$splitFAsizeL.txt";
	system "$smtBin faidx $refFA" unless (-e "$refFA.fai");
	if (!-e $regionFile && system "python $frDir/fasta_generate_regions.py $refFA.fai $splitFAsizeL > $regionFile\n"){
		print "python $frDir/fasta_generate_regions.py $refFA.fai $splitFAsizeL > $regionFile\n";
		die "Can't gemerate region file:\n";
	}
	open my $handle, '<', $regionFile;
	chomp(@curReg = <$handle>);
	close $handle;
	return (\@curReg);
}


sub getConsStats{
	my ($cnonsFas) = @_;
	return ("?,?") unless (-e $cnonsFas);
	my $str=`grep '>' $cnonsFas`;
	my @spl = split /\n/,$str;
	my $cov=0;my $repl=0;
	foreach (@spl){
		m/COV=(\d+) REPL=(\d+)/;$cov+=$1;$repl+=$2;
	}
	#die "$repl\n";
	return ($cov,$repl);
}
sub lowfreqSNP(){
	my ($odir,$fileAR,$ref,$name,$MAP,$BQ,$ID) = @_;
	my $mergF = "$odir/$name.merge.vcf";
	if (-e "$mergF.gz" ){return $mergF.".gz","";}
	system "rm -rf $odir;mkdir -p $odir";
	my @files = @{$fileAR};
	my $callDescr = "lofreq.maq$MAP.bq$MAP";
	my $refAR = getRegionsBam($splitFAsize);
	my @curReg = @{$refAR};
	my $lfOPT = "-m $MAP -Q $BQ "; #lf call .. $lfOPT #-J 50
	#die "@curReg\n\n";
	system "mkdir -p $qsubDir" unless (-d $qsubDir);
	#  -r | --region STR            Limit calls to this region (chrom:start-end) [null]
	#my $vcf = $outDir. "lofreq/$name.". ".lofreq.maq20.vcf";
	my $cnt=0; my @vcfs; my @depends;
	foreach my $f (@files){
		#last;
		$f =~ m/^.+\/([^\/]+)-/; my $tNam = $1;
		my $numThr = 1;my $useMem = "10G"; my $doPar = 0;
		if ($bamSiz[$cnt]> 100){$doPar=1;}#$numThr=20;$useMem="3G";}
		my $outF = "$odir/$name.$tNam.$callDescr.vcf";
		
		push(@vcfs,"$outF.gz");
		next if (-e "$outF.gz");
		my $tmpOut = "$odir/$tNam.first_pass";
		#my $cmd = "$lfBin call -m 20 -f $ref -o $outF $f\n";
#new 2.1 version
		#--no-default-filter -r "scaffold_2:962678-1925355"
		#my $cmd = "$lfBin call-parallel --pp-threads $numThr -q 30 -m 20 -f $ref -o $outF $f\n";
		my $cmd = "";
		my $postcmd = "";
		my $depX = "";
		if ($doPar){
			my @ddps; my $dcnt=0;
			foreach(@curReg){
				$cmd="$lfBin call -f $ref -r $_ $lfOPT -o $tmpOut.$dcnt.vcf $f\n"; #-q 20 -Q 20
				
				#die $cmd;
				my ($dep,$qcmd) = qsubSystem($qsubDir."lofreq$cnt.$dcnt.sh",$cmd,1,"10G","${ID}LF$cnt.$dcnt","","",1,[],$QSBoptHR);push(@ddps,$dep);
				$dcnt++;
			}
			$postcmd = "cat $tmpOut.* | /g/bork3/x86_64/bin/python /g/bork3/home/hildebra/bin/vcflib/bin/vcffirstheader | /g/bork3/home/hildebra/bin/vcflib/bin/./vcfstreamsort -a | /g/bork3/home/hildebra/bin/vcflib/bin/./vcfuniq > $outF\n";
			$postcmd .= "rm $tmpOut.*\n";
			$postcmd .= "bgzip $outF; tabix -p vcf $outF.gz\n";
			my ($dep,$qcmd) = qsubSystem($qsubDir."lofreq$cnt.post.sh",$postcmd,1,"10G","${ID}LF$cnt",join(";",@ddps),"",1,[],$QSBoptHR);
			$depX = $dep;
		} else { #single core call for small bams..
			$cmd="$lfBin call $lfOPT -f $ref -o $outF $f\n";
			$postcmd = "bgzip $outF; tabix -p vcf $outF.gz\n";
			my ($dep,$qcmd) = qsubSystem($qsubDir."lofreq$cnt.sh",$cmd."\n".$postcmd,1,"10G","${ID}LF$cnt","","",1,[],$QSBoptHR);
			$depX = $dep;
		}
		$cmd .= "$f\n";
		#"$bcBin convert -Ou $tmpOut | $bcBin norm -r $refFA -m - -O b -o $outF -\n";
#die $postcmd."\n";
#		my $postcmd = " $vcfTools/./vcfbreakmulti $tmpOut \\
		#| $vtBin normalize -r $refFA - > $outF\nrm $tmpOut\n";
		#die "$cmd\n";
		#next;
		$cnt++;
		push(@depends,$depX) unless ($depX eq "");
	}
	my $mrgCmd = "$bcBin merge $odir/$name.*.$callDescr.vcf.gz > $mergF.1\n";
	$mrgCmd .= "$vtBin normalize -r $ref $mergF.1 >  $mergF";
	$mrgCmd .= "\nbgzip $mergF; tabix -p vcf $mergF.gz;\n";
	#my $mrgCmd = "$concatScr $odir/$name.*.$callDescr.vcf > $name.merge.vcf";
	#die $mrgCmd."\n";
	my ($dep,$qcmd) = qsubSystem($qsubDir."mergeVCF.sh",$mrgCmd,1,"40G","LFmr$ID",join(";",@depends),"",1,[],$QSBoptHR);
	#print $qcmd."\n";
	return ($mergF.".gz",$dep);
}


sub freebayes_genotype(){#recall variants at specified positions..
	my ($odir,$mergeVCF,$ref,$name,$dep,$MAP,$BQ,$LC,$tag) = @_;
	#my @files = @{$fileAR};
	system "mkdir -p $odir" unless (-d $odir);
	my $fbOpt = "--ploidy 1 --pooled-discrete --pooled-continuous --report-genotype-likelihood-max --hwe-priors-off --min-repeat-entropy 1 --use-best-n-alleles 4";
	my $cmd = "$frbBin -f $ref -@ $mergeVCF -l $fbOpt -p 1 -m $MAP -q $BQ -L $bamListFiles > $odir/FB_reGT_${tag}_${LC}_$MAP"."_$BQ.vcf\n";
	#die "$cmd\n";
	my ($depX,$qcmd) = qsubSystem($qsubDir."FB_genotype.sh",$cmd,1,"40G","FBGT$tag",$dep,"",1,[],$QSBoptHR);
}


sub sinvictCall(){
	my ($bamLst,$odir,$LOC,$MAQ,$BQ,$name,$ID) = @_;
	my @bams = @{$bamLst};
	my $odirX = $odir."SV_${LOC}_${MAQ}_$BQ/";
	$odir = $odirX."/RC/";
	system "mkdir -p $odir" unless (-d $odir);
	my $callDescr = "sinvict.maq$MAQ.bq$BQ";
	my $cnt=0; my @vcfs; my @depends;
	foreach my $f (@bams){
		#last;
		$f =~ m/^.+\/([^\/]+)-/; my $tNam = $1;
		my $numThr = 1;my $useMem = "10G"; my $doPar = 0;
		if ($bamSiz[$cnt]> 100){$doPar=1;}#$numThr=20;$useMem="3G";}
		my $outF = "$odir/$name.$tNam.$callDescr.vcf";
		my $outFRC = "$odir/$name.$tNam.$callDescr.readcount";
		next if (-e $outFRC);
		my $cmd = "set -e\n";
		$cmd  .= "$br2dCntBin -q $MAQ -b $BQ -w 10 -f $refFA $f > $outFRC\n"; #--site-list  file containing a list of regions to report readcounts within.
		#die $cmd."\n";
		push(@vcfs,"$outF.gz");
		next if (-e "$outF.gz");
		my $tmpOut = "$odir/$tNam.first_pass";
		my ($dep,$qcmd) = qsubSystem($qsubDir."bamRDcnt$cnt.sh",$cmd."\n",1,"10G","${ID}BR$cnt","","",1,[],$QSBoptHR);
		push(@depends,$dep) unless ($dep eq "");
		$cnt++;
	}
	#my $odirX = $odir."SiVt/"; system "mkdir -p $odirX" unless (-d $odirX);
	#RCs need to be on clean dir with no subdirs etc
	my $cmd = "$sinvictBin --min-depth 3 --qscore-cutoff 90 -t $odir -o $odirX\n";
	my ($dep,$qcmd) = qsubSystem($qsubDir."sinvct.sh",$cmd."\n",1,"10G","${ID}ST",join(";",@depends),"",1,[],$QSBoptHR);
}

sub createConsensus($ $ $ $ $ $ $){
	my ($tar,$oVcfCons,$ofasCons,$x,$overwrite,$tarCtgNum,$splitFAsize)  = @_;
	my $cmd = ""; my $rdep="";
	my $hereCtgs = 0;
	my $short=0; #if 1, this is a fast job, no qsub..
	if (-e $oVcfCons){
		#only works for fastas
		#my $tmps = `grep -c '>' $oVcfCons`; $tmps=~m/^(\d+)/;$hereCtgs=$1;
	}
	#print "$tarCtgNum != $hereCtgs\n";
	#die;
	#if ($tarCtgNum != $hereCtgs){$overwrite=1;}
	if (!-e $oVcfCons && !-e "$oVcfCons.gz"){$overwrite=1;}

	#	my $cmd = "$smtBin bam2fq $tar | $seqtkBin seq -A >$ofasCons";#only gives reads back
#	$cmd = "$bam2cnsBin --bam $tar  --coverage 10 --qual-weighted --qv-offset 28 --prefix $ofasCons\n";#--ignore-weak-reads=20 
#	$cmd = "$smtBin mpileup -B  -C 20 -d 5000 -Q 28 -v -u -f $refFA  $tar | $vcfcnsScr  >$ofasCons\n";
	#bcftools consensus -f /g/bork3/home/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa test
	$cmd = "$frbBin -f $refFA  $frAllOpts -m 30 -q 30 ";
	my @allDeps2;
	#die $cmd;
#	system $cmd;
	#implement in parallel as too slow in single core mode :/
	my $tmpOut = "$tmpdir/${name[$idx]}.$x.cons.vcf";
	my @curReg = ("1");
	my $myParL=0;
	if ($splitFAsize>0){$myParL=1;}
	if ($myParL){ #no, don't redo freebayes part
		my $refAR = getRegionsBam($splitFAsize);
		@curReg = @{$refAR}; 
	}
	if (!$redoFreeBayes && -e "$oVcfCons.gz"){
		@curReg = ("1");$myParL=0;
	}
	my $qsubDir = $odir."/qsubs/";
	system "mkdir -p $qsubDir" unless (-d $qsubDir);
	for (my $i=0;$i<@curReg;$i++){
		my $cmd2 = $cmd ;
		if ($myParL){
			$cmd2 .= " --region '$curReg[$i]' $tar > $tmpOut.$i \n";
		} else {
			$cmd2 .= " $tar > $oVcfCons\n";
			$cmd2 .= "bgzip $oVcfCons; tabix -p vcf $oVcfCons.gz\n";
			if (!$redoFreeBayes && -e "$oVcfCons.gz"){$cmd2 = "" ; $short =1;}
			$cmd2 .= "zcat $oVcfCons.gz | $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat\n\n";
			#| $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat\n";
		}
		#if (-s $tmpOut){[print " tmp exists ";next;}
		#die $cmd2;
	
		if ($short && $overwrite){
			#die "$cmd2\n";
			systemW $cmd2;
		}elsif ($overwrite ){
			#die "overwr\n";
			my ($dep,$qcmd) = qsubSystem($qsubDir."FB_Cons$x.$i.sh",$cmd2,1,"5G","FBC$x.$i","","",1,[],$QSBoptHR);
			push (@allDeps2,$dep);
			$rdep =$dep;
		}
	}
	
	my $postcmd ="";

	if ($myParL && $overwrite ){
		my $vcf1stHd = "/g/bork3/x86_64/bin/python /g/bork3/home/hildebra/bin/vcflib/bin/vcffirstheader";
		my $vcfStrSrt = "/g/bork3/home/hildebra/bin/vcflib/bin/./vcfstreamsort";
		$postcmd .= "cat $tmpOut.* | $vcf1stHd | $vcfStrSrt -a > $oVcfCons\n";
		$postcmd .= "bgzip $oVcfCons; tabix -p vcf $oVcfCons.gz\n";
		$postcmd .= "zcat $oVcfCons.gz | $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat\n\n";
		#| $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat\n\n";
		$postcmd .= "rm $tmpOut.*\n";
		#die "$postcmd\n";
		my ($dep,$qcmd) = qsubSystem($qsubDir."postCmd.sh",$postcmd,1,"40G","Cons$x",join(";",@allDeps2),"",1,[],$QSBoptHR);
		$rdep =$dep;
		#die;
	}
	if(-e $oVcfCons && !-e "$oVcfCons.gz"){
		my $fix .= "bgzip $oVcfCons; tabix -p vcf $oVcfCons.gz\n";
		system "$fix\n";
	}
	#die "$ofasCons\n$oVcfCons\n";
	if(  !$overwrite && -e "$oVcfCons.gz" && (!-e $ofasCons || -z $ofasCons)){
		my $fix .= "zcat $oVcfCons.gz | $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat\n\n";
		#die "$fix\n";
		if (system "$fix\n"){system "rm $oVcfCons*"; print "RM\n";}
	}
	#die;
	return $rdep;
}

sub freebayes_multi(){ #freebayes in multicore mode
	my ($bamLst,$odir,$LOC,$MAQ,$BQ, $reportAll) = @_;
	$odir = $odir."FB_p1_p2_l_${LOC}_${MAQ}_$BQ/";
	system "mkdir -p $odir" unless (-d $odir);
	my $cmd = "set -e\n";
	#my $MAQ = 20; my $BQ = 20; 
	my $mac = 4; my $minFreq = 0.001;
	#my $outF = "$odir/$name[$idx].${LOC}.maf_$minFreq.mac_$mac.cov_20.mq_$MAQ.bq_${BQ}_vt.vanilla.vcf";
	my $reportAllTag = "";
	$reportAllTag = ".allsites" if ($reportAll);
	my $outF = "$odir/$name[$idx].${LOC}.mq_$MAQ.bq_${BQ}_vt$reportAllTag.vanilla.vcf";
	return if (-e $outF);
	my $tmpOut = "$odir/$name[$idx].${LOC}.$MAQ.$BQ.first_pass.vcf";
	my @curReg = ("1");
	if ($doPar || $myPar){
		my $refAR = getRegionsBam($splitFAsize);
		@curReg = @{$refAR};
	}
	my $qsubDir = $odir."/qsubs/";
	system "mkdir -p $qsubDir" unless (-d $qsubDir);
	if ($doPar){@curReg = ("1");}

	my $t = 0; my $threadCollect=0;
	for (my $i=0;$i<@curReg;$i++){
		$cmd = "set -e\n";
		if ($doPar){
			my $regionFile = "$odir/regions_par.txt"; #possibly buggy
			$cmd .= "bash $frDir/fb-parallel $regionFile $ncore --fasta-reference $refFA --bam-list $bamListFiles ";
		} else {
			$cmd .= "$frbBin --fasta-reference $refFA --bam-list $bamListFiles ";
		}
		if ($reportAll){
			$cmd .= "$frAllOpts -m $MAQ -q $BQ ";
		} elsif (0){#old options, single core
			$cmd .= "--pooled-discrete --pooled-continuous --report-genotype-likelihood-max --ploidy 1 \\
			--min-alternate-fraction 0.001 --min-alternate-count 3 --min-mapping-quality 30 --min-base-quality 30 --min-coverage 20 --hwe-priors-off \\
			--allele-balance-priors-off --min-repeat-entropy 1 --haplotype-length 20  ";
		} elsif (0) {
			$cmd .= "--pooled-continuous --pooled-discrete --hwe-priors-off --ploidy 1 --min-alternate-count $mac --min-alternate-fraction $minFreq \\
			--min-mapping-quality $MAQ --min-base-quality $BQ --min-repeat-entropy 1 --use-best-n-alleles 4 ";
		} else { #latest version from seb, "trial"
			$cmd .= " --pooled-continuous --pooled-discrete --hwe-priors-off --min-repeat-entropy 1 --haplotype-length 12 --report-genotype-likelihood-max --min-alternate-count $mac --min-alternate-fraction $minFreq \\
			--min-mapping-quality $MAQ --min-base-quality $BQ --no-partial-observations --ploidy 1";
		}
		if ($myPar){
			$cmd .= " --region '$curReg[$i]'  > $tmpOut.$i \n";
		} else {
			$cmd .= " > $tmpOut\n";
		}
		if (-s $tmpOut){next;}
		#print $cmd."\n";#die;
		#system $cmd;
		if ($doQsub){
			my ($dep,$qcmd) = qsubSystem($qsubDir."frRun$i.sh",$cmd,1,"10G","FR$i","","",1,[],$QSBoptHR);
			push (@allDeps,$dep);
		} else {
			if ($threadCollect){
				$thrs[$t]->join();
			}
			$thrs[$t] = threads->create(sub{system $cmd;});
			$t ++;
			if ($t > $ncore){#collect old threads
				$t=0; $threadCollect=1;
			}
		}
	}

	my $postcmd ="";
	if ($myPar){
		$postcmd .= "cat $tmpOut.* | /g/bork3/x86_64/bin/python /g/bork3/home/hildebra/bin/vcflib/bin/vcffirstheader | /g/bork3/home/hildebra/bin/vcflib/bin/./vcfstreamsort -a | /g/bork3/home/hildebra/bin/vcflib/bin/./vcfuniq > $tmpOut\n\n";
		$postcmd .= "rm $tmpOut.*\n"
	}
	my $postFilter = "| $vcfFilBin -f \"QUAL > 20 & QUAL / AO > 2 & SAF > 1 & SAR > 1 & RPR > 1 & RPL > 1 ";
	$postFilter = "" if ($reportAll);
	$postcmd.="cat $tmpOut $postFilter \\
	| $vtBin normalize -r $refFA -q - 2> /dev/null \\
	| $vcfTools/vcfuniqalleles | $vcfTools/vcfbreakmulti > $outF\n\ntouch $outF.sto";
	#$postcmd .= " $vcfTools/./vcfbreakmulti $tmpOut \\
	#| $vtBin normalize -r $refFA - > $outF\n";
	#/g/korbel/waszak/src/freebayes-v1.0.2-6-g3ce827/bin/freebayes —fasta-reference /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/SebSNP/TEC2/annual/MM3.ctgs.fna —pooled-continuous —pooled-discrete —hwe-priors-off —ploidy 1 —min-alternate-count 2 —min-alternate-fraction 0.001 —min-mapping-
	#quality 10 —min-base-quality 10 —min-repeat-entropy 1 —use-best-n-alleles 4 —bam-list time.series.txt
	#print "\n\nfreebayes post command:\n\n";
	#print $postcmd."\n";
	if ($doQsub){
		my ($dep,$qcmd) = qsubSystem($qsubDir."postCmd.sh",$postcmd,1,"40G","PPfr",join(";",@allDeps),"",1,[],$QSBoptHR);
	} else {
		for (my $t=0;$t<$ncore;$t++){
			#just go over all & collect..
			$thrs[$t]->join();
		}
		system $postcmd;
	}
}


#consensus:
#cat reference.fa | bcftools consensus calls.vcf.gz > consensus.fa


#add later? ------------------------
#| vcffilter -f "QUAL > 20" 
#| vcffilter -f "TYPE = snp"
# --use-best-n-alleles 4 <- speed up
# | vcfkeepinfo - AO RO TYPE
#documentation ------------------------
#haplotypes virus: https://groups.google.com/forum/#!topic/freebayes/BppUf80Oi-s
#freebayes phasing: https://groups.google.com/forum/#!topic/freebayes/fyhho8_H7J0
#freebays exp evo : https://groups.google.com/forum/#!topic/freebayes/Q-TFF8ollC4
#to add bam read groups: https://github.com/ekg/bamaddrg





#this is between sample fixed SNPs call
my $outF="";
my $cmd = "$frbBin \\
--fasta-reference $refFA \\
--report-genotype-likelihood-max \\
--ploidy 1 \\
--min-alternate-fraction 0.2 \\
--min-alternate-count 2 \\
--hwe-priors-off \\
--allele-balance-priors-off \\
--min-repeat-entropy 1 \\
--max-coverage 300 \\
--use-best-n-alleles 4 \\
--bam-list $bamListFiles \\
> raw.vcf

cat raw.vcf \\
| vcffilter \
  -f 'QUAL > 10 & QUAL / AO > 2 & SAF > 1 & SAR > 1 & RPR > 1 & RPL > 1' \
  -s \
| vt normalize -r ref.fasta -q - 2> /dev/null \
| vcfuniqalleles \
| vcfbreakmulti \
> $outF\n";

#/g/korbel/waszak/src/freebayes-v1.0.2-6-g3ce827/bin/freebayes —fasta-reference /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/SebSNP/TEC2/annual/MM3.ctgs.fna —pooled-continuous —pooled-discrete —hwe-priors-off —ploidy 1 —min-alternate-count 2 —min-alternate-fraction 0.001 —min-mapping-quality 10 —min-base-quality 10 —min-repeat-entropy 1 —use-best-n-alleles 4 —bam-list time.series.txt



die;
#-------------- old stuff ----------------
	my $refFa = $refDir."/assemblies/metag/scaffolds.fasta.filt";
	my @bamFiles = ($refDir."/mapping/Align_ment-smd.bam");
	my $outDir = $refDir."/SNVs/global/";

if (!$global){
	$refFa = $ARGV[1];
	$outDir = $ARGV[2];
	
	opendir(DIR, $refDir) or die $!;
	@bamFiles = sort ( grep { /\.bam$/ && -e "$refDir/$_" } readdir(DIR) );	
}
#die "@bamF\n";
system("mkdir -p $outDir");
for (my $i=0; $i<@bamFiles; $i++){
	my $cmd = "";
	print $bamFiles[$i]."\n";
	my $bamF = $refDir."/".$bamFiles[$i];
	my $file = $bamFiles[$i]; $file=~s/-smd\.bam//;
	unless (-e $refFa.".fai"){$cmd .= "$smtBin faidx $refFa\n";}
	unless (-e $bamF.".bai"){$cmd .= "$smtBin index $bamF\n";}

	my $outF = "$outDir/$file";
	next if (-e "$outF.vt.norm.uniq.vcf" && -e "$outF.fb.vcf" && -e "$outF.pileup.gz");
	# -q Qual threshold -R sum of these -C min reads in alternate allele
	#--pooled-continuous / --ploidy 1 : p-c will output all SNPs
	#-kwVa (--no-population-priors --hwe-priors-off --binomial-obs-priors-off --allele-balance-priors-off)
	$cmd .= "$frbBin -f $refFa --min-coverage 10 -C 3  --report-genotype-likelihood-max --pooled-continuous --pooled-discrete --max-complex-gap 50 -j -m 30 -q 20 -= -kwVa -F 0.001 $bamF > $outF.fb.vcf\n";
	#ALTERNATE
	$cmd.= "$vtBin discover -b $bamF -v snps -f 0.001 -e 2 -q 14 -o $outF.vt.vcf -r $refFa -s '1' > $outF.vt.vcf.log\n";
	#tabix file first
	#ALTERNATE 2 pileup[
	$cmd.= "$smtBin mpileup -f $refFa $bamF > $outF.pileup\n";
	#Alternate 3 varscan
	$cmd .= "$varscanBin pileup2snp $outF.pileup --min-coverage 8 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --p-value 1e-02 > $outF.vs2.txt\n";

	#$cmd.="$vtBin normalize LChrom.vcf -r LChrom.fasta  -o LChrom.norm.vcf\n";
	#vt normalize to get clear variants
	$cmd .= $vtBin." normalize -r $refFa -o $outF.fb.norm.vcf $outF.fb.vcf > $outF.normalize.fb.rep\n";
	$cmd .= $vtBin." uniq $outF.fb.norm.vcf -o $outF.fb.norm.uniq.vcf\n";
	$cmd .= $vtBin." normalize -r $refFa -o $outF.vt.norm.vcf $outF.vt.vcf > $outF.normalize.vt.log\n";
	$cmd .= $vtBin." uniq $outF.vt.norm.vcf -o $outF.vt.norm.uniq.vcf\n";
	#peek
	$cmd.= "gzip $outF.pileup\n";
	#die $cmd."\n";
	system $cmd;
}


sub countRdHits(){
	my ($inBam,$MAQ) = @_;
	my $ret =0;
	if (-e "$inBam.$MAQ.rdCnt"){
		$ret = `cat $inBam.$MAQ.rdCnt`;
	} else {
		$ret = `$smtBin view -c -q $MAQ $inBam`;
		open O2,">$inBam.$MAQ.rdCnt" or die "Can't open $inBam.$MAQ.rdCnt\n"; print O2 $ret; close O2;
	}
	chomp $ret;
	
	return $ret;
}
sub correctBamHD(){
	my ($inBam,$smplTag) = @_;
	die "undefined contig marker" if ($smplTag eq "");
	#die if (system "rm $inBam.bai; $smtBin index $inBam" );
	return "$inBam" if (-e "$inBam.hd");
	my $tmp= `$smtBin view -H $inBam `;
	#die "$tmp\n";
	my @XX = split '\n',$tmp;
	open I,">$inBam.sam";
	foreach my $l (@XX){
		if ($l =~ m/^\@SQ/){
			next unless ($l =~ m/SN:$smplTag/);
		}
		print I $l."\n";
	}
	close I;
	#system "$smtBin reheader $inBam.sam $inBam > $inBam.tmp"; 
	system "$smtBin view $inBam >> $inBam.sam; $smtBin  view -b $inBam.sam > $inBam; $smtBin  index $inBam";
	system "rm $inBam.sam;touch $inBam.hd";
	
	#die "$inBam\n";
	return "$inBam";
}

 
sub pilonCons{
	my ($refFa,$bam,$oD,$prefix) = @_;
	my $pilonJar = "";
	my $pBin = "java -Xmx16G -jar $pilonJar";
	my $cmd = "$pBin --genome $refFa --bam $bam --outdir $oD --output $prefix  --vcf --changes --fix snps";
}

sub help{
	print "Script for SNP calling preparatory phase\n";
	print "Arguments:\n";
	print "-map path to map with samples (MATAFILER format)\n";
	print "-inDir path to dir with mappings of metagenomics reads against ref Genome, bam format\n";
	print "-refGenome fasta formated reference Genome\n";
	print "-name string defining name of output\n";
	print "-outDir where to save output\n";
}




