#!/usr/bin/perl
#builds from a contig file several stats to seperate contigs into single species
#$subparts =~ a=gene abundance g=GC e=essential k=kmer m=metabat s=microsats
#./separateContigs.pl /g/scb/bork/hildebra/SNP/GNMass/alien-11-374-0/
use warnings;
use strict;
#use Scalar::Util qw(looks_like_number);
use Getopt::Long qw( GetOptions );


sub readFasta;
sub readGFF;
sub geneAbundance; sub runMaxBin;
sub findMicrSat;
use Mods::GenoMetaAss qw(systemW is_integer reverse_complement_IUPAC);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::phyloTools qw( getE100);


die "Not enough input args\n" if (@ARGV < 2);
my $inD =  "";#$ARGV[0];
my $assD = "";#$ARGV[1];
my $subparts = "";#$ARGV[2];
my $readLength = 150;#$ARGV[3]; #primary read length
my $readLengthSup = 9000; # support reads length, eg PacBio would be ~9000
my $tmpD = "";#$inD."/tmp/";$tmpD =~ s/\$/\\\$/g;$tmpD = $ARGV[4] if (@ARGV > 4);#die "$tmpD\n";
my $Nthreads = 1; #$Nthreads = $ARGV[5] if (@ARGV > 5);


GetOptions(
	"threads=i" => \$Nthreads,
	"inD=s" => \$inD,
	"assD=s" => \$assD,
	"tmpD=s" => \$tmpD,
	"subparts=s" => \$subparts,
	"readLength=i" => \$readLength,
	"readLengthSup=i" => \$readLengthSup,
	
);
if ($tmpD eq ""){
	$tmpD = $inD."/tmp/";
	$tmpD =~ s/\$/\\\$/g;
}

my $rdCovBin =getProgPaths("readCov");
my $pigzBin = getProgPaths("pigz");

my $inScaffs = "$assD/scaffolds.fasta.filt";
my $proteins = "$assD/genePred/proteins.shrtHD.faa";
my $genesNT = "$assD/genePred/genes.shrtHD.fna";
my $genesGFF = "$assD/genePred/genes.gff";
my $genesAA = "$assD/genePred/proteins.shrtHD.faa";
my $rawINPUTf = $inD."input.txt";
my $outDab = $inD."/assemblies/metag/ContigStats/";
my $outD = $assD."/ContigStats/";
my $cmd = "";
my $cleanUp = 0;#1=do overwrite
if ($cleanUp){
	system "rm -r $outD";
}
system "mkdir -p $outD";
system "mkdir -p $outDab";

#figure out mapping names
die "Mapping is not done yet (or not copied)\n" if (!-e "$inD/mapping/done.sto");
my $SmplNm = `cat $inD/mapping/done.sto`;
$SmplNm =~ s/-smd.bam\n?//;
#my $inBAM = $inD."mapping/$SmplNm-smd.bam";
#die "$coverage\n";
#system "cp $inScaffs $outD";
######################### ini calcs
#$cmd = "";
#$cmd .= "cut -f1 -d \" \" $proteins > $proteins.shrtHD.faa\n" unless (-e "$proteins.shrtHD.faa");
#$cmd .= "cut -f1 -d \" \" $genesNT > $genesNT.shrtHD.fna\n" unless (-e "$genesNT.shrtHD.fna");
#systemW $cmd unless ($cmd eq "");
###############################   Coverage per gene  ###############################


my $coverage = $inD."mapping/$SmplNm-smd.bam.coverage";
geneAbundance($coverage,0,$readLength) if ($subparts =~ m/a/);
my $covSup = $inD."mapping/$SmplNm.sup-smd.bam.coverage";
geneAbundance($covSup,1,$readLengthSup) if ($subparts =~ m/a/ && (-s $covSup || -s "$covSup.gz" ) );
#print "$covSup\n";

###############################   GC content  ###############################
if ((!-e "$outD/scaff.GC.gz" || !-e "$outD/scaff.pergene.GC.gz"|| !-e "$outD/scaff.pergene.GC3.gz") && $subparts =~ m/g/){
	if (-s "$outD/scaff.GC" && -s "$outD/scaff.pergene.GC" && -s "$outD/scaff.pergene.GC3"){
		systemW "$pigzBin -p $Nthreads $outD/scaff.*";
	} else {
		print "Analyzing GC content\n";
		my $GCP = getProgPaths("calcGC_scr");
		systemW "$GCP $inScaffs $outD/scaff.GC";
		systemW "$GCP $genesNT $outD/scaff.pergene.GC $genesGFF"; 
		systemW "$pigzBin -p $Nthreads $outD/scaff.*";
		#systemW $cmd."\n";
	}
	
} elsif ($subparts =~ m/g/) {
	print "GC content already calculated\n";
}

###############################   essential GTDB proteins  ###############################
my $oDess = "$outD/ess100genes/";my $outDFMG = "$outD/FMG/"; my $outDGTDB = "$outD/GTDBmg/";
if ($subparts =~ m/G/ && !-e "$outDGTDB/marker_genes_meta.tsv"){ #GTDB markers, Jogi
	system("mkdir -p $outDGTDB $tmpD/extract_gtdb_mg_tmp/");
	#$FMGd/extract_gtdb_mg.py
	my $GTDBmgScr =  getProgPaths("GTDBmg_scr");
	#python /hpc-home/hildebra/dev/python/extract_gtdb_mg/extract_gtdb_mg.py  --threads 4 --aa /hpc-home/hildebra/grp/data/projects/postDrama//AssmblGrp_s27/metag//genePred/proteins.shrtHD.faa --overwrite_tmp --tmp /nbi/local/ssd/50677730/MATAFILER//IN217/extract_gtdb_mg_tmp/ -o /hpc-home/hildebra/grp/data/projects/postDrama//AssmblGrp_s27/metag//ContigStats//GTDBmg// --list_only
	$cmd = "$GTDBmgScr  --threads $Nthreads --aa $proteins --nt $genesNT --overwrite_tmp --tmp $tmpD/extract_gtdb_mg_tmp/ --list_only -s --meta -o $outDGTDB/\n";
	#$cmd .= "mv $outDGTDB/marker_genes.tsv $outDGTDB/GTDBids.txt\n";
	print $cmd;
	systemW $cmd;
	system "touch $outDGTDB/marker_genes_meta.tsv" if (!-e "$outDGTDB/marker_genes_meta.tsv");
} elsif ($subparts !~ m/G/) {
	print "No GTDB core genes requested\n";
} else {
	print "GTDB MG already exists: $outDGTDB/marker_genes_meta.tsv\n";
}


#if ($subparts =~ m/eF/ && (!-e "$outDFMG/FMGids.txt" || !-e "$oDess/e100split.sto")){
	###############################   fetchMG  ###############################
if ( $subparts =~ m/F/  && !-e "$outDFMG/FMGids.txt"){
	my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
	my $FMGrwkScr = getProgPaths("FMGrwk_scr");
	print "Searching for fetch MG core proteins in predicted genes..    ";
		system("mkdir -p $outDFMG");
		$cmd = "perl $FMGd/fetchMG.pl -m extraction -o $outDFMG -l $FMGd/lib -t $Nthreads -d $genesNT $proteins"; # -x $FMGd/bin  -b $FMGd/lib/MG_BitScoreCutoffs.uncalibrated.txt
		$cmd .= "; $FMGrwkScr $outDFMG";
		if ( !-s "$outDFMG/COG0012.faa" && !-s "$outDFMG/COG0016.faa"){
			systemW $cmd;
		}
		system "rm  -rf $assD/genePred/*.cidx $outDFMG/COG0* $outDFMG/temp $outDFMG/hmmResults";
		system "touch $outDFMG/FMGids.txt" unless (-e "$outDFMG/FMGids.txt");
} elsif ($subparts !~ m/F/) {#(!-e "$outDFMG/FMGids.txt") {
	print "No FetchMG essential proteins requested\n";
} else {
	print "FetchMG essential proteins already exists: $outDFMG/FMGids.txt\n";
}

	############################### essential 100 proteins ###############################
if ( $subparts =~ m/E/ && !-e "$oDess/e100split.sto"){
	getE100($oDess,$proteins,$genesNT,$Nthreads);
	sleep (3);
	systemW("touch $oDess/e100split.sto;");
	#required for maxbin
	#rm $oDess/assembly.hmm*");
	#die ($protIDs[0]."\n");

	print "Done\n";
}

###############################   growth prediction of sample  ###############################
if ($subparts =~ m/g/){
	#currently deactivated..
	#my $growthBin =getProgPaths("growthP");
	#systemW("cat $oDess/*.fna >$oDess/alle100.fna") unless (-e "$oDess/alle100.fna");
	#$cmd = "$growthBin -f $oDess/alle100.fna -g $genesNT -c 0 -T $Nthreads -m";
}

###############################   kmer content  ###############################
#print $cmd."\n";
if (!-s "$outD/scaff.4kmer.gz" && $subparts =~ m/k/){
	print "Analyzing kmer content\n";
	my $KCPbin = getProgPaths("kmerFreqs");#"perl /g/bork5/hildebra/bin/multi-metagenome-master/R.data.generation/calc.kmerfreq.pl";
	systemW "$KCPbin -i $inScaffs -m 100 -o $outD/scaff.4kmer";
	systemW "$pigzBin -p $Nthreads  $outD/scaff.4kmer";
}
if (!-s "$outD/scaff.pergene.4kmer.gz" && $subparts =~ m/4/){
	my $KCPbin = getProgPaths("kmerFreqs");#"perl /g/bork5/hildebra/bin/multi-metagenome-master/R.data.generation/calc.kmerfreq.pl";
	systemW "$KCPbin -i $genesNT -m 100 -o $outD/scaff.pergene.4kmer";
	systemW "$pigzBin -p $Nthreads  $outD/scaff.pergene.4kmer\n";
}
if (!-s "$outD/scaff.pergene.4kmer.pm5.gz" && $subparts =~ m/4/){
	my $KperGenes = getProgPaths("kmerGeneWin");
	systemW "$KperGenes $outD/scaff.pergene.4kmer.gz 5\n";
}
sleep(1);
###############################   binning  ###############################
if ( !-s $inD."Binning/MaxBin/MB.summary" && $subparts =~ m/m/){
	#runMaxBin(); #outdated, don't use any longer..
}
if ( ( !-s "$inD/Binning/MetaBat/MeBa.sto" || !-s "$inD/Binning/MetaBat/$SmplNm.cm" ) &&  $subparts =~ m/m/){
#deactivated, needs to be in MATAFILER routine to take assembly groups into account
#	my $compoundBinningScr = getProgPaths("cmpBinScr");#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/compoundBinning.pl";

#systemW "perl $compoundBinningScr $inD $tmpD";
}

if (-e "$outD/microsat.txt" && $subparts =~ m/s/ && int(-s "$outD/microsat.txt") == 0 ){
	findMicrSat($inScaffs,"$outD/microsat.txt");
}

print "all done\n";
exit(0);

sub findMicrSat{
	my ($scaffs,$rep) = @_;
	my $mrepsB = getProgPaths("mreps");#"/g/bork3/home/hildebra/bin/mreps/./mreps";
	my $cmd = "$mrepsB -fasta -minperiod 2 -maxperiod 13 -maxsize 200 $scaffs > $rep\n";
	#die $cmd;
	print "Searching for microsattelites..";
	systemW $cmd;
	print "Done\n";
}



sub runMaxBin {
	print "Running MaxBin..\n";
	my $maxBin = getProgPaths("maxBin");#"perl /g/bork5/hildebra/bin/MaxBin-1.4.2/run_MaxBin.pl";
	my $rwkMB = getProgPaths("maxBinFrmt");# "perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/maxBin_rewrk.pl";
	my $covCtgs = "$outDab/Coverage.median.percontig";
	my $essHMM = "$outD/ess100genes/assembly.hmm.orfs.txt";
	my $outF = "MB";
	my $outD2 = $inD."Binning/MaxBin/";

	print "Binning contigs with MaxBin\n";
	system ("mkdir -p $outD2");

	my $cmd  = "$maxBin -contig $inScaffs -out $outD2/$outF -abund $covCtgs -thread $Nthreads -HMM $essHMM -AA $genesAA -min_contig_length 400\n";
	$cmd .= "$rwkMB $outD2 $outF";
	system $cmd;
}

sub calcGeneCov($ $ $){
	my ($arr,$hr,$curChr) = @_;
	#print "geneCov\n";
	my %gff = %{$hr}; my $retStr = "";
	my @cov = @{$arr};
	my $cnt=1;
	unless (exists($gff{$curChr})){
		print "$curChr chromosome doesn't exist\n" ;
		return "";
	}
	foreach my $ge (keys %{$gff{$curChr}} ){
		my $st = $gff{$curChr}{$ge}{start}; my $en = $gff{$curChr}{$ge}{stop};
		if ($en > @cov){print "$en in $curChr ".@cov."\n";}
		my $gcov = 0; for ($st .. $en){ $gcov += $cov[$_];}
		$gcov /= ($en-$st);
		#die $gcov."\n";
		$retStr .= $curChr."_$cnt\t$gcov\n";
		$cnt++;
	}
	#die();
	return $retStr;
}

sub geneAbundance{
	my ($inF, $isSupport, $readL) = @_;
	#die "$inF\n";
	#my $hr = readGFF($inD."assemblies/metag/genePred/genes.gff");
	my $outF = $inF . ".pergene"; 	my $outF2 = $inF . ".percontig";	
	my $outF3 = $inF . ".window";	my $outF4 = $inF . ".geneStats";
	my $outF5 = $inF . ".count_pergene"; my $outF6 = $inF . ".median.pergene";
	my $outF7 = $inF . ".median.percontig";
	my $fileEnd2 = "";
	my $oPrefix = "Coverage";
	$oPrefix = "Cov.sup" if ($isSupport);
	my $outFfin = $outDab . "$oPrefix.pergene$fileEnd2"; 	my $outF2fin = $outDab . "$oPrefix.percontig$fileEnd2";	
	my $outF3fin = $outDab . "$oPrefix.window$fileEnd2";	my $outF4fin = $outDab . "GeneStats.txt";	
	my $outF5fin = $outDab . "$oPrefix.count_pergene$fileEnd2";	my $outF6fin = $outDab . "$oPrefix.median.pergene$fileEnd2";
	my $outF7fin = $outDab . "$oPrefix.median.percontig$fileEnd2";	
	my $stone = "$outDab/$oPrefix.stone";
	#print $stone."\n\n";
	my $inFG = $inF.".gz";
	if ( (-e $outF7fin || -e "$outF7fin.gz") && (-s $outFfin || -s "$outFfin.gz") && -e $stone){
		print "Coverage $isSupport was already calculated in $outDab\n";
		#some cleanup operations.. good to run, if already here..
		if (-s $inFG && -s $inF){system "rm -f $inF";}
		if (-s $outFfin && !-s "$outFfin.gz"){
			print "gzipping existing output..\n";
			systemW "$pigzBin -p $Nthreads $outDab/${oPrefix}.*";
		}
		#print "$outFfin\n";
		return;
	}
	system "echo \"$assD\" > $inD/assemblies/metag/assembly.txt" unless (-e "$inD/assemblies/metag/assembly.txt");
	print "Calculating coverage of assemblies..\n";
	if (-s $outFfin && -s $outF2fin&& -s $outF3fin&& -s $outF4fin && -s $outF5fin && -s $outF6fin){print "Gene abundance was already calculated\n";return;}
	my $clnCmd = "";
	if (!-e $inFG && -s $inF){system "$pigzBin -p  $inF $Nthreads";}
	if (-e $inFG && !-s $inF){system "rm -f $inF";}
	
	#no longer needed, rdCov can also read in .gz 
	#if (-e $inFG && !-e $inF){systemW("gunzip -c $inFG > $inF"); $clnCmd .= "rm -f $inF\n";
	#} elsif (-e $inF && !-e $inFG){$clnCmd .= "$pigzBin -p $Nthreads $inF\n";}
	
	#open I,"<$inF" or die "Can't open $inF\n";	open O,">$outF" or die "Can't open output $outF\n";	open O2,">$outF2" or die "Can't open output $outF2\n";
	#1st delinearize	my $curChr = ""; my $cont=1; my $fresh = 1; my $expLength = 0;	my @chrCov=();	my $start = time; my $ctgCnt = 0; my $ccov = 0; #contig coverage
	#while (my $line = <I>){		chomp $line; 		my @spl = split(/\t/,$line);		my $chr = $spl[0];		if ($curChr ne $chr ){			unless  ($curChr eq ""){				#print "1";				$curChr =~ m/_length_(\d+)_/;				#if (length(@chrCov)<= $1){@chrCov = (@chrCov,(0)x($1-@chrCov+1)); }								#for (0 .. (@chrCov-1)){ $ccov += $chrCov[$_];} $ccov /= (@chrCov-1);				print O2 "$curChr\t$ccov\n"; $ccov = 0;				#print @chrCov." @chrCov $1\n";				#next if ($1 < 500); #because no gene calls on these				print O calcGeneCov(\@chrCov,$hr,$curChr);				$ctgCnt++;				if ($ctgCnt == 500){					my $duration2 = time - $start;					die "Time for 500 : ". $duration2." s\n";				}			}			#die $outF."\n";			$curChr = $chr;  $fresh = 1;			$curChr =~ m/_length_(\d+)_/;			$expLength = $1; @chrCov = ((0)x($expLength+1));			#my @test = (int($spl[1])..(int($spl[2])-1)); print "@test X $spl[3]\n";		}		#print $spl[3]."\n";		#if ($fresh){$fresh = 0;	if ($spl[1] != 0){@chrCov = (@chrCov,(0)x($spl[1]));}}
		#@chrCov = push(@chrCov,($spl[3])x($spl[2]-$spl[1]));		my $fill = $spl[3];		#$ccov += $spl[3] * ($spl[2] - $spl[1]-1);		foreach ($spl[1] .. ($spl[2]-1)){ $chrCov[$_] = $fill }		#if (!is_integer($spl[2])){die $spl[2]."  2 \n";}		#if (!is_integer($spl[1])){die $spl[1]."  1 \n";}		#if (!is_integer($spl[3])){die $spl[3]."  3 \n";}		#$chrCov[$spl[1]..$spl[2]-1] = $spl[3];		#$chrCov[11..22] = int($spl[3]);				#die @chrCov." @chrCov\n" if ($spl[3] == 6);	}
	#close I; close O; close O2;
	#my $duration = time - $start;
	#print "Gene abundances calculated in $duration s\n";
	#die ("$clnCmd");
	
	#print "$rdCovBin $inF $assD/genePred/genes.gff $readLength";
	my $cmd = "rm -f ${inFG}.*;$rdCovBin $inFG $assD/genePred/genes.gff $readL;\n";
	systemW $cmd;
	#$cmd .= "gzip -f $outF $outF2 $outF3\nmv $outF.gz $outFfin\n mv $outF2.gz $outF2fin\nmv $outF3.gz $outF3fin\n";
	
	$cmd = "rm -f ${inFG}.*.gz;\n$pigzBin -p $Nthreads $outF $outF2 $outF3 $outF4 $outF5 $outF6 $outF7; \n"; $fileEnd2 = ".gz";
	systemW $cmd;
	
	$cmd =  "mv $outF$fileEnd2 $outFfin$fileEnd2\nmv $outF2$fileEnd2 $outF2fin$fileEnd2\nmv $outF3$fileEnd2 $outF3fin$fileEnd2\nmv $outF4$fileEnd2 $outF4fin$fileEnd2\nmv $outF5$fileEnd2 $outF5fin$fileEnd2\n";
	$cmd .= "mv $outF6$fileEnd2 $outF6fin$fileEnd2\nmv $outF7$fileEnd2 $outF7fin$fileEnd2\n";
	#die $cmd;
	systemW $cmd;
	
	#die $cmd."\n";
	#systemW $cmd;
	#print "cont";
	systemW($clnCmd);
	systemW "touch $stone\n";
	#print $stone."\n\n";
}

sub readGFF($){
	my ($inF) = @_;
	my %gff;
	open I,"<",$inF or die "Can't find ".$inF."\n";
	my $curChrCnt = 0;
	while (my $line = <I>){
		next if ($line =~ m/^#/);
		chomp $line; 
		#die $line."\n";
		my @spl = split (/\t/,$line);
		my $chr = $spl[0];
		if (exists($gff{$chr})){
			$curChrCnt ++;
		} else {
			$curChrCnt = 0;
		}
		#die $spl[3]."\n";
		$gff{$chr}{$curChrCnt}{start} = $spl[3];
		$gff{$chr}{$curChrCnt}{stop} = $spl[4];
	}
	close I;
	print "\nFinished reading gff\n";
	return \%gff;
}

