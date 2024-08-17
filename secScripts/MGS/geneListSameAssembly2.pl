#!/usr/bin/env perl
#takes a list of genes and finds contigs in each sample, on which these genes occur
#USAGE: ./geneListSameAssembly2.pl -m [mode] -iS [geneListFile.txt] -G [path to gene catalog] -o [output path]
#./geneListSameAssembly2.pl -m stat -iM /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can.Rhcl.mgs.srt -o /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/Binning/MetaBat/contigs/ -G /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/
#./geneListSameAssembly2.pl -m stat -iM /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Binning/MetaBat2/MB2.clusters -o /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Binning/MetaBat2/perSmpl/ -G /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/
use warnings;
use strict;
use Getopt::Long qw( GetOptions );
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(getDirsPerAssmblGrp readTabByKey systemW readClstrRev readMapS getAssemblPath readFasta);

sub readInfile;
sub getContigsMGS;

my $samBin = getProgPaths("samtools"); #"/g/bork5/hildebra/bin/samtools-1.2/samtools";

my $mods = "stat";
my $inF = ""; 
my $GCd = ""; 
my $oDir = "";
my $multiInF = "";
my $refSmpl =""; 

GetOptions(
	"mode=s" => \$mods,
	"o=s" => \$oDir,
	"iS=s" => \$inF, #single MGS, as used in TEC2
	"iM=s" => \$multiInF, #multiple MGS
	"G=s" => \$GCd,
	"refSmpl=s" => \$refSmpl, #reference sample
);

die "Not enough input args given!\n" if (($inF eq "" && $multiInF eq "") || $oDir eq "" || $GCd eq "");
$oDir .= "/" unless ($oDir =~ m/\/$/);$GCd .= "/" unless ($GCd =~ m/\/$/);
if (!-f $inF && !-f $multiInF){die "Can;t find input gene file $inF$multiInF\nAborting\n";}
my $refSmplDir = "";

my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;

my ($hrD,$hrm) = getDirsPerAssmblGrp($mapF);
my %map = %{$hrm};
my %DOs = %{$hrD};
#my @DOkeys = @{$DOs{16}{SmplID}};print "@DOkeys\n$DOkeys[-1]\n";print "$DOs{16}\n";die;

my $tamocDir = $map{outDir}; $tamocDir =~ s/,.*//;

if ($refSmpl ne ""){
	if (!exists $map{$refSmpl}{wrdir}){die "Can't find refsample $refSmpl in map!\n";}
	$refSmplDir = $map{$refSmpl}{wrdir};
	print "Processing contigs of ".$refSmplDir."\n";
}

$oDir .= "/" unless ($oDir =~ /\/$/);
system "mkdir -p $oDir";

my %totGenes;
#read gene list
readInfile();

my ($hr1s,$hr2s) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",1,\%totGenes);
$hr1s = undef; #my %gene2cl = %{$hr1s};
my %cl2gene = %{$hr2s};
undef %totGenes; #free some mem
 my %GenInMGS;
#my @tmp = keys %cl2gene; #print "$tmp[1] $tmp[1024]\n";
#print $cl2gene{1058}."\n".$cl2gene{"1058"}."\n";
my %Contigs; my $cnt=0;
foreach my $cMGS (sort keys %GenInMGS){
	next if ($cMGS eq "Bin");
	my $oDir1 = $oDir."/$cMGS/";
	if (-d $oDir1){
		system "rm -r $oDir1; mkdir -p $oDir1;";
	} else {
		system "mkdir -p $oDir1" ;
	}
	#reset file
	open S,">$oDir1/SummaryStats.txt" or die $!;
	print S "Smpl\tSmplGrp\tFinalGene2CtgF\t#genes\t#genesAdded\t#FinalGenes\t#FinalGenesAdded\t#MB2Genes\t#Ctgs\t#CtgsFound\t#CtgsNotFnd\t#Ns\n";
	close S;
	my @tgenes = @{$GenInMGS{$cMGS}};
	print "Processing ".scalar @tgenes ." genes of $cMGS .. \n"; 
	my $hr = getContigsMGS(\@tgenes);
	$Contigs{$cMGS} = $hr;
	$cnt++;
	#last if ($cnt ==3 );
}
print "Done pre-processing MGS\n";
undef %cl2gene;


#foreach my $refSamp (@samples){
foreach my $aGr (sort keys %DOs){
	#$aGr = "16"; #DEBUG
	#open this object again each for loop to reduce mem usage...
	my %contigs2Extract; #just store which contigs get extracted, so fasta only need to be read once..
	my %Sline;
	my @smpls = @{$DOs{$aGr}{SmplID}}; my $refSamp = $smpls[-1];
	#my @DOkeys = @{$DOs{16}{SmplID}};#die "@smpls\n@DOkeys\n";

	my $refSamp2 = $refSamp; if (exists(  $map{altNms}{$refSamp}  )){$refSamp2 = $map{altNms}{$refSamp}; }
	if (!exists($map{$refSamp2})){die "can't find sample $refSamp\n";}
	#my $SmplDir = $map{$refSamp2}{wrdir};
	#print "$refSamp2\n$SmplDir\n";
	my $metaGD = getAssemblPath($map{$refSamp2}{wrdir});
	print " --- $aGr - @smpls ---\n";
	if ( !-e "$metaGD/Binning/MB2/$refSamp.cm"){ #-e "$metaGD/Binning/MB2/MeBa.sto" &&
		print "next"; delete $DOs{$aGr}; next;
	}
	my ($hrCM,$hrMC,$hrQ) = readMB2("$metaGD/Binning/MB2/$refSamp");
	print "Read MB2 bins..\n";
	#print "$metaGD\n";
	my %geneCnts = readTabByKey("$metaGD/genePred/genes.per.ctg");
	#go over bins within this sample
	foreach my $cMGS (sort keys %Contigs){
		next unless (exists($Contigs{$cMGS}{$refSamp}));
		my %Conts = %{$Contigs{$cMGS}{$refSamp}};
		my @allCtgs = keys %Conts;
		my %extrCtgs;
		#extend by MB2?
		my $hr = extend2MB2(\%Conts,$refSamp,$hrCM,$hrQ,$hrMC,\%geneCnts);
		my %locMB2ctgs = %{$hr};

		#die;

		
		#grand summary file
		#stats and determine if this contig is even taken..
		#count number of genes in Sample associated to MGS
		my $genesOnCtgs=0; my $addedGenes = 0;my $genesOnCtgsAF=0; my $addedGenesAF = 0; my $frac=0;
		my $MB2added=0;
		#foreach my $spCtg (@allCtgs){
		my $idx =0; my $ctgLine = "";
#		while ($idx <= $#allCtgs){
#			my $spCtg = $allCtgs[$idx]; 
		foreach my $spCtg (@allCtgs){
			my $ctgKey = $refSamp."__".$spCtg."=";
			#print "$ctgKey ";
			#@{$Conts{$refSamp}{$spCtg}} = sort {$a <=> $b}(@{$Conts{$refSamp}{$spCtg}});
			my @tmp = sort {$a <=> $b}(@{$Conts{$spCtg}}); #sort the genes in a contig
			die "Can't find contig $ctgKey in length ref\n" if (!-exists $geneCnts{$ctgKey});
			my $maxGenes = $geneCnts{$ctgKey};
			$frac = @tmp / ($maxGenes+1);
			#if (@tmp != ($maxGenes) ){print STDERR @tmp ." != $maxGenes\n";}
			$genesOnCtgs += scalar @tmp;
			$addedGenes += ($maxGenes+1) - scalar(@tmp) ;
			my $rmTag = "";
			if ($frac < 0.5 && !exists($locMB2ctgs{$ctgKey})){#splice @allCtgs,$idx,1; $rmTag="X";
			} else {
				$MB2added += scalar @tmp if ($frac <0.5);
				$idx ++; $genesOnCtgsAF += scalar @tmp; $addedGenesAF += ($maxGenes+1) - scalar(@tmp);
				#delete $locMB2ctgs{$spCtg} if (exists($locMB2ctgs{$spCtg}));
				$extrCtgs{$ctgKey} = undef;
			}
			$ctgLine .= "$rmTag$ctgKey\t".@tmp."\t" . ($maxGenes+1) ."\t".$frac."\t". join(",",@tmp) . "\n";
			
		}
		foreach my $ct (keys %locMB2ctgs){$extrCtgs{$ct} = undef;}
		print  scalar (@allCtgs)." Contigs, $genesOnCtgs($genesOnCtgsAF) genes, $addedGenes($addedGenesAF) missed genes in $refSamp, $cMGS\n";
		if (scalar @allCtgs == 0 || ($genesOnCtgsAF + $addedGenesAF) == 0 ){next;}
		my $oDir1 = $oDir."/$cMGS/";
		next if ($ctgLine eq ""); #nothing found...
		open O,">$oDir1/$refSamp.ctgsList.txt";
		print O $ctgLine;
		close O;
		#print S "$refSamp\t$qual\t$genesOnCtgs\t$addedGenes\t$genesOnCtgsAF\t$addedGenesAF\t$MB2added\t".scalar(@allCtgs)."\n";
		my $qual = ($genesOnCtgsAF+$addedGenesAF) / scalar(@allCtgs);
		$Sline{$aGr}{$cMGS} = "$refSamp\t$qual\t$genesOnCtgs\t$addedGenes\t$genesOnCtgsAF\t$addedGenesAF\t$MB2added\t".scalar(@allCtgs)."";
		#do this later, so fasta is only read once?
		#actual contig extraction from fasta
		#for (my $i=0;$i<@allCtgs;$i++){$allCtgs[$i] = $refSamp."__".$allCtgs[$i]."=";}
		$contigs2Extract{$aGr}{$cMGS} = \%extrCtgs;
		#my $scaffs = $metaGD . "/scaffolds.fasta";
		#my $cmd = "$samBin faidx $scaffs ". join (" ", @allCtgs) . " > $oDir1/$refSamp.ctgs.fna";
		#print $cmd."\n";
		#system $cmd;
		#die;

	}
	next if (scalar( keys %{$contigs2Extract{$aGr}}) == 0);
	#last;
	#so far only summary files created, automatic filtering is later required..
#} 

print "========== Extracting single Fasta's per sample =========\n";

#collect all fastas from SNP calling..
#foreach my $aGr (sort keys %contigs2Extract){
	#print "$refSamp\n";
	#$aGr = "16"; #DEBUG
	#only need to get assembly dir once
#	my $refSampA = $smpls[-1];
#	my $refSamp2A = $refSampA; if (exists(  $map{altNms}{$refSampA}  )){$refSamp2A = $map{altNms}{$refSampA}; }
#	my $SmplDirA = $map{$refSamp2A}{wrdir};
#	my $metaGD = getAssemblPath($SmplDirA);
	foreach my $refSampB (@smpls){
		my $refSamp2B = $refSampB; if (exists(  $map{altNms}{$refSampB}  )){$refSamp2B = $map{altNms}{$refSampB}; }
		if (!exists($map{$refSamp2B})){die "can't find sample $refSampB\n";}
		my $SmplDir = $map{$refSamp2B}{wrdir};
		#print "$refSampB ";

		my $cD = $map{$refSamp2B}{wrdir};
		my $tarF = $cD."/SNP/contig.SNPc.MPI.fna.gz";
		if (!-e $tarF){print "Can't find $tarF!\n"; next;}
		print "$tarF\n";
		my $Ncount=0;
		#my $tarF2 = $cD."/SNP/proteins.shrtHD.SNPc.MPI.faa.gz";
		my $hr = readFasta($tarF,1,"\\s"); my %FNA = %{$hr};
		foreach my $cMGS (keys %{$contigs2Extract{$aGr}}){
			my %loc = %{$contigs2Extract{$aGr}{$cMGS}};
			foreach my $k (keys %loc){ $k =~m/_L=(\d+)=/; $loc{$k} = $1;}
			my @ctgs = sort {$loc{$b} <=> $loc{$a} } keys %loc;
			my $oDir1 = $oDir."/$cMGS/";
			my $savFNA = "$oDir1/$refSampB.ctgs.fna";
			my $fnd=0; my $notFnd=0;
			my $fasout="";
			foreach my $ct (@ctgs){
				#die "can't find fasta $ct in $tarF\n" unless (exists($FNA{$ct}));
				if (!exists($FNA{$ct})){
					#print "miss: $ct ";
					$notFnd++;
					next;
				}
				my $cpy = $FNA{$ct};
				$fasout .= ">$ct\n$cpy\n";
				$fnd++;
				$Ncount += $cpy =~ tr/N//;
				$Ncount += $cpy =~ tr/n//;
			}
			next if ($fnd==0); #just don't create output for this..
			open O,">$savFNA" or die "Can't open fasta out $savFNA\n";
			print O $fasout;
			close O;
			print "$refSampB:$cMGS: $notFnd/$fnd   ";
			#print to Summary stats
			open S,">>$oDir1/SummaryStats.txt";
			print S "$refSampB\t$Sline{$aGr}{$cMGS}\t$fnd\t$notFnd\t$Ncount\n";
			close S;

			#last;
		}
		print "\n";
	}
}
print "Zipping fastas..\n";
#zipping of fnas..
foreach my $cMGS (sort keys %Contigs){
	my $oDir1 = $oDir."/$cMGS/";
	systemW "gzip $oDir1/*.ctgs.fna";
}

print "Done with extracting contigs\n$oDir\n";
exit(0);


sub extend2MB2{
	my ($hr,$refSamp,$hrCM,$hrQ,$hrMC,$hrGC) = @_;
	my %Conts = %{$hr};
	my @allCtgs = keys %Conts;
	my %ctg2MB2 = %{$hrCM}; my %MB2Q = %{$hrQ}; my %MB2ctgs = %{$hrMC};
	my %geneCnts = %{$hrGC};

	my %MB2cnts; my %MB2ctgcnts;
	foreach my $spCtg (@allCtgs){
		my $ctgKey = $refSamp."__".$spCtg."=";
		if (exists($ctg2MB2{$ctgKey})){
			$MB2cnts{$ctg2MB2{$ctgKey}} += $geneCnts{$ctgKey};
			$MB2ctgcnts{$ctg2MB2{$ctgKey}} ++;
		}
	}
	if (scalar keys %MB2cnts == 0){return {};}
	my @sMB = sort {$MB2cnts{$b} <=> $MB2cnts{$a}} keys %MB2cnts;
	my $sumM=0;foreach(@sMB){$sumM += $MB2cnts{$_};}
	my %locMB2ctgs;
	#print "@sMB  -  $MB2cnts{$sMB[0]} $MB2cnts{$sMB[1]} $MB2cnts{$sMB[2]}\n". $MB2Q{$sMB[0]}{compl} ."\n$MB2Q{0}{compl} $MB2Q{215}{compl}\n";
	if ($MB2cnts{$sMB[0]} > 500 && 
			( ($MB2ctgcnts{$sMB[0]} > scalar(@{$MB2ctgs{$sMB[0]}}) * 0.9  || $MB2cnts{$sMB[0]} > $sumM * 0.9) ) #pretty unique hit rate..
			&& $MB2Q{$sMB[0]}{compl} > 0.5 && $MB2Q{$sMB[0]}{conta}< 0.05 ){
		my @MB2ctg = @{$MB2ctgs{$sMB[0]}};
		my $addedGenes=0; my $totGenes=0;
		foreach my $ct1 (@MB2ctg){
			$ct1 =~ m/^.*__(.*)=$/;
			my $ct2 = $1;
			$locMB2ctgs{$ct1} = undef;
			$totGenes++;
			next if (exists($Conts{$ct2}));
			$addedGenes++;
			#print $ct2." ";
			#push (@allCtgs,
		}
		print "MB2 added $addedGenes / $totGenes contigs\n"; #if ($addedGenes>0);
	}
	return \%locMB2ctgs;
}

sub readInfile{
	if ($inF ne ""){#single MGS mode
		my @tgenes;#genes in GC that are somehow binned
		open I,"<$inF";
		while (my $line = <I>){
			$line =~ s/\r\n?/\n/; chomp $line; 
			push (@tgenes,$line);
		}
		close I;
		if (@tgenes == 1){#only one entry, check for , separation
			@tgenes = split /[,\s+]/,$tgenes[0];
		}
		foreach my $gg (@tgenes){$totGenes{$gg} = undef;}
		$GenInMGS{1} = \@tgenes;
	} else {#read multi MGS guide
		open I,"<$multiInF" or die "Can't open $multiInF\n";
		while (<I>){
			chomp; my @spl = split /\t/;
			my $curMGS = $spl[0];
			my @tgenes = split /[,\s+]/,$spl[1];
			push(@{$GenInMGS{$curMGS}}, @tgenes);
			#ref of what genes are actually used..
		}
		close I;
		foreach my $mgs (keys %GenInMGS){
			foreach my $gg (@{$GenInMGS{$mgs}}){$totGenes{$gg} = undef;}
		}
	}
	print "Read infile, ".scalar(keys(%GenInMGS))." clusters\n";
}



sub getContigsMGS{
	my ($ar) = @_;
	my @tgenes = @{$ar};
	my %Conts ;
	my $geneCnt=0;my $geneCnt1=0;
	foreach my $genNum (@tgenes){
		#MM3__C10_L=109019;_51
		#$genNum = int $genNum;
		die "can't find gene \"$genNum\"\n" unless (exists $cl2gene{$genNum});
		my $genLst = $cl2gene{$genNum};
		foreach my $gen (split (/,/,$genLst)){
			$gen =~ s/^>//; $geneCnt++;
			#print $gen."\t$genNum\n";
			my @splX = split(/__/,$gen); #this is the sample
			my $smpl = $splX[0];
			my @spl = split(/=_/,$splX[1]); #this is the contig
			#print $smpl."\n" unless (exists ($Conts{$splX[0]}));
			push (@{$Conts{$smpl}{$spl[0]}},$spl[1]);
		}
		$geneCnt1++;
	}
	#my @m15= keys %{$Conts{MM15}}; print @m15. " @m15\n";
	print "Done with $geneCnt / $geneCnt1 genes\n";
	return \%Conts;
}


sub readMB2{
	my ($inF) = @_;
	$inF =~ m/(^.*\/)[^\/]+$/;
	my $inD = $1;
	my $inQ = $inF.".cm";
	my %MB1; my %MB2; my %MBq;
	die "Can't find MB2 qual $inQ\n" unless (-e $inQ);
	open I,"<$inF" or die $!;
	while (<I>){chomp; my@spl=split /\t/;
		next if ($spl[1] eq "0");
		$MB1{$spl[0]} = $spl[1];
		push(@{$MB2{$spl[1]}}, $spl[0]);
	}
	close I;
	open I,"<$inQ" or die $!;
	while (<I>){chomp; my@spl=split /\t/;
		next if ($spl[0] eq "Bin Id");
		$MBq{$spl[0]}{compl} = $spl[11];
		$MBq{$spl[0]}{conta} = $spl[12];
		$MBq{$spl[0]}{hetero} = $spl[13];
	}
	return (\%MB1,\%MB2,\%MBq);
}





