#!/usr/bin/env perl
#extracts abundance of all essential marker genes, as well as FMGs
#./extrAllE100GC.pl /g/scb/bork/hildebra/SNP/GCs/GNM3_ABR
#./extrAllE100GC.pl /g/scb/bork/hildebra/SNP/GCs/
use warnings;
use strict;
sub processSubGenes;sub getEgenes;
sub getGeneSeqsSubGenes;
sub MGintoCats;

use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(readMapS systemW readClstrRev readFasta gzipopen getAssemblPath);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive);

my $smtBin = getProgPaths("samtools");#/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $rarBin = getProgPaths("rare");#

die "No input args \n" if (@ARGV < 1);

my $GCd = $ARGV[0];
my $tmpD = $ARGV[1];
my $oldNameFolders = -1;#$ARGV[1];
my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;
print "Reading map(s): $mapF\n";
my ($hrm,$hr2X) = readMapS($mapF,$oldNameFolders);
my %map = %{$hrm};
my @samples = @{$map{opt}{smpl_order}};
my $Nsmpls = scalar(@samples);
$tmpD.="/FMGextras/";
system "mkdir -p $tmpD";
$| = 1;

my $doSubmit = 1;my $submSys = ""; my $qsubDir = "$GCd/LOGandSUB/";
my $QSBoptHR = emptyQsubOpt($doSubmit,"",$submSys);#,"bash"
$QSBoptHR->{qsubDir} = $qsubDir;

print "Extraction of marker genes started..\n";


#my $GCd = "$inD/GeneCatalog/";
#getGeneSeqsSubGenes("FMG");die();

my $cntSkips=0;my $cnt=0;
my $bonSplit = 5;
my %genesE1h; 
my $genesFMGfilesHR = {}; my %genesFMGstreams;
my $genesGTDBfilesHR = {};
my %seenAssembls; 
if (!-e "$GCd/Mattrix.FMG.mat" && !-e "$GCd/Mattrix.FMG.mat.gz"){
	#$genesE1h{1}{gg} = "falk";
	#die $genesE1h{1}{gg};
	print "Reading reference FMG/GTDBs from assemblies..\n";
	foreach my $smpl(@samples){
		$cnt++;
		if ($cnt % $bonSplit == 0 || $cnt == $Nsmpls){
			print "$cnt/$Nsmpls\n" ;
			if ($cnt > 40){$bonSplit = 20;}
			if ($cnt > 250){$bonSplit = 150;}
			if ($cnt > 1050){$bonSplit = 400;}
		}
		my $SmplName = $map{$smpl}{SmplID};
		my $metaGD = getAssemblPath($map{$smpl}{wrdir});
		next if (exists($seenAssembls{$metaGD}));
		$seenAssembls{$metaGD} = 1;

		#print $SmplName."\n";
		#read in ess 100 genes in assembly
		#DEACTIVATE as e100 matrix is not used currently
		if (0){
			my $e1f = "$metaGD/ContigStats/ess100genes/ess100.id.txt";
			if (-f $e1f){
				#open my $I,"<$e1f" or die "Can't open e100 file $e1f\n";
				my ($I,$OK) = gzipopen($e1f,"e100",1);
				while (my $l = <$I>){
					chomp $l;next if (length($l) < 5);
					my @spl = split /\s/,$l;
					$genesE1h{$spl[1]}{$spl[0]}=1; 
				}	close $I;
			} else {print "Can't find e100 file $e1f\n";}
		}
		#read in FMG genes in assembly
		my $FMGf = "$metaGD/ContigStats/FMG/FMGids.txt";
		print STDERR " Extracting FMG genes..\n";
		$genesFMGfilesHR = MGintoCats($FMGf, $genesFMGfilesHR);
		print STDERR " Extracting GTDB genes..\n";
		my $GTDBf = "$metaGD/ContigStats/GTDBmg/marker_genes_meta.tsv";
		$genesGTDBfilesHR = MGintoCats($GTDBf, $genesGTDBfilesHR);
		
	}
	print "Read ref dataset\n";
	print "skipped $cntSkips samples\n" if ($cntSkips > 0);
}
#foreach my $cat (keys %genesFMGstreams){
	#close $genesFMGstreams{$cat};
#}
if (!-e "$GCd/Mattrix.FMG.mat" && !-e "$GCd/Mattrix.FMG.mat.gz"){
	
	my $allGs = {};
	$allGs = getEgenes($genesFMGfilesHR,$allGs);
	$allGs = getEgenes($genesGTDBfilesHR,$allGs);
	
	print STDERR "Subset of gene cats to read: " . scalar(keys %{$allGs}) . "\n";
	
	print STDERR "Reading cluster index\n";
	my ($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",2,$allGs);
	print STDERR "Done reading idx\n";

	#print "@e1cat\n";
	#Sort matrix genes into FMGs and extract matrix subset that contains these genes
	print "Creating FMG gene matrix\n";
	processSubGenes($genesFMGfilesHR,"FMG",$hr1);
	#Sort matrix genes into FMGs and extract matrix subset that contains these genes
	print "Creating GTDB gene matrix\n";
	processSubGenes($genesGTDBfilesHR,"GTDBmg",$hr1);
	undef $hr1 ;
	#don't create for e100
	#print "Creating e100 gene matrix\n";
	#processSubGenes(\%genesE1h,"e100");
}
#get actual FNAs & FAAs of all FMG genes
print "extracting FNA & FAA's of FMG genes\n";
getGeneSeqsSubGenes("FMG");
print "extracting FNA & FAA's of GTDB genes\n";
getGeneSeqsSubGenes("GTDBmg");


print "Done FMG\n";
print "Finished\n";
exit(0);



sub MGintoCats{
	my ($FMGf, $hr) = @_;
	my %genesFMGfilesL = %{$hr};
	if (-f $FMGf){
		#open my $I,"<$FMGf" or die "Can't open FMG file $FMGf\n";
		#print "$FMGf\n";
		my %genesFMG; #local for each file
		my $gCnt=0;
		my ($I,$OK) = gzipopen($FMGf,"FMG infile",1);
		while (my $l = <$I>){
			chomp $l;next if (length($l) < 5);
			my @spl = split /\s+/,$l;
			$genesFMG{$spl[1]}{$spl[0]}=1; 
			$gCnt++;
		}	
		close $I;
		#print "N=$gCnt ";
		#store between runs..
		foreach my $cat (keys %genesFMG){
			if (!exists($genesFMGfilesL{$cat})){
				my $tmpF = "$tmpD/cat.$cat.idx";
				#print STDERR "Deleting $tmpF\n";
				system "rm -f $tmpF";
				$genesFMGfilesL{$cat} = $tmpF;
				#open ($genesFMGstreams{$cat},">$tmpF") or die "can't open $tmpF\n";
			}
			open OOX,">>",$genesFMGfilesL{$cat} or die "Can't opebn write tmp $genesFMGfilesL{$cat}\n";
			foreach my $ge (keys %{$genesFMG{$cat}}){
				print OOX "$ge\n";
			}
			close OOX;
		}
	} else {print "Can't find FMG file $FMGf\n";}
	return \%genesFMGfilesL;
}

sub getGeneSeqsSubGenes(){
	my ($tag) = @_;
	my $subF = "$GCd/$tag.subset.cats";
	my $fmgOD = "$GCd/$tag/";
	#die "TODO getGeneSeqsSubGenes\n";
	system "mkdir -p $fmgOD"; 
	my $hr = readFasta("$GCd/compl.incompl.95.fna",1);
	my %FNA = %{$hr};
	
	print "Read FNA \n";
	
	open I,"<$subF" or die "can't open $subF\n"; 
	while (my $line=<I>){
		chomp $line;
		my @spl = split /\t/,$line;
		my @spl2 = split /,/,$spl[2];
		#die "\n@spl2\n";
		my $ofile = $fmgOD."/$spl[0]";
		if (0){
			system "$smtBin faidx $GCd/compl.incompl.95.fna ". join (" ", @spl2) . " > $ofile.fna";
		} else {
			open O1,"> $ofile.fna" or die $!;
			foreach my $k (@spl2){
				print O1 ">".$k."\n$FNA{$k}\n";
			}
			close O1; 
		}
	} 
	close I;
	%FNA = ();
	
	#split to lessen mem burden
	$hr = readFasta("$GCd/compl.incompl.95.prot.faa",1);
	my %FAA = %{$hr};
	print "Read FAA \n";
	open I,"<$subF" or die "can't open $subF\n"; 
	while (my $line=<I>){
		chomp $line;
		my @spl = split /\t/,$line;
		my @spl2 = split /,/,$spl[2];
		#die "\n@spl2\n";
		my $ofile = $fmgOD."/$spl[0]";
		if (0){
			system "$smtBin faidx $GCd/compl.incompl.95.prot.faa ". join (" ", @spl2) . " > $ofile.faa";
		} else {
			open O2,"> $ofile.faa" or die $!;
			foreach my $k (@spl2){
				print O2 ">".$k."\n$FAA{$k}\n";
			}
			close O2;
		}
	} 
	close I;


}


sub getEgenes{
	my ($ghr,$r1) = @_;
	my %genesF = %{$ghr};
	my @e1cat = keys %genesF;

	foreach my $e1c (@e1cat){
		#my @spG = keys %{$genes{$e1c}};#my @spG;
		open I,"<$genesF{$e1c}" or die "cant open tmp file $genesF{$e1c}\n";
		while (<I>){chomp; ${$r1}{$_}=1;}
		close I;
	}
	return $r1;
}

sub processSubGenes{
	my ($ghr,$tag,$gene2cl) = @_;
	my %genesF = %{$ghr};
	#my %gene2cl = %{$hr1};
	my $subF = "$GCd/$tag.subset.cats";
	my @e1cat = keys %genesF;
	my %selC; open O ,">$subF" or die "Can't open $subF";
	my $catCnt=0;my $geneCnt=0;my $fail_abort=0;
	foreach my $e1c (@e1cat){
		#my @spG = keys %{$genes{$e1c}};
		my @spG;
		open I,"<$genesF{$e1c}" or die "cant open tmp file $genesF{$e1c}\n";
		while (<I>){chomp; push @spG,$_;}
		close I;
		#die "N=".@spG."  $spG[0]   $genesF{$e1c}\n";
		my %selD;   #print "\n";
		foreach my $gene (@spG){
			unless (exists ${$gene2cl}{">".$gene}){#can also be gene too short..
				print $gene. " missing from GC\n" ;
				$fail_abort++;
				next;
			}
			#print $gene."\n";
			$selD{$$gene2cl{">".$gene}} = 1;
		}
		my @keysHds = keys %selD;
		my $size = @keysHds;
		if ($size == 0){die "@spG"."\n";}
		my $addD = $keysHds[0];
		if ($size > 1 ){$addD = join(",",sort(@keysHds));}
		#print $addD  ."    @keysHds\n";
		print O $e1c."\t$size\t".$addD."\n";
		$geneCnt += $size;
		my %tmp = (%selC, %selD); %selC = %tmp;
		$catCnt++;
	#die ( @spG."\n");
	}
	close O;
	if ($geneCnt <= 40 ){die "Error in extrAllE100GC.pl: found only $geneCnt marker genes!\n";}
	#my $sedStr = join("p;",@rows);
	#system O "sed -n '1p;$sedStr"."p' $GCd/Matrix.mat > $GCd/e100subset.mat";
	my @rows = keys %selC;
	if ($fail_abort){print "Missing $fail_abort essential genes.. aborting\n";}
	push(@rows,1);
	print "Selected ".@rows." rows.. writing to $GCd/$tag.subset.mat\n";
	open O,">$GCd/$tag.lines" or die "Can't open output lines file\n";;print O join("\n",@rows);close O;
	@rows = ();
	#%gene2cl = ();
	my $pigzBin = getProgPaths("pigz");
	my $cores = 6;
	my $cmdX= "$rarBin lineExtr -i $GCd/Matrix.mat.gz -o $GCd/Matrix.$tag.mat -reference $GCd/$tag.lines -t $cores\n "; #-checkRowName2Idx 
	$cmdX .= "$rarBin lineExtr -i $GCd/Mat.cov.mat.gz -o $GCd/Mat.cov.$tag.mat -reference $GCd/$tag.lines -t $cores\n ";
	$cmdX .= "$rarBin lineExtr -i $GCd/Mat.med.mat.gz -o $GCd/Mat.med.$tag.mat  -reference $GCd/$tag.lines -t $cores\n ";
	$cmdX .= "$pigzBin -p 6 $GCd/Mat.cov.$tag.mat $GCd/Mat.med.$tag.mat\n";
	$cmdX .= "rm $GCd/$tag.lines\n";
	print $cmdX."\n";
	#die $cmdX."\n";
	#systemW $cmdX;
	#die "Couldn't create gene abundance matrix $GCd/Matrix.$tag.mat , fatal error\n" unless (-e "$GCd/Matrix.$tag.mat");
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
	my ($dep1,$qcmd1) = qsubSystem($qsubDir."Matrix.sub.$tag.sh",$cmdX,$cores,int(100/$cores)."G","extE100$tag","","",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} =$tmpSHDD;

	
	#system "rm $GCd/$tag.lines";
	print "Submitted $tag\n";
	print "Found $catCnt gene categories, with total of $geneCnt members\n";
}
