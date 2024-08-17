#!/usr/bin/env perl
#get FMG genes within a single MB run from a specific sample (1st arg), specific Bin number (arg 2), which can be named (ARG3) and get this from all subfolders folders (arg4) to output folder (arg 5) using the GeneCatalog folder (arg 6)
#run modes: 
#		"allGeneIDs" - gets all the genes associated to TEC and an abundance matrix thereof
#		"singleSample" - only extract fasta seqs for each FMG/e100; does not use GC & doesn't build phylo tree
#		"FMGcrossSmpls" - normal run mode extracting fastas of FMGs/e100 for each sample where present and builds a *tree* in the end
#				has an extra option %id as last arg; if provided, scans for relatives of each FMG within %id radius and adds to tree
#				also does "allGeneIDs" functionality now, extracting the complete matrix
#		"contigsOnly" - extracts the corresponding contigs for a TEC; pretty simple task
#Examples
#./getMarkersMultiSmpl.pl FMGcrossSmpls alien-11-380-0/ 2 TEC6 /g/scb/bork/hildebra/SNP/GNMass3/ /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4.1/TEC6/ /g/scb/bork/hildebra/SNP/GCs/GNM3_ABR
#./getMarkersMultiSmpl.pl onlyIDs alien-11-376-0/ 2 TEC2 /g/scb/bork/hildebra/SNP/GNMass3/ /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/ /g/scb/bork/hildebra/SNP/GCs/T2_GNM3_ABR
#./getMarkersMultiSmpl.pl allGeneIDs alien-11-376-0/ 2 TEC2 /g/scb/bork/hildebra/SNP/GNMass3/ /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/ /g/scb/bork/hildebra/SNP/GCs/T2_GNM3_ABR 90
#./getMarkersMultiSmpl.pl FMGcrossSmpls simulated_metaG3SA_2 2 TEC2 /g/scb/bork/hildebra/SNP/SimuA/ /g/scb/bork/hildebra/SNP/SimuA/xtest 
#./getMarkersMultiSmpl.pl singleSample alien-11-374-0/ 1 TEC2 /g/scb/bork/hildebra/SNP/GNMass3/ 
#./getMarkersMultiSmpl.pl specificIDs /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/specific_genes/Tec2.ABR.txt - ABR /g/scb/bork/hildebra/SNP/GNMass3/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/specific_genes/ABR/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/
#./getMarkersMultiSmpl.pl specificIDs /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/R_filt/T2_4_tree.txt A TR /g/scb/bork/hildebra/SNP/GNMass3/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/specific_genes/T2_tree/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/
#./getMarkersMultiSmpl.pl specificIDs /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/ABR_all/ABR.list S ABRa /g/scb/bork/hildebra/SNP/GNMass3/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/specific_genes/ABR_all/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/
#./getMarkersMultiSmpl.pl specificIDs /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T6/R_filt/T6_4tree.txt A TR /g/scb/bork/hildebra/SNP/GNMass3/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/specific_genes/T6_tree/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/
#./getMarkersMultiSmpl.pl specificIDs /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T3/R_filt/T3_4tree.txt A TR /g/scb/bork/hildebra/SNP/GNMass3/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/specific_genes/T3_tree/ /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/
#866196,2029841,429290


use warnings;
use strict;

use Mods::GenoMetaAss qw(readGFF readClstrRev reverse_complement_IUPAC readMapS renameFastHD convertNT2AA);
use Mods::IO_Tamoc_progs qw(getProgPaths);


sub extractGenesPerSample;
sub extractFNAFAA2genes;
sub readFasta;
sub buildTree;

my $buildTees = 0;

my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $rareBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";
my $buildTreeScr = getProgPaths("buildTree_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl";


my $arg5 = $ARGV[0];
my $bigF = $ARGV[4];
my $specificInD = $ARGV[1]; 
my $numTec = $ARGV[2];
my $tag = $ARGV[3];
my $outD = $ARGV[5];
#dir with gene catalog
my $GCd = $ARGV[6];#"$bigF/GeneCatalog/";
my $pid_relative = 0;
$pid_relative = $ARGV[7] if (@ARGV > 7) ;

my $inD = $bigF.$specificInD;

my $mapF = "";
$mapF = $ARGV[8] if (@ARGV > 8) ;
#die "$mapF";

#my $arg5 ="";
#$arg5 = $ARGV[6] if (@ARGV > 6);
print "Mode: $arg5\n";
my $ncore = 12;

my $mode = "x";#bat
my $getMOTUs = 0; my $contigMode=0; my $getAllGenes =0; my $FMGcrossSmpls = 0;
my $GCmod = 1;
if ($arg5 eq "singleSample" || $arg5 eq "contigsOnly"){
	if ($arg5 eq "contigsOnly"){$contigMode=1;}
	$GCmod=0;
	if ($mode eq "bat" ){$outD = "$inD/Binning/MetaBat/FMG/genes";
	} else { $outD = "$inD/Binning/MaxBin/FMG/genes"; }
} elsif ($arg5 eq "allGeneIDs"){
	$getAllGenes = 1;
} elsif ($arg5 eq "FMGcrossSmpls"){
	$getAllGenes = 1; $FMGcrossSmpls = 1;
} elsif ($arg5 eq "specificIDs"){
	$mode = $arg5;
}
system "mkdir -p $outD" unless (-d $outD);



#read GC indexes
my %gene2cl; my %cl2gene;
#if ( $arg5 eq "geneLnk" ){
if ($GCmod){
	my ($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx");
	%gene2cl = %{$hr1}; %cl2gene = %{$hr2};
}
my %COGid;my %allFMGs; 

#die $cl2gene{"1"}.$cl2gene{"429290"}."\n";



#------------------------------  specific IDs   -------------------------------------
if ($mode eq "specificIDs"){
	die "ARG1 needs to be a file with gene cat ids ($inD)\n" unless (-f $specificInD);
	print "Specific Gene ID from gene catalog run mode..\n";
	my @geneList;
	open I, "<$specificInD" or die "cant open $specificInD\n";
	while (<I>){chomp;my @spl = split(/\t/);push (@geneList,$spl[1]); $COGid{$spl[1]} = $spl[0];}	close I;
	#get abundance matrix
#	my $cmdX= "$rareBin lineExtr $GCd/Matrix.mat $outD/$tag.genes.matX $specificInD\n";
#	$cmdX.= "head -n 1 $GCd/Matrix.mat >  $outD/$tag.genes.mat;cat  $outD/$tag.genes.matX >>  $outD/$tag.genes.mat; rm  $outD/$tag.genes.matX;";
#	system $cmdX;
	
	#get sequences
	my ($hr1,$hr2,$hr3);
	($hr1) = extractGenesPerSample(\@geneList,"genes",0);
	exit(0) if ($arg5 eq "onlyIDs");
	my $outD2 = $outD."/rename$tag/"; system "rm -r $outD2; mkdir -p $outD2";
	#and extract the corresponding fna/ faa from every other dir.. return (\%COGgenes,\%nFilsFNA, \%nFilsFAA);
	($hr1,$hr2,$hr3) = extractFNAFAA2genes($hr1,$outD2,"gene");
	print "\nGene extraction finished\n";
	#	exit(0);
	
	my %nFilsFAA = %{$hr3}; my %nFilsFNA =%{$hr2}; my %COGgenes = %{$hr1};
	
	#my $str = join " ",keys(%nFilsFAA); system "cat " . $str . " > $outD2/allFAAs.faa";
	# $str = join " ",keys(%nFilsFNA); system "cat " . $str . " > $outD2/allFNAs.fna";
	open O,">$outD2/categories4ete.txt"; foreach (keys %COGgenes){print O $COGgenes{$_}."\n";}; close O;
	 
	#single gene trees
	my @fnaFs = sort (keys (%nFilsFNA)); my @faaFs = sort (keys (%nFilsFAA));
	system "rm -f  $outD2/allFAAs.faa $outD2/allFNAs.fna";
	foreach my $f (@fnaFs){system "cat $f >> $outD2/allFNAs.fna";}
	foreach my $f (@faaFs){system "cat $f >> $outD2/allFAAs.faa";}
	my $eteCmd = "";
		$ncore = 40;
	if ($numTec =~ m/S/){
		system "mkdir -p $outD/trees";
		for (my $i=0;$i < @fnaFs; $i++){
			#last;
			
			buildTree($fnaFs[$i],$faaFs[$i],"","$outD/trees/t_$nFilsFNA{$fnaFs[$i]}",0.99,$ncore);
			#$eteCmd = "ete3 build -n $fnaFs[$i] -a $faaFs[$i] -w clustalo_default-none-none-none --cpu $ncore -o $outD/tree$i --clearall --nt-switch 0.0 --noimg  --tools-dir /g/bork3/home/hildebra/bin/ete/ext_apps-latest"; #--no-seq-checks
			#die $eteCmd."\n";
			#print $i . " ";
			#system $eteCmd;
		}
	}
	if ($numTec =~ m/A/){
		system "mkdir -p $outD/tree";
		buildTree("$outD2/allFNAs.fna","$outD2/allFAAs.faa","$outD2/categories4ete.txt",$outD,0.85,$ncore);
	}
	#print "\n".$eteCmd."\n";
	
	exit(0);
}


#------------------------------  specific IDs  END  -------------------------------------




#load FMG markers in current smaple .. 
my $metaGD = `cat $inD/assemblies/metag/assembly.txt`; chomp $metaGD;
open I, "<$metaGD/ContigStats/FMG/FMGids.txt" or die "cant open $metaGD/ContigStats/FMG/FMGids.txt\n";
while (my $line = <I>){
	#MM1__C104459_L=563;_1 COG0552
	my @spl = split(/\s+/,$line);
	$COGid{$spl[0]} = $spl[1];
	$spl[0] =~ m/^(.*L=\d+;)(_\d+)/;
	if (exists($allFMGs{$1})){
		$allFMGs{$1} .= "_XX_".$spl[0];
	} else {
		$allFMGs{$1} = $spl[0];
	}
}
close I;

my $MBctgs = "$inD/Binning/MaxBin/MB.ctg.bin.txt";
my $startCnt=1; my $binMode = "MaxBin";
if ($mode eq "bat"){
	my $nms =  `cat $inD/Binning/MetaBat/MeBa.sto`; chomp $nms;
	$MBctgs = "$inD/Binning/MetaBat/$nms.fasta.fna";
	$binMode = "metaBat";
}


my @FMGs;my @noFMGs; my @contigs; my $ctgCnt =0;
#load binning file and search for FMG matches .. or just store the contig names
open I,"<$MBctgs" or die "cant open $MBctgs";
while (<I>){
	my @spl = split(/\s+/);
	if ($spl[1] eq $numTec){
		my $cC = $spl[0];
		
		if ($contigMode || $getAllGenes){push(@contigs,$cC); }
		
		$ctgCnt++;
		if (exists($allFMGs{$cC})){
			my @spl = split (/_XX_/, $allFMGs{$cC} );
			#print $spl[0]."\n";
			push(@FMGs,@spl);
		} else {
			push(@noFMGs,$cC);
		}
	}
}
close I;


#my @xx = keys %gene2cl; die $xx[0].'  '.$xx[0]."\n";

if ($getAllGenes){
	#scan cl2gene for matches to contig
	my @geneList; my $totGenFnd=0;
	foreach my $ct (@contigs){
		#die $ct."\n";
		my $geCnt = 1; my $testN = ">".$ct."_".$geCnt;
		print $testN." ";
		while (exists( $gene2cl{$testN} )){
			push(@geneList,$gene2cl{$testN});
			$geCnt++; $testN = ">".$ct."_".$geCnt; 
		}
		$totGenFnd += $geCnt;
		print "contig $ct has $geCnt genes\n";
	}
	print "Found $totGenFnd genes";
	my $genesA = "$outD/$tag.$binMode.gene.list";
	open O,">$genesA" or die "Can't open gene list file $genesA\n";
	foreach (@geneList){print O $_."\n";}
	close O;
	#and extract the corresponding gene profiles for FMGs
	my $cmdX= "$rareBin lineExtr -i $GCd/Matrix.mat -o $outD/$tag.$binMode.genes.matX -reference $genesA\n";
	$cmdX.= "head -n 1 $GCd/Matrix.mat >  $outD/$tag.$binMode.genes.mat;cat  $outD/$tag.$binMode.genes.matX >>  $outD/$tag.$binMode.genes.mat; rm  $outD/$tag.$binMode.genes.matX;";
	system $cmdX;
#	print $cmdX;
	print "Stored in $outD/$tag.$binMode.genes*\n";
	exit(0) unless($FMGcrossSmpls);
}

if ($contigMode){
	my $tarS = "$metaGD/scaffolds.fasta";
	my $srch = join('\' \'' ,@contigs);
	$srch =~ s/'>/'/g;
	#die $srch."\n";
	system "$samBin faidx $tarS '".$srch. "' > $outD/$tag.ctgs.fna" ;
	print "Wrote contigs to $outD/$tag.ctgs.fna\n";
	exit(0);
}


print "No FMGs in Ctg".join("  ",@noFMGs)."\n\n" if ( @noFMGs > 0 );
print "Found ".@FMGs." FMGs for $tag in $ctgCnt contigs, using mode $mode and Bin # $numTec\n";

#now start pulling these out of the GC..
my ($hr1,$hr2,$hr3);
($hr1) = extractGenesPerSample(\@FMGs,"FMG",1);
exit(0) if ($arg5 eq "onlyIDs");

#die();
my $outD2 = $outD."/rename$tag/"; system "rm -r $outD2; mkdir -p $outD2";
#and extract the corresponding fna/ faa from every other dir.. return (\%COGgenes,\%nFilsFNA, \%nFilsFAA);
($hr1,$hr2,$hr3) = extractFNAFAA2genes($hr1,$outD2,"FMG");
my %nFilsFAA = %{$hr3};my %nFilsFNA = %{$hr2};my %COGgenes = %{$hr1};

#do blast against FMG genes from gene cat to find relatives of target species
if ($pid_relative > 0){
	#define gene cat set of FMGs
	#$GCd/FMG/COG X.fna
}



#tree building
if ($buildTees){
	my $str = join " ",keys(%nFilsFAA); system "cat " . $str . " > $outD2/allFAAs.faa";
	 $str = join " ",keys(%nFilsFNA); system "cat " . $str . " > $outD2/allFNAs.fna";
	 open O,">$outD2/categories4ete.txt"; foreach (keys %COGgenes){print O $COGgenes{$_}."\n";}; close O;
	 
	system "mkdir -p $outD/tree";
	#-m cog_all-alg_concat_default-phyml_default    -m sptree_raxml_all   brh_cog_all-alg_concat_default-phyml_default
	buildTree("$outD2/allFNAs.fna","$outD2/allFAAs.faa","$outD2/categories4ete.txt",$outD,0.8,40);
	print " Done.\n$outD/tree\n";
} else {
	print "No tree is being built\n";
}

print "    ---FINISHED---\n$outD\n";
exit(0);



#######################################################
#######################################################

sub buildTree(){
	my ($fnFna, $aaFna,$cogCats,$outD,$ntFrac,$ncore) = @_;
	#my $ncore= 40;
	my $Ete = 0; #run ete or not?
	system "mkdir -p $outD" unless (-d $outD);
	$cogCats = "''" if ($cogCats eq "");
	my $cmd = "$buildTreeScr -fna $fnFna -aa $aaFna -cats $cogCats -outD $outD -cores $ncore -useEte $Ete -NTfilt $ntFrac > $outD/tree_build.log";
	#die $cmd;
	print "$cmd\n";
	if (system $cmd ) {"Failed $cmd\n";}
}



sub extractFNAFAA2genes(){
	my ($hr1,$outD2,$name) = @_;
	my %allIDs = %{$hr1};
	$mapF = $bigF."LOGandSUB/inmap.txt" if ($mapF eq "");
	my ($hrm,$hr2) = readMapS($mapF);
	my %map = %{$hrm};
	my @protFils; my @ntFils;
	#print keys $map{altNms}."\n";
	foreach my $sd ( keys %allIDs ){
		print "$sd - ";
		my $sd2 = $sd;
		unless (exists $allIDs{$sd}){print "skip\n";next;}
		if (exists(  $map{altNms}{$sd}  )){$sd2 = $map{altNms}{$sd}; print "Nx "}
		unless (exists ($map{$sd2}) ) {
			print "Can't find map entry for $sd\n"; next;
		}
		my $cD = $map{$sd2}{wrdir}."/";
		die "Assembly path missing $cD\n" unless (-e "$cD/assemblies/metag/assembly.txt");
		my $metaGD = `cat $cD/assemblies/metag/assembly.txt`; chomp $metaGD;
		my $tar = $metaGD."genePred/genes.shrtHD.fna";
		#print $tar."\n";
		my $srch = $allIDs{$sd}; $srch =~ s/'>/'/g;
		#print $srch."\n";
		print "Genes - ";
		system "$samBin faidx $tar ".$srch. " > $outD/$sd.$tag.$name.fna" ;
		#print $srch."\n";
		#die "$tar \n $outD/$sd.$numTec.$name.fna\n";
		renameFastHD("$outD/$sd.$tag.$name.fna",\%COGid,$sd);
		push(@ntFils,"$outD/$sd.$tag.$name.fna");
		
		print "Prot\n";
		$tar = $metaGD."genePred/proteins.shrtHD.faa";
		system "$samBin faidx $tar ".$srch. " > $outD/$sd.$tag.$name.prot.faa" ;	
		renameFastHD("$outD/$sd.$tag.$name.prot.faa",\%COGid,$sd);
		push(@protFils,"$outD/$sd.$tag.$name.prot.faa");
		
		#----------------   mOTU
		next unless ($getMOTUs);
		#mOTU - extract 100 bp +-
		my $gffHref= readGFF("$metaGD/genePred/genes.gff"); my %gff = %{$gffHref};
		my $fnaFllRef = readFasta("$metaGD/scaffolds.fasta.filt"); my %fna = %{$fnaFllRef};
		my @siIDs = split /' '/,$srch ;
		open OX,">$outD/$sd.$tag.$name.mOTU.ext.fna" ;
		foreach my $tID (@siIDs){
			$tID =~ s/'//g;
			#print $tID."\n".$gff{">".$tID}."\n";
			my $cGff = $gff{">".$tID};  my @spl = split /\t/,$cGff;
			my $curcog = $COGid{">".$tID};
			#if($spl[6] eq "-"){
			#reverse_complement_IUPAC 
			#print substr($fna{$spl[0]},$spl[3],($spl[4]-$spl[3]))."\n$curcog\n";
			if ($spl[4] >= (length($fna{$spl[0]}) - 100) ){die "too short $tID\n$cGff\n";}
			my $nSeq = substr($fna{$spl[0]},($spl[3]-100),($spl[4]-$spl[3]+200) );
			my $sta = 100; my $stop = length ($nSeq) - 100;
			print OX ">$curcog.$sta.$stop.$spl[6]\n$nSeq\n";
			#die;
		}
		#die $srch."\n";
		close OX;
		last;
	}
	

	#die;

	#rewrite to specific format (Jaimie) and ete
	my @samples = @{$map{opt}{smpl_order}}; 
	my %COGgenes; my %nFilsFNA ;  my %nFilsFAA;
	for (my $i=0;$i<@ntFils;$i++){
		my $hr = readFasta($ntFils[$i]);
		my %fas = %{$hr};
#		print $ntFils[$i]."\n";
		$hr = readFasta($protFils[$i]);
		my %faa = %{$hr};
		
		#check for double extractions and choose the longer gene
		my %lPerGene;
		my @keyss = sort keys %fas;
		foreach my $k (@keyss){
			my @spl = split(/_/,$k);
			if ($spl[1] =~ m/(.*)\.\.\d+$/){
				my $curL = length($fas{$k});
				#replace old (shorter gene)
				die $k." @keyss\n" unless (exists $lPerGene{$1});
				#print " $1 $lPerGene{$1}  ..  $curL  $k   @keyss\n";
				if ($lPerGene{$1} < $curL ){  $lPerGene{$1} = $curL; $k =~ m/(.*)\.\.\d+$/; $fas{$1} = $fas{$k}; $faa{$1} = $faa{$k};}
				#keep longer gene
				delete $fas{$k};delete $faa{$k};
			} else {
				#print " X $spl[1] ";
				$lPerGene{$spl[1]} = length($fas{$k});
			}
		}
		foreach my $k (keys %fas){
			my @spl = split(/_/,$k);
			if (!exists($COGgenes{$spl[1]})){$COGgenes{$spl[1]}=$k;}else{$COGgenes{$spl[1]}.= "\t".$k;  }
			open O,">>$outD2/$spl[1].fna";		print O ">".$k."\n".$fas{$k}."\n";		close O;
			$nFilsFNA{"$outD2/$spl[1].fna"} = $spl[1];
		}
		foreach my $k (keys %faa){
			my @spl = split(/_/,$k);
			open O,">>$outD2/$spl[1].faa";		print O ">".$k."\n".$faa{$k}."\n";		close O;
			$nFilsFAA{"$outD2/$spl[1].faa"} = $spl[1];
		}
	}
	return (\%COGgenes,\%nFilsFNA, \%nFilsFAA);
}

#just gets the gene name in each sample, based on gene ID
sub extractGenesPerSample($ $ $ $){
	my ($FMGsAR,$name, $trackGeneNumber,$minNumG) = @_;
	my @FMGs = @{$FMGsAR};
	#my %COGid = %{$COGidHR};
	my %allIDs; my %bigTracker; my %double;
	my $namesF = $tag."Names.txt";
	my $namAssoF = $tag."Genes.txt"; 
	open OR,">$outD/$namesF" or die "Can;t open $name names out file\n";
	open ORa,">$outD/$namAssoF" or die "Can;t open $name Genes out file\n";
	foreach my $v (@FMGs){
		#this parts translates gene name to gene number
		my @ssp;
		if ($trackGeneNumber){ #get gene number first
			@ssp = ($v,">".$v);
			if ($GCmod){ #get from gene catalog homologs in different Samples
				#@ssp = split(/\s+/,`grep '$v' $GCd/compl.incompl.95.fna.clstr.idx`);
				die "can't find gene $v in gene cat\n" unless (exists($gene2cl{">".$v}));
				my $genNum = $gene2cl{">".$v};
				@ssp = ($genNum,$cl2gene{$genNum});
				#die "@ssp\n$v\n";
			}
		} else { #gene number already given in @FMGs
			#print $v."\n".$cl2gene{"429290"}."\n$cl2gene{$v}\n";

			die "can't find gene number $v\n" unless (exists($cl2gene{$v}) );
			@ssp = ($v,$cl2gene{$v});
		}
		print OR $ssp[0]."\n";
		my @tmpArr = split /,/,$ssp[1];
		print ORa $ssp[0]."\t".@tmpArr."\t".join(",",sort(@tmpArr))."\n";
		my @allC = split(/,/, $ssp[1]);
		#and this gene number info info is reformated into genes per sample format
		foreach my $t (@allC){ #go through to find sample ID
			chomp $t;
			die "can't find $name $v in COGid hash\n" unless (exists $COGid{$v});
			$COGid{$t} = $COGid{$v};
			#$t =~ s/^>//;
			#print $t."\t$COGid{$t}\n";
			$t =~ m/>(.*)__/;
			
			#print $1."\n";
			if ( exists($allIDs{ $1 }) ){
				$allIDs{ $1 } .= " '".$t."'";
			} else {
				$allIDs{ $1 } = "'".$t."'";
			}
			if (exists ($bigTracker{$1} {$COGid {$t} }) ){
				#print $1."\n";
				if (exists ($double{ $bigTracker{$1} {$COGid {$t} } } ) ){
					$double{ $bigTracker{$1} {$COGid {$t} } } .= " ".$t;
				} else {
					$double{ $bigTracker{$1} {$COGid {$t} } } = $t;
					
				}
			} else {
				$bigTracker{$1} {$COGid {$t} } = $t;
			}
		}
		#last;
		#rename to T$numTec_COG_MMX
	}
	close OR;close ORa; #namAssoF 
	my $cmdX= "$rareBin lineExtr -i $GCd/Matrix.mat -o $outD/$tag.matX -reference $outD/$namesF\n";
	$cmdX.= "head -n 1 $GCd/Matrix.mat > $outD/$tag.mat;cat $outD/$tag.matX >> $outD/$tag.mat; rm $outD/$tag.matX;";
	system $cmdX;
	
	print "Detected homologs in GC\n";
	my @kDbl = keys %double;
	if (@kDbl > 1){
		print "Found ".@kDbl." samples with >1 $name\n";
	}
	return (\%allIDs);
}





sub readMap2{
	my ($inF) = @_;
	my $cnt =-1; my %ret; my @order = ();
	open I,"<$inF" or die "Can't open map $inF\n";
	while (<I>){
		chomp;	$cnt++;
		if ($cnt == 0|| m/^#/){next;} #print"next"; maybe later check for col labels etc
		my @spl = split(/\t/);
		#die $spl[0]." ".$spl[1]."\n";
		$ret{$spl[0]}{dir} = $spl[1];
		$ret{$spl[0]}{SmplID} = $spl[0];
		push(@order,$spl[0]);
	}
	close I;
	$ret{opt}{smpl_order} = \@order;
	return \%ret;
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
