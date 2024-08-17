package Mods::phyloTools;
use warnings;
#use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 

#use Mods::GenoMetaAss qw(qsubSystem);

use Exporter qw(import);
our @EXPORT_OK = qw(runRaxMLng runRaxML readFMGdir prep40MGgenomes prepNOGSETgenomes getE100 getGenoGenes getFMG renameFMGs 
			runFasttree runQItree fixHDs4Phylo getGenoName calcDisPos2 getTreeLeafs filterMSA MSA);
use Mods::GenoMetaAss qw(systemW readFasta renameFastHD gzipwrite gzipopen);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::FuncTools qw(assignFuncPerGene readGene2Func);


sub zorroFilter{
	my ($inMSA) = @_;
	my ($HR,$OK) = readFasta($inMSA);
	my %FA = %{$HR};
	my @FAhds = keys %FA;

	my $zorroBin = getProgPaths("zorro");#"/g/bork3/home/hildebra/bin/zorro/./zorro_linux_x86_64";#"/g/bork3/home/hildebra/bin/zorro/bin/./zorro";
	my $fastTBin = getProgPaths("fasttree");

	my $cmd = "$zorroBin -treeprog $fastTBin -Nsample ". scalar(@FAhds)*2 . " $inMSA ";
	#die $cmd."\n";
	my $tSco = `$cmd`;
	my @sco = split /\n/,$tSco;
	my $GoodPos=0;my $BadPos=0;
	for (my $i=0;$i<@sco;$i++){
		if ($sco[$i] < 6){
			$BadPos++;
			foreach my $k (@FAhds){
				my $cur = substr($FA{$k},$i,1);
				substr($FA{$k},$i,1) = "N" if ($cur ne "-");
			}
		} else {
			$GoodPos++;
		}
	}
	if ($BadPos > 0 ){print "Zorro Filter: $BadPos vs $GoodPos;; ";}
	writeFasta(\%FA,"$inMSA");
	return ($BadPos / ($GoodPos+$BadPos ));
	#die "$inMSA\n";
}

#builds an MSA, including filtering etc
sub MSA{
	my ($tmpInMSA,$tmpOutMSA2,$ncore,$clustalUse,$continue,$numSeq)  = @_;
	my $cmd = "";
	if ($clustalUse==1){
		my $clustaloBin = getProgPaths("clustalo");#= "/g/bork3/home/hildebra/bin/clustalo/clustalo-1.2.0-Ubuntu-x86_64";
		$cmd = $clustaloBin." -i $tmpInMSA -o $tmpOutMSA2 --outfmt=fasta --threads=$ncore --force\n"; #--guidetree-out $tmpTree 
	} elsif ($clustalUse == 0) {
		my $msapBin = getProgPaths("msaprobs");
		$cmd = "$msapBin -num_threads $ncore $tmpInMSA > $tmpOutMSA2";
	} elsif ($clustalUse == 2) {
		my $mafftBin = getProgPaths("mafft");#= "/g/bork3/home/hildebra/bin/clustalo/clustalo-1.2.0-Ubuntu-x86_64";
		$cmd = "$mafftBin --thread $ncore --quiet $tmpInMSA > $tmpOutMSA2";
	} elsif ($clustalUse == 3) {
		die "guidance: rework phyloTools.pm\n";
		my $guid2Path = "/g/bork3/home/hildebra/bin/guidance.v2.02/";
		#guidance has some strange results..
		$cmd = " $guid2Path/www/Guidance/guidance.pl --seqFile $tmpInMSA --msaProgram MAFFT --seqType aa ";
		#$cmd .= " --dataset $spl2[1].$cnt --mafft $mafftBin --outDir $tmpD --proc_num $ncore\n";
	} elsif ($clustalUse == 4) {
		#my $nseqs = 0;
		#`grep -c '^>' $tmpInMSA`;
		my $algo = "-align";
		$algo = "-super5" if ($numSeq > 300); #just too slow otherwise..
		my $MUSCLE5Bin = getProgPaths("MUSCLE5");
		$cmd = "$MUSCLE5Bin $algo $tmpInMSA -output $tmpOutMSA2 -threads $ncore "; #

	} else {
		die "msa option ($clustalUse) unknown)\n";
	}
	#die $cmd;
	systemW $cmd unless (-s $tmpOutMSA2 && $continue);
	return $tmpOutMSA2;
}

sub filterMSA{ #pretty useless atm.. not really used
	my ($tmpInMSA,$tmpOutMSA2,$ncore,$postFilter,$useAA4tree)  = @_;
	my $zScore = 0;
	my $cmd = "";
	if ($postFilter eq "macse"){ #gives strange results..
		#-out_NT output_NT.fasta -out_AA output_AA.fasta
		my $macseBin = "java -jar /g/bork3/home/hildebra/bin/macse_v2.03.jar";
		my $outTag = "-out_AA"; $outTag = "-out_NT" if ($useAA4tree);
		$cmd = "$macseBin -prog refineAlignment -align $tmpOutMSA2 $outTag $tmpOutMSA2.2";
	} elsif ($postFilter eq "zorro"){
		$zScore = zorroFilter($tmpOutMSA2);
	}
	if ($zScore > 0.4){
		$tmpOutMSA2 = "";
	}

	return $tmpOutMSA2;

}


sub getTreeLeafs($){
	my ($nwkFile)  = @_;
	my $nwk = `cat $nwkFile`;
	my %nwLfs;
	my @nwLfs1 = split /[\(\),;]+/,$nwk;
	for (my $i=0;$i<@nwLfs1;$i++){
		my @tmp = split /:/,$nwLfs1[$i];
		next if (@tmp ==0 || $tmp[0] eq "" );
		$nwLfs{$tmp[0]} = 1;
	}
	return(\%nwLfs);
}
sub fixHDs4Phylo ($){
	#routine to check that headers of fastas don't contain ":", ",", ")", "(", ";", "]", "[", "'"
	my ($inF) = @_;
	my $reqFix=0;
	if ($inF eq ""){return "";}
	my $hr = readFasta($inF,1); my %FAA = %{$hr};
	foreach my $hd (keys %FAA){
		if ($hd =~ m//){
			$reqFix=1;last;
		}
	}
	my $outF = $inF;
	if ($reqFix){
		$outF .= ".fix";
		print "Fixing headers in input file (to $outF)\n";
		my %newHDs;
		open O,">$outF";
		foreach my $hd (keys %FAA){
			my $hd2 = substr $hd,0,40; #cut to raxml length
			$hd2 =~ s/[:,\}\{;\]\[']/|/g;
			$newHDs{$hd2} ++;
			$hd2 .= $newHDs{$hd2};
			print O ">$hd2\n$FAA{$hd}\n";
		}
	}
	return $outF;
}

sub runQItree{
	my ($hr) = @_; my %treeOpts = %{$hr};
	my ($inMSA,$treeOut,$ncore,$outgr,$bootStrap,$useAA,$fast,$autoModel,$partiF) = 
		($treeOpts{inMSA},$treeOpts{IQtreeout},$treeOpts{ncore},$treeOpts{outgr},$treeOpts{bootStrap},$treeOpts{useAA},
		$treeOpts{iqtreeFast},$treeOpts{autoModel},$treeOpts{partition});
	
	#die "AA use $useAA\n";
	my $constraintTree = $treeOpts{constraintTree};
	die ("Constraint tree $constraintTree does not exist") if ($constraintTree ne "" && !-e $constraintTree);
	my $iqTree  = getProgPaths("iqtree");
	my $vcheck = `$iqTree --version`;
	unless ($vcheck =~ m/version 2/){die "Needs iqtree version 2\n:$vcheck\n";}
	$treeOut =~ s/\.nwk$//;
	my $treNM = "IQtree";
	my $cmd = "$iqTree -s $inMSA -T $ncore -pre $treeOut -seed 678 "; #-nt AUTO -ntmax $ncore
	$cmd .= " -p $partiF " unless ($partiF eq "");
	$cmd .= "-o $outgr " unless ($outgr eq "" && $outgr !~ m/,/);
	$cmd .= "-g $constraintTree " unless ($constraintTree eq "");
	$cmd .= "--quiet " if (exists($treeOpts{silent}) && $treeOpts{silent});
	unless ($fast == 0){$cmd .= "--fast "; print "IQtree - fast\n"; $treNM .= "_fast";}
	if ($autoModel){$treNM .= "_autoMOD";}
	if ($useAA){
		if ($autoModel){
			$cmd .= "-m TEST  "; 
		} else{
			$cmd .= "-m LG+F+G "; #needs to be HKY for nts
		}
	} else {
		if ($autoModel){
			$cmd .= "-m TEST ";#-mset HKY,HKY+F,HKY+F+I,HKY+F+I+G4,JC,F81,K2P,K3P,K81uf,GTR "; 
		} else {
			$cmd .= "-m GTR+F+I+G4 "; #default model, as spotted on 40 MG for phylo tree..
			#$cmd .= "-m HKY+F+G "; 
		}
	}
	if ($bootStrap >0){
		if ($bootStrap < 1000){
			print "standard non parametric bootstrap ($bootStrap). Use >1000 bootstraps to do ultrafast bootstrap\n";
			die  "normal bootstrap requires >= 100 iterations ($bootStrap given)\n" if ($bootStrap <100);
			$cmd .= "-b $bootStrap " ;
		} elsif ($fast == 0) {
			print " ultrafast bootstrap $bootStrap\n";
			$cmd .= "-B $bootStrap " ;
			#die "$fast\n";
		} else {
			print " SH-like approximate likelihood ratio test $bootStrap\n";
			$cmd .= "-alrt $bootStrap " ;
		}
		#also consider -b >=100 for std bootstrap
	} else {
		$cmd .= "--alrt 1000 ";
	}
	#TODO: include booster for better bootstrap values
	#booster -a tbe -i 40MG.IQtree.treefile -b 40MG.IQtree.boottrees -@ 8 -o 40MG.IQtree_booster.tre
	$cmd .= "";
	#die $cmd ;#if ($constraintTree ne "");
	systemW $cmd;
	#$treNM .= ".nwk";
	#"mv $treeOut/IQtree_fast_allsites.treefile $treeOut/$treNM";
}
sub runFasttree{
	my ($inMSA,$treeOut,$isAA,$ncore) = @_;
	my $fsttreeBin  = getProgPaths("fasttree");
	my $ntFlag = "";
	$ntFlag = "-nt -gtr" if (!$isAA);
	my $cmd = "$fsttreeBin $ntFlag $inMSA > $treeOut\n";
	systemW $cmd;
}


sub getGenoName($){my $GenomeN = $_[0];$GenomeN =~ s/\.f[n]?a$//; $GenomeN =~ s/^.*\///;  $GenomeN =~ s/-/_/; return $GenomeN;}

#routine to format marker genes with naming etc to work with buildTree script
sub prep40MGgenomes{ 
	#\@refGenos,$rewrite,$ncore,$fnFna1,$aaFna1,$cogCats1)
	my @refGenos = @{$_[0]};my $finalD=$_[1]; my $tag= $_[2];
	my $rewrite = $_[3]; my $ncore=$_[4];
	my ($fnFna1,$aaFna1,$cogCats1) = ("","","");
	my $outGrp = "";
	$outGrp = $_[5] if (@_ > 5);
	if (@_ > 6){$fnFna1= $_[6];$aaFna1= $_[7];$cogCats1= $_[8];}
	print "Using outgroup $outGrp\n" if ($outGrp ne "");
	
	my $buildTreeScr = getProgPaths("buildTree_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl";
	my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
	my $FMGrwkScr = getProgPaths("FMGrwk_scr");

	#system "$renameCtgScr $refGenos $GenomeN";
	system "mkdir -p $finalD" unless (-d $finalD);
	my $fnFna = $fnFna1; if ($fnFna eq ""){$fnFna = $finalD."$tag.fna";} else {$fnFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $aaFna = $aaFna1; if ($aaFna eq ""){$aaFna = $finalD."$tag.faa";} else {$aaFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $cogCats = $cogCats1; if ($cogCats eq ""){$cogCats = $finalD."$tag.cats.txt";} else {$cogCats =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	if (-e $fnFna1){system("cat $fnFna1 > $fnFna");} else {system "rm -f $fnFna";}#shortFNAhd($fnFna);}
	if (-e $aaFna1){system("cat $aaFna1 > $aaFna");} else {system "rm -f $aaFna";}
	my %COGgenes; my %allGenomes;
	#assumes ref genome, no genes called
	for (my $i=0;$i<@refGenos;$i++){
		my $refG = $refGenos[$i];
		next if ($refG eq "");
		my $GenomeN = getGenoName($refGenos[$i]);
		
		if (exists($allGenomes{$GenomeN})){
			print "Genome $GenomeN already exists; assumming double entry and skipping genome\n";next;
		} else {
			$allGenomes{$GenomeN} = 1;
		}
		#die "$GenomeN\n";
		my ($ntGenes,$proteins) = getGenoGenes($refG);
		#my $proteins = $refG; my $ntGenes = $refG; $proteins =~ s/\.[^\.]+$/\.genes\.faa/;
		#$ntGenes =~ s/\.[^\.]+$/\.genes\.fna/;
		#die "Can't find ref genome $refG\n" unless (-e $refG || (-e $ntGenes && -e $proteins));
		#my $prodigal_cmd .= "$prodigalBin -i $refG -a $proteins -d $ntGenes -f gff -p single > /dev/null\n";
		#die $prodigal_cmd."\n$rewrite || !-e $proteins || !-e $ntGenes\n";
		#system $prodigal_cmd if ($rewrite || !-e $proteins || !-e $ntGenes);
		#fetch mg
		my $cmd ;
		my $GenomeDir = $proteins; $GenomeDir =~ s/[^\/]+$//;
		my $outDFMG = $GenomeDir."FMGs$GenomeN/";
		#getFMG($outDFMG,$ntGenes,$proteins,$ncore,$rewrite+3);
		if (!-s "$outDFMG/COG0012.faa" && !-s "$outDFMG/COG0016.faa"){
			$cmd = "perl $FMGd/fetchMG.pl -m extraction -o $outDFMG -l $FMGd/lib -t $ncore -d $ntGenes $proteins"; # -x $FMGd/bin  -b $FMGd/lib/MG_BitScoreCutoffs.uncalibrated.txt
			system $cmd;
			system "rm  -rf $GenomeDir/*.cidx $outDFMG/temp $outDFMG/hmmResults";
		}
		my $hr = renameFMGs($outDFMG,$GenomeN,"allFMG",1);
		#system "$FMGrwkScr $outDFMG" ;
		my %COGg = %{$hr};
		foreach my $k (keys %COGg){
			push (@{$COGgenes{$k}},$COGg{$k});
		}
		#open I,"<$outDFMG/FMGids.txt";
		#while (<I>){
		#	chomp;my @spl = split /\s/;
		#	if (!exists($COGgenes{$spl[1]})){$COGgenes{$spl[1]}=$spl[0];}else{$COGgenes{$spl[1]}.= "\t".$spl[0];  }
		#}
		#close I;
		#now add this info to the existing ete input files
		system("cat $outDFMG/*.fna >> $fnFna");#shortFNAhd($fnFna);}
		system("cat $outDFMG/*.faa >> $aaFna") ;#shortFNAhd($aaFna);
	}
	#open genecat file 
	my $cogcatstr = "";my $catCnt=0;
	if (-e $cogCats1){
		open I,"<$cogCats1";
		while(<I>){
			chomp;my @spl = split /\t/;
			$spl[0] =~ m/_([^_]+)$/;
			die "Can't find category $1 in ref gene set\n" unless (exists($COGgenes{$1}));
			push (@spl,@{$COGgenes{$1}});
			$cogcatstr .= join ("\t",@spl) . "\n";
			$catCnt++;
		}
		close I;
	} else { #make anew
		foreach my $k (keys %COGgenes){
			$cogcatstr .= join ("\t",@{$COGgenes{$k}}) . "\n";
			$catCnt++;
		}
	}
	open O,">$cogCats";
	print O $cogcatstr;
	close O;
	print "Found $catCnt / ".scalar(keys %COGgenes)." gene cats\n";
	#and run tree building
	my $ogrpStr="";
	if ($outGrp ne ""){
		$ogrpStr = "-outgroup $outGrp";
	}
	my $cmd = "$buildTreeScr -fna $fnFna -aa $aaFna -cats $cogCats -outD $finalD -cores $ncore -NTfilt 0.15 -bootstrap 100  $ogrpStr > $finalD/tree_build.log";
	return $cmd;
}
#similar to prep40MGgenomes, but extended to more genes (as given by NOG annotation)
sub prepNOGSETgenomes{ 
	#\@refGenos,$rewrite,$ncore,$fnFna1,$aaFna1,$cogCats1)
	my @refGenos = @{$_[0]};my @NOGs = @{$_[1]};
	my $finalD=$_[2]; my $tag= $_[3];
	my $rewrite = $_[4]; my $ncore=$_[5];
	my ($fnFna1,$aaFna1,$cogCats1) = ("","","");
	my $outGrp = "";
	$outGrp = $_[6] if (@_ > 6);
	#these files will be added on top of the newly derrived genes
	if (@_ > 7){$fnFna1= $_[7];$aaFna1= $_[8];$cogCats1= $_[8];}
	print "Collating ".@NOGs." NOG genes for tree building\n";
	print "Using outgroup $outGrp\n" if ($outGrp ne "");
	
	my $buildTreeScr = getProgPaths("buildTree_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/buildTree.pl";
	my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
	my $FMGrwkScr = getProgPaths("FMGrwk_scr");

	#system "$renameCtgScr $refGenos $GenomeN";
	system "mkdir -p $finalD" unless (-d $finalD);
	my $tmpP = $finalD."tmp/";
	system "mkdir -p $tmpP" unless (-d $tmpP);
	
	my $fnFna = $fnFna1; if ($fnFna eq ""){$fnFna = $finalD."$tag.fna";} else {$fnFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $aaFna = $aaFna1; if ($aaFna eq ""){$aaFna = $finalD."$tag.faa";} else {$aaFna =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	my $cogCats = $cogCats1; if ($cogCats eq ""){$cogCats = $finalD."$tag.cats.txt";} else {$cogCats =~ s/\.([^\.]+)$/\.$tag\.$1/;}
	if ( -e $fnFna1){system("cat $fnFna1 > $fnFna");} else {system "rm -f $fnFna";}#shortFNAhd($fnFna);}
	if ( -e $aaFna1){system("cat $aaFna1 > $aaFna");} else {system "rm -f $aaFna";}

	my %COGgenes; my %allGenomes;
	

	#assumes ref genome, no genes called
	for (my $i=0;$i<@refGenos;$i++){
		my $refG = $refGenos[$i];
		next if ($refG eq "");
		my $GenomeN = getGenoName($refGenos[$i]);
		
		if (exists($allGenomes{$GenomeN})){
			print "Genome $GenomeN already exists; assumming double entry and skipping genome\n";next;
		} else {
			$allGenomes{$GenomeN} = 1;
		}
		#die "$GenomeN\n";
		my ($ntGenes,$proteins) = getGenoGenes($refG);
		#die "$proteins\n";
		#functional annotation via NOG
		my $GenomeDir = $proteins; $GenomeDir =~ s/[^\/]+$//;
		my $outDNOG = $GenomeDir."NOGs$GenomeN/";
		my %optsDia = (eval=>1e-12,percID=>35,minPercSbjCov=>0.6,fastaSplits=>1,ncore=>$ncore,splitPath=>$tmpP,bacNOG=>1);
		my ($outF,$jdep) = assignFuncPerGene($proteins,$outDNOG,$tmpP,"NOG",\%optsDia);
		my $tarAnno = "${outF}geneAss.gz";
		#print "XX $tarAnno\n";
		
		my %COG2Gn;
		my ($I,$OK) = gzipopen($tarAnno,"Func anno",1); 
		while (<$I>){
			chomp;
			my @spl = split /\t/;
			$COG2Gn{$spl[1]} = $spl[0];
			#print;
		}
		close $I;
		#print "$tarAnno\n";
		#next; #DEBUG
		#now load AA and NT file, and rename genes to fit
		my $hr = readFasta($ntGenes,1);
		my %FAS = %{$hr};
		$hr = readFasta($proteins,1);
		my %FASAA = %{$hr};
		#print "Read NT AA\n";
		my %COGg ;
		open O,">>$fnFna" or die "Can't open $fnFna\n";
		open O2,">>$aaFna" or die "Can't open $aaFna\n";
		my $hitcnt=0;
		my $cnt=0;
		foreach my $tag (@NOGs){
			my $newName = "${GenomeN}_${tag}";
			#print $tag."\n" if ($cnt<5);; $cnt++;
			#if ($tag =~ m//)
			if (exists($COG2Gn{$tag})){
				#print "$COG2Gn{$tag}\n";
				die "Can't find $COG2Gn{$tag} in nt's\n" unless(exists $FAS{$COG2Gn{$tag}});
				die "Can't find $COG2Gn{$tag} in AA's\n" unless(exists $FASAA{$COG2Gn{$tag}});
				print O ">$newName\n$FAS{$COG2Gn{$tag}}\n";
				print O2 ">$newName\n$FASAA{$COG2Gn{$tag}}\n";
				$COGg{$tag} = $newName;
				$hitcnt++;
			}
		}
		close O;
		close O2;
		print "Found $hitcnt core genes in genome $GenomeN\n";
		
		foreach my $k (keys %COGg){
			push (@{$COGgenes{$k}},$COGg{$k});
		}
		#die "$fnFna\n";
		#last if ($i==10);
	}

	#open genecat file 
	my $cogcatstr = "";my $catCnt=0;
	if (-e $cogCats1){
		open I,"<$cogCats1" or die "Can't open $cogCats1\n";
		
		while(<I>){
			chomp;my @spl = split /\t/;
			$spl[0] =~ m/_([^_]+)$/;
			die "Can't find category $1 in ref gene set\n" unless (exists($COGgenes{$1}));
			push (@spl,@{$COGgenes{$1}});
			$cogcatstr .= join ("\t",@spl) . "\n";
			$catCnt++;
		}
		close I;
	} else { #make anew
		foreach my $k (keys %COGgenes){
			$cogcatstr .= join ("\t",@{$COGgenes{$k}}) . "\n";
			$catCnt++;
		}
	}
	open O,">$cogCats";
	print O $cogcatstr;
	close O;
	print "Found $catCnt / ".scalar(keys %COGgenes)." gene cats\n";
	#and run tree building
	#die "$cogCats\n";
	my $ogrpStr="";
	if ($outGrp ne ""){
		$ogrpStr = "-outgroup $outGrp";
	}
	my $cmd = "$buildTreeScr -fna $fnFna -aa $aaFna -cats $cogCats -outD $finalD -cores $ncore -NTfilt 0.15 -bootstrap 100  $ogrpStr > $finalD/tree_build.log";
	system "rm -r $tmpP" ;
	return $cmd;
}


sub runRaxMLng{
	my ($hr) = @_; my %treeOpts = %{$hr};
	my ($mAli,$outTree,$ncore,$outgroup,$bootStrap,$useAA,$autoModel,$cont) = 
		($treeOpts{inMSA},$treeOpts{RAXNGtreeout},$treeOpts{ncore},$treeOpts{outgr},$treeOpts{bootStrap},$treeOpts{useAA},
		$treeOpts{autoModel}, $treeOpts{cont});
	my $cmd = "";
	my ($raTmpF,$raTmpF2,$raTmpF3) = ("RXMtmp","R2XM","RX3M"); #my $raxFile = "RXMall";
#	my $raxFile2 = "RXMsyn";
	my $testModel = getProgPaths("modeltest");
	my $raxmlBin = getProgPaths("raxmlng");
	my $trDist = getProgPaths("treeDistScr");
	my $model = "GTR+G4+I";
	if ($useAA){
		$model = "LG+I+G";
	}
	my $outTreeBasic = $outTree;
	$outTreeBasic =~ s/\.[^\.]+$//;
	
	#big modeltest section..
	my $ntaaS = "-d nt"; $ntaaS = "-d aa" if ($useAA);
	if ($autoModel){
		$cmd = "$testModel -p $ncore -i $mAli -t ml -T raxml $ntaaS > $outTreeBasic.MT\n";
		print "starting model test..\n";
		systemW "$cmd\n";
		my $tmp = `cat $outTreeBasic.MT`;
		$tmp =~ m/> raxml-ng --msa.*--model (\S+)/;
		$model = "$1";
		#unlink("$outTreeBasic.MT");
		print "Best model is $model\n";
	}
	
	print "using model $model ($useAA)\n";
	
	my $contS = ""; $contS = "--redo " unless ($cont);
	my $outgrS = ""; $outgrS = "--outgroup $outgroup" if ($outgroup ne "");
	$cmd = "$raxmlBin --msa $mAli --force --model $model --prefix $outTree $contS $outgrS --threads $ncore --site-repeats on \n"; #--seed 52352
	$cmd .= "mv $outTree.raxml.bestTree $outTree; mv $outTree.raxml.log $outTreeBasic.log\nrm $outTree.raxml*";
	systemW "$cmd\n";
	return $cmd;
}

#de novo aligns pairwise via vsearch and calcs id (iddef 2)
sub calcDisPos2{
	my ($MSA,$opID, $isNT) = @_;
	my $ncore = 1;	$ncore = $_[3] if (@_ > 3);
	my $tmpD = "$MSA"; $tmpD =~ m/^(.*)\/[^\/]+$/; $tmpD = $1."/";
	$tmpD = $_[4] if (@_ > 4);
	my $doWr = 1;
	$doWr = 0 if ($opID eq "");
	die "can't find input file $MSA\n" unless (-e $MSA);
	system "mkdir -p $tmpD" unless (-d $tmpD);
	my $vsBin  = getProgPaths("vsearch");
	#my $swBin = getProgPaths("swipe");
	#my $mkBldbBin = getProgPaths("makeblastdb");
	my $diaBin = getProgPaths("diamond");
	my $cmd="";
	$cmd = "$vsBin --threads $ncore --quiet --allpairs_global $MSA --acceptall --blast6out $tmpD/tmp.b6 --iddef 2 \n";
	my $prog = 0; #NT
	if ($isNT == 0){
		$prog = 1 ; #AA
#		$cmd = "$mkBldbBin -in $MSA -dbtype 'prot'\n";
#		$cmd .= "$swBin -p $prog -a $ncore -q $MSA  -d $MSA -c 0 -e 999 -m 0 -o $tmpD/tmp.b6  -b 9999999999999 \n" ;
		$cmd = "$diaBin makedb --in $MSA -d $MSA.db -p $ncore --quiet\n";
		$cmd .= "$diaBin blastp -f 6 -k 99999999 --min-score 0 --unal 1 --masking 0 --compress 0 --quiet -d $MSA.db -q $MSA -e 1e-4 --sensitive -o $tmpD/tmp.b6 -p $ncore\n";


		#die "$cmd\n";
	}
	#$cmd = $clustaloBin." -i $MSA -o $MSA.tmp --outfmt=fasta --percent-id --use-kimura --distmat-out $opID --threads=$ncore --force --full\n";
	#$cmd .= "rm -f $MSA.tmp\n";
	#die $cmd."\n";
	systemW $cmd;
	#die $cmd;
	my %perID; my %perIDpSmpl; my %perIDpSmplSum;
	open I,"<$tmpD/tmp.b6" or die "Can't find $tmpD/tmp.b6\n";
	while (<I>){
		chomp;
		my @spl = split /\t/;
		#if (exists($perID{$spl[0]}{$spl[1]})){die "Entry in dis mat already exists: $spl[0] $spl[1]\n$cmd\n$tmpD/tmp.b6\n";}
		$perID{$spl[0]}{$spl[1]} = $spl[2];
		$perID{$spl[1]}{$spl[0]} = $spl[2];
	}
	close I;
	if ($doWr){
		open O,">$opID" or die "Can't open dist mat output $opID\n";
	}
	my @smpls = sort keys %perID;
	print O "percID\t".join("\t",@smpls)."" if ($doWr);
	my $totalPID = 0;my $totalPIDcnt=0;
	foreach my $k1 (@smpls){
		print O "\n$k1\t" if ($doWr);
		foreach my $k2 (@smpls){
			if ($k1 eq $k2){print O "\t100" if ($doWr);
			} else {
				if (!exists($perID{$k1}{$k2})){print O "\t0" if ($doWr);}#die "Can;t locate keys $k1 $k2\n";}
				else {print O "\t".$perID{$k1}{$k2} if ($doWr); 
					$totalPID+= $perID{$k1}{$k2}; $totalPIDcnt++;
					$perIDpSmpl{$k1} += $perID{$k1}{$k2}; $perIDpSmplSum{$k1}++;
				}
			}
		}
		#DEBUG
		
		#if ($k1 =~ m/PNPM68/){my $ovA=0;foreach my $k2 (@smpls){$ovA += $perID{$k1}{$k2} if (exists($perID{$k1}{$k2}));} print "PNPM68 ".$ovA/scalar(@smpls)."\n";}
	}
	close O if ($doWr);
	foreach my $k (keys %perIDpSmpl){
		$perIDpSmpl{$k} /= $perIDpSmplSum{$k};
	}
	#die "$tmpD/tmp.b6\n";
	unlink "$tmpD/tmp.b6";
	unlink "$MSA.db.dmnd" if (!$isNT);
	my $avgID = 0;
	$avgID = ($totalPID/$totalPIDcnt) if ($totalPIDcnt > 0);
	return $avgID,\%perIDpSmpl,\%perID;
}


sub runRaxML{
	my ($hr) = @_; my %treeOpts = %{$hr};
	my ($mAli,$outTree,$ncore,$outgroup,$bootStrap,$useAA,$autoModel,$cont) = 
		($treeOpts{inMSA},$treeOpts{RAXtreeout},$treeOpts{ncore},$treeOpts{outgr},$treeOpts{bootStrap},$treeOpts{useAA},
		$treeOpts{autoModel}, $treeOpts{cont});
	#my ($mAli,$bootStrap,$outgroup,$outTree,$ncore) = @_;
	#my $cont = 0; my $useNT=1; my $useAA=0;
	#$cont = $_[5] if (@_ > 5);
	#$useNT = $_[6] if (@_ > 6);
	#die "$cont\n";
	my ($raTmpF,$raTmpF2,$raTmpF3) = ("RXMtmp","R2XM","RX3M"); #my $raxFile = "RXMall";
#	my $raxFile2 = "RXMsyn";
	my $testModel = getProgPaths("modeltest");
	my $raxmlBin = getProgPaths("raxml");
	my $trDist = getProgPaths("treeDistScr");
	my $raxD = $outTree;$raxD=~m/([^\/]+$)/; my $trNm = $1;$trNm=~ s/\..*$// if ($trNm =~ m/\./);$raxD =~ s/[^\/]+$/${trNm}tmp\//;
	my $raxLogD= $outTree; $raxLogD =~ s/[^\/]+$/${trNm}Logs\//;
	#die "$raxD";
	system "rm -rf $raxD;" unless ($cont);
	system "mkdir -p $raxD";
	system "mkdir -p  $raxLogD" unless (-d $raxLogD);
	my $outGrpRXML = "";
	if ($outgroup ne ""){$outGrpRXML = "-o $outgroup";}
	my $raxDef = " --silent -m GTRGAMMA -p 312413 ";
	if ($useAA){
		$raxDef = " --silent -m PROTGAMMALG -p 312413 ";
	}
	
	#raxml - on all sites
	#die "$outTree\n";
	if (!$cont || !-e $outTree){
		my $tcmd="";
		$tcmd =  "$raxmlBin -T$ncore -f d -s $mAli $raxDef -n $raTmpF -w $raxD $outGrpRXML > $raxLogD/ini.log\n";
		#die $tcmd."\n";
		my $expTree1 = "$raxD/RAxML_bestTree.$raTmpF";
		$expTree1 = "$raxD/RAxML_result.$raTmpF" if ($useAA);
		if (!-e $expTree1){print "Calculating ML tree..\n" ; system "rm -f $raxD/*$raTmpF";system $tcmd;}
		if (!-e $expTree1 && !$useAA){#failed.. prob optimization problem, use other tree instead
			print "Calculating ML tree with GTRGAMMAI model..\n" ;
			$raxDef = " --silent -m GTRGAMMAI -p 3512413 ";	system "rm -f $raxD/*$raTmpF";
			systemW "$raxmlBin -T$ncore -f d -s $mAli $raxDef -n $raTmpF -w $raxD $outGrpRXML > $raxLogD/optimized.log\n";
		}elsif (!-e $expTree1){
			die "Failed $tcmd\n";
		}
		#decide which support vals to calc
		if ($bootStrap==0){
			
			print "Calculating SH support tree..\n";
			systemW "$raxmlBin -T$ncore -f J -s $mAli $raxDef -n $raTmpF3 -w $raxD -t $expTree1 $outGrpRXML > $raxLogD/optimized.log\n";
			system "mv $raxD/RAxML_fastTreeSH_Support.$raTmpF3 $outTree";
		} else {
			my $done=0; my $partBoots = "";
			if (-e "$raxD/RAxML_bootstrap.$raTmpF2" || -e "$raxD/RAxML_bootstrap.prev.x"){
				if (-e "$raxD/RAxML_bootstrap.$raTmpF2"){
					my $tmp = `wc -l $raxD/RAxML_bootstrap.$raTmpF2`; $tmp =~ m/^(\d+) /; $done += $1;
				}
				if (-e "$raxD/RAxML_bootstrap.prev.x"){
					my $tmp = `wc -l $raxD/RAxML_bootstrap.prev.x`; $tmp =~ m/^(\d+) /; $done += $1;
					$partBoots="$raxD/RAxML_bootstrap.prev.x";
				}
				if ($done < $bootStrap){system "cat $raxD/RAxML_bootstrap.$raTmpF2 >> $raxD/RAxML_bootstrap.prev.x;rm $raxD/RAxML_bootstrap.$raTmpF2" if (-e "$raxD/RAxML_bootstrap.$raTmpF2"); }
				
			}
			#die "$done < $bootStrap\n";
			if (!$cont || $done < $bootStrap){
				system "rm -f $raxD/*$raTmpF2"; 
				#if ($partBoots ne ""){system "cat $partBoots >> $raxD/RAxML_bootstrap.$raTmpF2";}
			}
			if ($done >= $bootStrap){
				system "mv $partBoots  $raxD/RAxML_bootstrap.$raTmpF2";
			}
			$bootStrap = $bootStrap - $done ;
			if (!-e "$raxD/RAxML_bootstrap.$raTmpF2"){
				my $bootCmd.= "-N $bootStrap -b ". int(rand(10000));
				print "Calculating bootstrap tree with $bootStrap bootstraps..\n";
				systemW  "$raxmlBin -T$ncore -f d  -s $mAli $raxDef -n $raTmpF2 -w $raxD $bootCmd $outGrpRXML > $raxLogD/ini.log\n";
				if ($partBoots ne ""){system "cat $partBoots >> $raxD/RAxML_bootstrap.$raTmpF2";}
			}
			#die;
			system "rm -f $raxD/*$raTmpF3";
			systemW  "$raxmlBin -T$ncore -f b  -s $mAli $raxDef -n $raTmpF3 -w $raxD -t $raxD/RAxML_bestTree.$raTmpF -z $raxD/RAxML_bootstrap.$raTmpF2 $outGrpRXML > $raxLogD/ini.log\n";
			system "mv $raxD/RAxML_bipartitions.$raTmpF3 $outTree";
		
		}
	} else {
		print "Tree already present\n";
	}
	#and create distance matrix
#	print "Calculating ML matrix from tree..\n";
#	system "rm -f $raxD/*all"; 
#	systemW  "$raxmlBin -T$ncore -f x -s $mAli $raxDef -n all -w $raxD -t $raxD/RAxML_bestTree.$raTmpF $outGrpRXML > $raxLogD/.log\n";
	#die "$outDist\n";
#	system "mv $raxD/RAxML_distances.all $outDist";
	
	print "Calculating tree distance matrix from tree..\n";
	my $outDist = $outTree; $outDist =~ s/\.[^\.]+$/\.dist/;
	#system "rm -f $raxD/*all"; 
	#this perl script has some errors!!
	#system  "$trDist $outTree > $outDist\n";
#	system "mv $raxD/RAxML_distances.all $outDist";

	
	#clean up
	system "rm -rf $raxD";
}

sub getGenoGenes{
	my $refG = $_[0];
	my $rewrite =0;
	$rewrite = $_[1] if (@_ >= 2);
	my $cores = 1;
	$cores = $_[2] if (@_ >= 3);
	my $pprodigalBin = getProgPaths("pprodigal");
	my $proteins = $refG; my $ntGenes = $refG; $proteins =~ s/\.[^\.]+$/\.genes\.faa/;
	$ntGenes =~ s/\.[^\.]+$/\.genes\.fna/;
	die "Can't find ref genome $refG\n" unless (-e $refG || (-e $ntGenes && -e $proteins));
	# ("$pprodigalBin -i $inputBac -o $tmpGene/genes$bacmark.gff -a $tmpGene/proteins.faa -d $tmpGene/genes.fna -f $output_format_prodigal -p meta -T $numThr \n");
	my $prodigal_cmd .= "$pprodigalBin -i $refG -a $proteins -d $ntGenes -f gff -p single -T $cores > /dev/null\n";
	#die $prodigal_cmd."\n$rewrite || !-e $proteins || !-e $ntGenes\n";
	system $prodigal_cmd if ($rewrite || !-e $proteins || !-e $ntGenes);
	return ($ntGenes,$proteins);
}
sub getFMG{
	my ($oDess,$proteins,$genesNT) = @_;
	my $redo=0;my $ncore=1;my $rename="";
	$ncore = $_[3] if (@_ >= 4);
	$redo = $_[4] if (@_ >= 5);
	$rename = $_[5] if (@_ >= 6);
	
	if ($proteins eq ""){die "getFMG:: need protein file\n$proteins\n";}
	if ($oDess eq ""){
		$proteins=~m/^(.*)\/([^\/]+)$/;
		$oDess = $1;
		if ($2 =~ m/(.*)\b\.genes\b?\.faa$/){
			$oDess .= "/FMGs$1/";
		} else { die "getFMG:: protein name needs to end on .faa\n$2\n";}
	}
	#die $proteins."\n".$oDess."\n";
	my $GenomeDir = $proteins; $GenomeDir =~ s/[^\/]+$//;
	system "mkdir -p $oDess" unless (-d $oDess);
	#die $cmd;
	if ( $redo  || (!-s "$oDess/COG0012.faa" && !-s "$oDess/COG0016.faa")){
		my $FMGd = getProgPaths("FMGdir");#"/g/bork5/hildebra/bin/fetchMG/";
		my $FMGrwkScr = getProgPaths("FMGrwk_scr");
		my $cmd = "";
		$cmd = "perl $FMGd/fetchMG.pl -m extraction -o $oDess -l $FMGd/lib -t $ncore -d $genesNT $proteins "; # -x $FMGd/bin  -b $FMGd/lib/MG_BitScoreCutoffs.uncalibrated.txt
		#print $cmd;
		systemW $cmd;
		system "rm  -rf $GenomeDir/*.cidx $oDess/temp $oDess/hmmResults";
		system "$FMGrwkScr $oDess" if ($redo);
	}
	#die "$redo\n";
	return $oDess;
}
#reads the FMGs nicely by cat, creating good header name
sub readFMGdir{
	my $inD = $_[0];
	my $denom = "unknown";
	if (@_ > 1){
		$denom = $_[1];
	} else {
		$inD =~ /\/([^\/]*)\/?$/;
		$denom = $1;
	}
	my $SaSe = "|"; #separator
	if (@_>2){
		$SaSe = $_[2];
	}
	
	#die "$SaSe\n";
	my %catGe;my %FMnt; my %FMaa;
	my $allFine=0;
	opendir(DIR, $inD) or die "Can't find: $inD\n";	
	my @fnas = sort ( grep { /COG\d+\.fna/  && -e "$inD/$_"} readdir(DIR) );	close(DIR);
	my @fnaT=("","");my $aaT="";my $cnt=0;
	foreach my $nf (@fnas){
		my ($hr,$OK) = readFasta("$inD/$nf",1);
		my %FNA = %{$hr};
		$nf =~ s/\.fna$/\.faa/;
		($hr,$OK) = readFasta("$inD/$nf",1);
		my %FAA = %{$hr};
		$nf =~ s/\.faa$//;
		#decide on longest gene variant
		my $long=0; my $lk="";
		foreach my $k (keys %FNA){
			my $lfaa=length($FAA{$k});
			if ($lfaa > $long){$lk=$k;$long=$lfaa;}
		}
		next if ($lk eq "");
		my $ng = "$denom$SaSe$nf";
		if ($SaSe eq "-1"){$ng="$denom.$cnt";}
		$FMnt{$ng} = $FNA{$lk};
		$FMaa{$ng} = $FAA{$lk};
		#die "$ng\n";
		$catGe{$ng} = $nf;
		$cnt++;
	}
	return (\%FMnt,\%FMaa,\%catGe);
}
sub renameFMGs{
	my ($oDess,$rename,$outN,$renameShrt) = @_;
	#create an "allFMG" file with the correct fasta header names
	my %catGe;
	my $allFine=0;
	opendir(DIR, $oDess) or die "Can't find: $oDess\n";	
	my @fnas = sort ( grep { /COG\d+\.fna/  && -e "$oDess/$_"} readdir(DIR) );	close(DIR);
	my @fnaT=("","");my $aaT="";
	foreach my $nf (@fnas){
		$nf =~ m/(.*)\.fna/;
		my $cat = $1;
		for (my $j=0;$j<2;$j++){
			$nf =~ s/\.fna/\.faa/ if ($j==1);
			my $hdFnd=0;
			open I ,"<$oDess$nf";
			while (my $li=<I>){
				if ($li =~ m/^>/){ 
					if ($hdFnd){print "Found > 1 seq in $nf\n";last;}
					$hdFnd=1;
					my $newH = $rename."_". $cat ;
					$fnaT[$j] .=  ">".$newH  . "\n";
					$catGe{$cat} = $newH;
					if ($li =~ m/>$newH/){$allFine=1;}
				} else {
					$fnaT[$j] .= $li;
				}
			}
			close I;
			if (!$allFine && $renameShrt==1){
				open O ,">$oDess$nf"; print O $fnaT[$j]; close O;
				$fnaT[$j]=""; 
			}
		}
	}
	if ($renameShrt==0){
		print "renaming FMGs to name $rename\n";
		open O,">$oDess/$outN.fna" or die "Can't open renameFMGs outfile $oDess/$outN.fna\n";
		print O $fnaT[0]; close O;
		open O,">$oDess/$outN.faa" or die "Can't open renameFMGs outfile $oDess/$outN.faa\n";
		print O $fnaT[0]; close O;
	}
	return (\%catGe);
}
sub getE100($ $ $ $){
	my ($oDess,$proteins,$genesNT,$ncore) = @_;
	my $hmmbin = getProgPaths("hmmsearch");#"/g/bork5/hildebra/bin/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch";
	my $essDB = getProgPaths("essentialHMM");#"/g/bork5/hildebra/bin/multi-metagenome-master/R.data.generation/essential.hmm";
	my $essEukDB = getProgPaths("essentialEUK");#"/g/bork3/home/hildebra/DB/HMMs/eukCore/eukCore.hmm"; #TODO
	systemW("mkdir -p $oDess");
	my $cmd = "";
	$cmd .= "$hmmbin --domtblout $oDess/assembly.hmm.orfs.txt -o $oDess/assembly.hmmsout.txt --cut_tc --cpu $ncore --notextw $essDB $proteins\n";
	$cmd .= "tail -n+4 $oDess/assembly.hmm.orfs.txt | sed \'s/\\s\\s*/ /g\' | sed \'s/^#.*//g\' | cut -f1,4 -d \" \" > $oDess/ess100.id.txt\n";

	if (!-s "$oDess/ess100.id.txt" && !-e "$oDess/e100split.sto"){
		print "\nDetecting 100 essential proteins\n";
		systemW $cmd ."\n";
	}
	my $protIDss = `cut -f1 -d " " $oDess/ess100.id.txt `;
	my $protClsTmp = `cut -f2 -d " " $oDess/ess100.id.txt `;
	my @protIDs = split("\n",$protIDss);
	my @protCls = split("\n",$protClsTmp);
	my $phr = readFasta("$proteins");
	my $ghr = readFasta("$genesNT");
	my %prots = %{$phr}; my %essProt;
	my %genes = %{$ghr};

	my %seen;
	#split 100 essentials into separate files for each gene class
	#my @unique = grep { ! $seen{$_}++ } @faculty;
	foreach ( grep { ! $seen{$_}++ } @protCls){
		open O,">$oDess/pe100_".$_.".faa";	close O;
		open O,">$oDess/ge100_".$_.".fna";	close O;
	}
	for (my$i=0;$i<@protIDs;$i++){
		my $id = $protIDs[$i];
		next if ($id eq "");
		if (!exists($prots{$id})){die "Can't find $id protein in file $proteins\n";}
		if (!exists($genes{$id})){die "Can't find $id protein in file $proteins\n";}
		open O,">>$oDess/pe100_".$protCls[$i].".faa" or die "Can't open $oDess/pe100_$protCls[$i].faa";
		print O ">$id\n$prots{$id}\n";
		close O;
		open O,">>$oDess/ge100_".$protCls[$i].".fna"or die "Can't open $oDess/ge100_$protCls[$i].faa";
		print O ">$id\n$genes{$id}\n";
		close O;
	}
}
