#!/usr/bin/perl
#relies on helpers\GC\annotateMGwSpecIs.pl to be run first to extract SpecI's..
#kind of the same as strain_within.pl, but this is more focused between samples..
#does not include MGS, so does actually not need to be run
#./FMGgenes2tree.pl [GCdir] [ncore] [extra guide file (e.g. MGS)]
#example ./FMGgenes2tree.pl /g/bork3/home/hildebra/data/SNP/GCs/alienGC2 12



use warnings;
use strict;
use threads ('yield',
                 'stack_size' => 64*4096,
                 'exit' => 'threads_only',
                 'stringify');

use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
use Mods::IO_Tamoc_progs qw(getProgPaths );
my $bts = getProgPaths("buildTree_scr");
my $neiTree = getProgPaths("neighborTree");
sub extractFNAFAA2genes;


#read pairs
my $GCd = $ARGV[0];
my $numCores = $ARGV[1];
my $treeFile = "";
$treeFile = $ARGV[2] if (@ARGV > 2);
my $justSubmit =0;
$justSubmit = $ARGV[3] if (@ARGV > 3);

my $xtraGuids = 
my $outD = "$GCd/FMG/phylo/";

my $QSBoptHR = emptyQsubOpt(1,"");
my %QSBopt = %{$QSBoptHR};

system "mkdir -p $outD" unless (-d $outD);

print "\n-------------------------------\nCreating within species strains for FMGs in $GCd\n";
print "Using tree $treeFile to create automatically outgroups\n" if ($treeFile ne "");
print "Only submission mode\n" if ($justSubmit);


my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;
#$mapF = $GCd."LOGandSUB/inmap.txt" if ($mapF eq "");
my ($hr1,$hr2) = readMapS($mapF,-1);
my %map = %{$hr1}; my %AsGrps = %{$hr2};
my %AGlist;
#get all samples in assembly group, but only last in mapgroup
my @samples = @{$map{opt}{smpl_order}};
my $fileAbsent = 0;
foreach my $smpl (@samples){
	next if ($map{$smpl}{AssGroup} eq "-1");
	my $cAssGrp = $map{$smpl}{AssGroup};
	my $cMapGrp = $map{$smpl}{MapGroup};
	#print "$smpl $cAssGrp $cMapGrp\n";
	$AsGrps{$cMapGrp}{CntMap} ++;
	if (!exists($AsGrps{$cMapGrp}{CntMap})){ die "Can;t find CntMap for $cMapGrp";}
	next if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} );

	push(@{$AGlist{$cAssGrp}} , $smpl);
	
	if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} ){
		next;
	}

	#check if SNP file is present
	next if ($justSubmit);
	my $cD = $map{$smpl}{wrdir}."/";
	my $tar = $cD."/SNP/genes.shrtHD.SNPc.MPI.fna.gz";
	my $tar2 = $cD."/SNP/proteins.shrtHD.SNPc.MPI.faa.gz";
	if (!-e $tar || !-e $tar2){
		print "Can't find SNP  file: $cD\n" ;
		$fileAbsent = 1;
	}
}
if ($fileAbsent){
	print "Not all required input present\nAborting FNGgenes2tree.pl\n" ;
	exit(0);
} else {
	print "All samples have SNP calls\n";
}
#foreach (sort keys %AGlist) {   print "$_ : @{$AGlist{$_}}\n";}die;

my %SIgenes; my %Gene2COG;
my %replN; 
my %allFNA; my %allFAA; #big hash with all genes in @allGenes
my %gene2genes;
my %cl2gene2;
my %FNAfmg; my %FAAfmg;


open I,"<$GCd/FMG/gene2specI.txt" or die "Can't open gene 2 specI file:\n$GCd/FMG/gene2specI.txt\n";
while (my $line = <I>){
	chomp $line;
	my @spl = split /\t/,$line;
	$SIgenes{$spl[1]}{$spl[2]} = $spl[0];
	$Gene2COG{$spl[0]} = $spl[2];
}
close I;


if (!$justSubmit){
	($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",0);my %cl2gene = %{$hr2}; $hr1 = {};
	#stores alt names (wihtout M4__ at end)
	#read binning based on SpecI's 
	my @allGenes;
	foreach my $gene (keys %Gene2COG){
		my $geneStr = $cl2gene{$gene};
		$geneStr =~ s/>//g;
		my @genegenes = split /,/,$geneStr;
		$cl2gene2{$gene } = \@genegenes;
		push(@allGenes,@genegenes);
	}
	undef %cl2gene; #lessen mem

	print "Extracting ".scalar @allGenes." genes from dirs\n";
	#and extract the corresponding fna/ faa from every other dir.. 
	extractFNAFAA2genes(\@allGenes);
	print "\nGene extraction finished\n";
	@allGenes = ();
	#	exit(0);

	#also read FMG gene seqs
	print "reading FMG ref genes..";
	my $hr = readFasta("$GCd/FMG/COG*.faa"); %FAAfmg = %{$hr};
	$hr = readFasta("$GCd/FMG/COG*.fna"); %FNAfmg = %{$hr};
	print "done\n";
}


#go through every SpecI;
my @specis = sort(keys(%SIgenes));
my $cnt=0; my $SaSe = '|';
foreach my $SI (@specis){ #loop creates per specI file structure to run buildTreeScript on..
	#PART I: create fasta files required by tree
	my $outD2 = "$outD/$SI/";
	my $multiSmpl=0;
	system "mkdir -p $outD2" unless (-d $outD2);
	my $FNAtf = "$outD2/allFNAs.fna"; my $FAAtf = "$outD2/allFAAs.faa";
	my $CATtf = "$outD2/cats4ete.txt"; my $Linkf = "$outD2/link2GC.txt";
	
	
	#include outgroup?
	my $outgS = "";my $OG = "";
	if ($treeFile ne ""){
		my $OG1 = `$neiTree $treeFile $SI`; chomp $OG1;
		#print "$neiTree $treeFile $SI\nX${OG1}X\n";
		my @sspl = split /\s/,$OG1; $OG = $sspl[0];my $cntX=1;
		next if (@sspl == 0);
		while (!exists($SIgenes{$OG})){$OG=$sspl[$cntX];$cntX++;last if ($cntX > @sspl);}
		unless (exists($SIgenes{$OG})){
			print "can't find speci $OG\n$OG1\n";
			$OG="";
		} else {
			$outgS = "-outgroup $OG " ;
		}
		#next;
	}
	my $Tcmd= "$bts -fna $FNAtf -aa $FAAtf -smplSep '\\$SaSe' -cats $CATtf -outD $outD2 $outgS -runIQtree 1 -runFastTree 0 -cores $numCores  ";
	$Tcmd .= "-AAtree 0 -bootstrap 000 -NTfiltCount 300 -NTfilt 0.1 -NTfiltPerGene 0.5 -runRaxMLng 1 -minOverlapMSA 2 -MSAprogram 2 -AutoModel 0 -iqFast 1 \n";
	if ($justSubmit){
		die "Can't find required input files:\n$FNAtf\n$FAAtf\n$CATtf\n" unless (-e $FNAtf && -e $FAAtf && -e $CATtf);
		my ($dep,$qcmd) = qsubSystem($outD2."treeCmd.sh",$Tcmd,$numCores,"1G","FT$cnt","","",1,[],$QSBoptHR);
		$cnt ++;
		#die "$outD2";
		next;
	}
	
	#and fasta/faa/cat files..
	open OC,">$CATtf" or die "Can't open cat file $CATtf\n";
	open OF,">$FNAtf" or die "Can't open NT file $FNAtf\n";
	open OA,">$FAAtf" or die "Can't open AA file $FAAtf\n";
	open OL,">$Linkf" or die "Can't open link file $Linkf\n";
	foreach my $cog (sort keys %{$SIgenes{$SI}}){
		my $tar = $SIgenes{$SI}{$cog};
		my @genes = @{$cl2gene2{$tar}};
		#print OC "$cog";#.join("\t",@genes)."\n";
		my $new=1; my %smplsSeen; my $multiSmpl1 = 0;
		print OL "$cog\t$tar\t".scalar @genes . "\t".join(",",@genes)."\n";
		if ($OG ne "" && exists($SIgenes{$OG}{$cog})){
			die "can't find gene $SIgenes{$OG}{$cog}" unless (exists($FNAfmg{$SIgenes{$OG}{$cog}}));
			my $ng = "$OG$SaSe$cog";
			print OF ">$ng\n$FNAfmg{$SIgenes{$OG}{$cog}}\n";
			print OA ">$ng\n$FAAfmg{$SIgenes{$OG}{$cog}}\n";
			if ($new){ print OC "$ng";$new=0;
			} else { print OC "\t$ng"; }
		}
		foreach my $gX (@genes){
			my @genes2 = ($gX);
			if (exists($gene2genes{$gX})){
				#here I add another layer that each gene is linked to its spawned childs..
				@genes2 = (@{$gene2genes{$gX}},$gX);
				#die "mrg $gX : @genes2\n";
			} 
			#for (my $i=0;$i<@genes2;$i++){
			#	$genes2[$i] =~ s/M\d+__//;
			#}
			foreach my $g (@genes2){
				#my $mm=0;$mm=1 if ($allFAA{$g} =~ m/X/);print "$allFAA{$g}\n" if ($mm);
				if (!exists($allFAA{$g})){die "$g : $cog : $SI have NT but not AA\n" if (exists($allFNA{$g}));print "miss: $g\n"; next;}
				my $strCpy = $allFAA{$g};
				my $AAlen = length($allFAA{$g});
				next if ($AAlen == 0);
				my $num1 = $strCpy =~ tr/[\-Xx]//;
				next if ($num1 >= ($AAlen-1));
				
				#filtering is handled by buildtree only
				#my $count = $allFAA{$g} =~ tr/X//;
				#print "$allFAA{$g}\n" if ($mm);
				#next if ($count >= $AAlen);
				$g =~ m/(^.*)__/; my $curSmpl=$1;
				if (exists($replN{$curSmpl})){
					#print "exist: $curSmpl = $replN{$curSmpl}\n";
					$curSmpl = $replN{$curSmpl};
				}
				unless(exists($smplsSeen{$curSmpl})){
					$smplsSeen{$curSmpl} =1;
				} else {
					$multiSmpl++;
					$curSmpl.="_".$multiSmpl1;
					$multiSmpl1++;
				}
				die "$g not found\n" unless (exists($allFNA{$g}));
				my $ng = "$curSmpl$SaSe$cog"; #must contain 2 informations: 1)sampleID 2)COG 
				print OF ">$ng\n$allFNA{$g}\n";
				print OA ">$ng\n$allFAA{$g}\n";
				if ($new){ print OC "$ng";$new=0;
				} else { print OC "\t$ng"; }
			}
		}
		print OC "\n";
	}
	close OC;close OA; close OF;
	#note done somewhere how many genes these actually are..
	system "echo \"$OG\" > $outD2/numG_". scalar(keys %{$SIgenes{$SI}})."_smpls_$multiSmpl";
	system "echo \"OG:$OG\" > $outD2/data.log";
	
	if ($multiSmpl){
		print "$SI: multiSmpls: $multiSmpl\n";
	} else {next;}

	#PART II: qsub tree build command
	
	#die "$cmd\n" if ($cnt ==10);
	my ($dep,$qcmd) = qsubSystem($outD2."treeCmd.sh",$Tcmd,$numCores,"1G","FT$cnt","","",1,[],$QSBoptHR);
	$cnt ++;
	#die $outD2."treeCmd.sh\n";

}

print "\nAll done\n";





sub readGenesSample_Singl{
	#go into curSpl dir and extract genes..
	my @subG = @{$_[0]};
	die "regex failed: $subG[0]\n" unless ($subG[0] =~ m/^(.*)__/);
	my $sd = $1;
	my $sd2 = $sd;
	#print "$sd ";
	if (exists(  $map{altNms}{$sd}  )){
		$sd2 = $map{altNms}{$sd}; 
		$replN{$sd} = $sd2;
	}
	unless (exists ($map{$sd2}) ) {
		print "Can't find map entry for $sd\n"; die;
	}
	
	#find out if other samples are in the same assmblGrp..
	my @subSds = ($sd2);
	my $cAssGrp = $map{$sd2}{AssGroup};
	if (exists($AGlist{$cAssGrp})){
		@subSds = @{$AGlist{$cAssGrp}};
	}
	print "YY @subSds : $sd2\n";
	foreach my $sd3 (@subSds){
		my $cD = $map{$sd3}{wrdir}."/";
		my $rename = 0;
		$rename = 1 if ($sd2 ne $sd3);
		print "r:$rename $sd3  ";
		die "Assembly path missing $cD\n" unless (-e "$cD/assemblies/metag/assembly.txt");
		my $metaGD = `cat $cD/assemblies/metag/assembly.txt`; chomp $metaGD;
		if ($metaGD eq ""){die "assembly not available: $cD $sd3";}
		#get NT's
		#my $tar = $metaGD."genePred/genes.shrtHD.fna";
		my $tar = $cD."/SNP/genes.shrtHD.SNPc.MPI.fna.gz";
		unless (-e $tar){
			print "\n=====================================\nCan't find nt file $tar\n=====================================\n";
			next;
		}
		#return(%locFAA);
		#print "$tar\n";
		my $hr = readFasta($tar,1);
		my %FNA = %{$hr};
		#my @kk = keys(%FNA);		die scalar(@kk);
		foreach my $ge (@subG){
			#die "Can't find gene $ge\n" unless (exists($FNA{$ge}));
			#die"YES $ge $tar" if ($ge =~ m/C1404_L=8071=_3/);
			next unless (exists($FNA{$ge}));
			my $tmp = $FNA{$ge};
			my $ge2= $ge;
			if ($rename){
				#die "${sd2}__/${sd3}__\n";
				$ge2 =~ s/${sd}__/${sd3}__/;
				push(@{$gene2genes{$ge}},$ge2); #only save subset..
			}
			$allFNA{$ge2} = $tmp; #allFNA
		}
		#get AA's
		#$tar = $metaGD."genePred/proteins.shrtHD.faa";
		$tar = $cD."/SNP/proteins.shrtHD.SNPc.MPI.faa.gz";
		$hr = readFasta($tar,1);
		%FNA = %{$hr};
		foreach my $ge (@subG){
			#die "Can't find gene $ge\n" unless (exists($FNA{$ge}));
			next unless (exists($FNA{$ge}));

			my $tmp = $FNA{$ge};
			my $ge2= $ge;
			if ($rename){
				if ($ge2 =~ s/${sd}__/${sd3}__/){
					#push(@{$gene2genes{$ge}},$ge2); #only save subset..
				} else {
					die "could not replace $sd with $sd3 in string $ge2\n";
				}
			}
			$allFAA{$ge2} = $tmp; #allFAA

		}
	}
}


#this routine hast to get genes out of each sample, that are needed
#and save them to be later written per specI
sub extractFNAFAA2genes{
	my ($ar) = @_;
	my %allIDs = %{$hr1};
	my @allGenes1 = sort(@{$ar});
	my $maxT = 3;
	
	my $curSmpl = ""; my @subG;
	my @thrs; my @thrsUse; my $t=0;
	foreach my $g (@allGenes1){
		$g =~ m/(^.*)__/;
		if ($1 ne $curSmpl){
			if ($curSmpl eq ""){
				$curSmpl = $1;
				push @subG, $g;
				next;
			}
			#print "$t  ";
			#readGenesSample(@subG);
			if (0){
				if ($thrsUse[$t]){
					$thrsUse[$t] = 0;my %lret = $thrs[$t]->join();
					@allFNA{keys %{$lret{NT}}} = values %{$lret{NT}};@allFAA{keys %{$lret{AA}}} = values %{$lret{AA}};
				}
				($thrs[$t]) = threads->create(\&readGenesSample,@subG); $thrsUse[$t] = 1;
			} else {
				readGenesSample_Singl(\@subG);
			}
			#read cycle for sample finished, move to next
			@subG = ($g); 
			$curSmpl = $1;
			$t++;
			if ($t>=$maxT){$t=0;}

		} else { #just add to count...
			push @subG, $g;
		}
	}
	readGenesSample_Singl(\@subG);
	if (0){
		for ($t=0;$t<$maxT;$t++){
			if ($thrsUse[$t]){
				$thrsUse[$t] = 0;
				my %lret = $thrs[$t]->join();
				@allFNA{keys %{$lret{NT}}} = values %{$lret{NT}};
				@allFAA{keys %{$lret{AA}}} = values %{$lret{AA}};
			}
		}
	}
}



sub readGenesSample{
	#go into curSpl dir and extract genes..
	my @subG = @_;#@{$_[0]};
	die "dont use this readGenesSample\n";
	die "regex failed: $subG[0]\n" unless ($subG[0] =~ m/^(.*)__/);
	my $sd = $1;
	my $sd2 = $sd;
	my %locFAA;  $locFAA{NT} = {}; $locFAA{AA} = {}; 
	print "$sd ";
	if (exists(  $map{altNms}{$sd}  )){$sd2 = $map{altNms}{$sd}; }
	unless (exists ($map{$sd2}) ) {
		print "Can't find map entry for $sd\n"; die;
	}
	my $cD = $map{$sd2}{wrdir}."/";
	die "Assembly path missing $cD\n" unless (-e "$cD/assemblies/metag/assembly.txt");
	my $metaGD = `cat $cD/assemblies/metag/assembly.txt`; chomp $metaGD;
	if ($metaGD eq ""){die "assembly not available: $cD $sd";}
	#get NT's
	#my $tar = $metaGD."genePred/genes.shrtHD.fna";
	my $tar = $cD."/SNP/genes.shrtHD.SNPc.MPI.fna.gz";
	die "Can't find nt file $tar\n" unless (-e $tar);
	#return(%locFAA);
	#print "$tar\n";
	my $hr = readFasta($tar,1);
	my %FNA = %{$hr};
	foreach my $ge (@subG){
		die "Can't find gene $ge\n" unless (exists($FNA{$ge}));
		$locFAA{NT}{$ge} = $FNA{$ge}; #allFNA
	}
	#get AA's
	#$tar = $metaGD."genePred/proteins.shrtHD.faa";
	$tar = $cD."/SNP/proteins.shrtHD.SNPc.MPI.faa.gz";
	$hr = readFasta($tar,1);
	%FNA = %{$hr};
	foreach my $ge (@subG){
		die "Can't find gene $ge\n" unless (exists($FNA{$ge}));
		$locFAA{AA}{$ge} = $FNA{$ge}; #allFAA
	}
	return (%locFAA);
}
