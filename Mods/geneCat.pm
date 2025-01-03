package Mods::geneCat;
use warnings;
use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(  systemW gzipopen readFasta);
use Mods::FuncTools qw( readGene2Func);
use Mods::math qw( meanArray);

use Exporter qw(import);
our @EXPORT_OK = qw(readGeneIdx readGeneIdxSpl correlation calculate_spearman_correlation read_matrix checkAntiOcc
				readGene2tax createGene2MGS sortFNA readMG_LCA
				attachProteins attachProteins2 attachProteins3 );


sub readMG_LCA{
	my $markerGdir = $_[0];
	my @MGcats = ();
	if (@_ > 1){
		my $MGcatsAref = $_[1] ;
		@MGcats = @{$MGcatsAref};
	}
	my $depth = 6;
	$depth = $_[2] if (@_ > 2);
	#read in MG cats anew..
	if (@MGcats == 0){
		die "automatic MG category detection not implemented in :::readMG_LCA\n";
	}
	my $LCAfiles=0; my $LCAentries=0; my %LCA; my $LCAhead="";
	foreach my $GTcat (@MGcats){
		my $GTLCA = "$markerGdir/$GTcat.LCA";
		next unless (-e $GTLCA);
		$LCAfiles++; my $lcnt=0;
		open I ,"<$GTLCA" or die $!;
		while (<I>){
			$lcnt++;
			chomp; 
			my @spl = split /\t/;
			my $id = shift @spl;
			if ($lcnt==1){
				$LCAhead = join(";",@spl);next;
			}
			#m/(^\S+)\s/;
			while (@spl < $depth){push(@spl,"?");}
			#$LCA{$id} = join(";",@spl);
			$LCA{$id} = \@spl;
			$LCAentries++;
		}
		close I;
	}
	print "Read $LCAentries LCA assignments from $LCAfiles .LCA files\n";
	return \%LCA;
}

sub attachProteins3{
	my ($curSmpl,$prF,$protIn,$hrGI,$SEP) = @_;
	my %gene2num = %{$hrGI};
	die "Protein file does not exist: $protIn\n" unless (-e $protIn);
	my $hr = readFasta($protIn);
	my %fas = %{$hr};
	#print "$protIn\n";
	#my $protStore = `cat $inT | xargs $samBin faidx $protIn`;
	#die length($protStore)."\n";
	#my @prots = @{$inT};
	my $tmpStr = "";
	foreach my $pr (keys %gene2num){
		#$pr = substr $pr,1; #s/^>//;
		#1 identify protein
		my $seq = "";
		my $protHd = $curSmpl.$SEP.$pr;
		if (!exists($fas{$protHd})){
			print "Can't find $protHd in $protIn\n";
		} else {
			$seq = $fas{$protHd};
		}
		unless(exists($gene2num{$pr})){die "can not identify $pr gene in index file while rewritign prot names\n$protIn\n";}
		#print "$gene2num{$spl[0]}\n $pr\n";
		foreach my $prX (@{$gene2num{$pr}}){
			$tmpStr .= ">".$prX."\n".$seq."\n";
		}
	}
	open Oe,">>$prF" or die "Can't open $prF\n";
	print Oe $tmpStr;
	close Oe;
	$tmpStr="";
}


sub attachProteins2{
	my ($inT,$prF,$protIn) = ($_[0],$_[1],$_[2]);
	
	my %gene2num; my $doRename=0;
	if (@_ > 3){
		$doRename=1;
		my $hr = $_[3];
		%gene2num = %{$hr};
	}
	die "Protein file does not exist: $protIn\n" unless (-e $protIn);
	my $hr = readFasta($protIn);
	my %fas = %{$hr};
	#print "$protIn\n";
	#my $protStore = `cat $inT | xargs $samBin faidx $protIn`;
	#die length($protStore)."\n";
	#my @prots = @{$inT};
	my $tmpStr = "";
	my @inTa = @{$inT};
	foreach my $pr (@inTa){
		#$pr = substr $pr,1; #s/^>//;
		#1 identify protein
		my $seq = "";
		if (!exists($fas{$pr})){
			print "Can't find $pr in $protIn\n";
		} else {
			$seq = $fas{$pr};
		}
		if ($doRename){
			unless(exists($gene2num{$pr})){die "can not identify $pr gene in index file while rewritign prot names\n$protIn\n";}
			#print "$gene2num{$spl[0]}\n $pr\n";
			foreach my $pr (@{$gene2num{$pr}}){
				$tmpStr .= ">".$pr."\n".$seq."\n";
			}
		} else {
			$tmpStr .= ">".$pr."\n".$seq."\n";
		}
	}
	open Oe,">>$prF" or die "Can't open $prF\n";
	print Oe $tmpStr;
	close Oe;
	$tmpStr="";
}

sub attachProteins{#"$basD/tmp.txt",$protF){
	my ($inT,$prF,$protIn) = ($_[0],$_[1],$_[2]);
	
	my %gene2num; my $doRename=0;
	if (@_ > 3){
		$doRename=1;
		my $hr = $_[3];
		%gene2num = %{$hr};
	}
	my $samBin = getProgPaths("samtools");
	#the old way
	#systemW("cat $basD/tmp.txt | xargs samtools faidx $protIn  >> $prF");
	die "Protein file does not exist: $protIn\n" unless (-e $protIn);
	#print "$protIn\n";
	my $protStore = `cat $inT | xargs $samBin faidx $protIn`;
	#die length($protStore)."\n";
	my @prots = split />/,$protStore;
	open O,">>$prF" or die "Can't open $prF\n";
	foreach my $pr (@prots){
		#1 identify protein
		next if (length($pr) <= 1);
		my @spl = split /\n/,$pr;
		if ($doRename){
			unless(exists($gene2num{$spl[0]})){die "can not identify $spl[0] gene in index file while rewritign prot names\n$protIn\n";}
			#print "$gene2num{$spl[0]}\n $pr\n";
			my $pr = shift @{$gene2num{$spl[0]}};
			$spl[0] = ">".$pr;
		} else {
			$spl[0] = ">$spl[0]";
		}
		print O join("\n",@spl)."\n";
	}
	close O;
	#die "cat $inT | xargs $samBin faidx $protIn\n";
}


				
				
sub sortFNA{
	my ($bdir,$tag,$srtMode,$tmpDir,$nc) = @_;
	my $rcmd="";
	if (!-e "$bdir/$tag.srt.fna"){
		if ($srtMode ){
			my $vsBin = getProgPaths("vsearch");

			if (!-e "$bdir/$tag.fna"){die "can't find $bdir/$tag.fna\n";}
			#$cmd .= "\nmv $bdir/compl.fna $bdir/compl.unsrt.fna\n";
			
			#$rcmd .= "awk '/^>/ {printf(\"\%s\%s\\t\",(N>0?\"\\n\":\"\"),\$0);N++;next;} {printf(\"%s\",\$0);} END {printf(\"\\n\");}'  $bdir/$tag.fna  | awk -F '\\t' '{printf(\"\%d\\t\%s\\n\",length(\$2),\$0);}' > $tmpDir/$tag.s1.fna \n";
			#$rcmd .= "sort -S 70\% -k1,1nr -T $tmpDir $tmpDir/$tag.s1.fna | cut -f 2- | tr \"\\t\" \"\\n\" > $bdir/$tag.srt.fna\n";
			#$rcmd .= "rm $tmpDir/$tag.s1.fna $bdir/$tag.fna\n\n";
			$rcmd .= "$vsBin --sortbylength $bdir/$tag.fna --output $bdir/$tag.srt.fna --threads $nc\n";
		} elsif (!$srtMode){
			$rcmd .= "mv $bdir/$tag.fna $bdir/$tag.srt.fna\n";
		}
	}
	return $rcmd;
}

#create gene2tax file required for MGS, based on typical list file, but add in COG info
sub createGene2MGS{
	my ($MGSfile,$GCd) = @_;
	my $outF = "$MGSfile.gene2MGS";
	return ($outF) if (-e $outF);
	
	#first read in COG assignments / gene
	my $hr = readGene2Func("$GCd","NOG"); my %COG = %{$hr};
	#keep some stats..
	my $COGcnt=0;my %MGS; my $COGnot=0;
	my %MGScnts;
	open I,"<$MGSfile" or die "Can't open MGS guide file: $MGSfile\n";
	open O,">$outF" or die "Can't open gen2MGS file: $outF\n";
	while (my $lin = <I>){
		chomp $lin; my @spl = split (/\t/,$lin,-1);
		if (@spl <2){die "incomplete entry in MGS guide file: @spl\n";}
		my @genes = split /,/,$spl[1];
		#not needed here..
		#$MGS{$spl[0]} = \@genes;
		my $cMGS =$spl[0];
		if (exists($MGS{$cMGS} )){ print "Found $cMGS double!\n";}
		$MGS{$cMGS} = 1;
		foreach my $x (@genes){
			#desired format: 23075792        specI_v2_0561   COG0201
			my $curCOG = "";
			if (exists($COG{$x})){
				$COGcnt++;
				$curCOG = $COG{$x};
				$MGScnts{$cMGS}++;
			} else {$COGnot++;}
			print O "$x\t$cMGS\t$curCOG\n";
		}
	}
	close I; close O;
	#report top lowest MGS
	print "Lowest represented MGS::\n";
	my @keys = sort { $MGScnts{$a} <=> $MGScnts{$b} } keys(%MGScnts); my $cntX=0;
	foreach my $kk(@keys){
		print "$kk: $MGScnts{$kk};   "; $cntX++; last if ($cntX > 5);
	}
	print "\n";
	
	print "Found ". scalar(keys%MGS) . " MGS, $COGcnt/".($COGnot+$COGcnt)." genes could be assigned to NOGs from $MGSfile\n";
	print "MGS2gene file created at $outF\n";
	return ($outF);
}


sub readGene2tax{
	my $inF = $_[0];
	my $limit = -1;
	$limit = $_[1] if (@_ > 1);
	my %SIgenes;my %Gene2COG;my %Gene2MGS;
	my %uniqs; my %cogPrio;
	#some stats
	my %totalTax; my $totalGenes=0; my $inclGenes=0;
	open I,"<$inF" or die "Can't open gene 2 tax (specI/MGS) file:\n$inF\n";
	my $curTax = ""; my $curTcnt=0;
	
	while (my $line = <I>){
		chomp $line;
		$totalGenes++;
		my @spl = split (/\t/,$line,-1);
		#only read a limited number of genes.. used for MGS to only take first few genes
		if ($limit>0 && exists($totalTax{$spl[1]}) && $totalTax{$spl[1]}  >= $limit){
			next;
		} 
	
		my $OG = $spl[2];
		if ($OG eq "" || $OG eq "-"){ #make artificial OG
			$uniqs{$spl[1]}++;
			$OG="uniq$uniqs{$spl[1]}";
			#die "$OG\n";
		}
		$inclGenes++;
		$totalTax{$spl[1]} ++;
		unless (exists($SIgenes{$spl[1]}{$OG})){#only register gene if COG is not already reserved..
			push(@{$cogPrio{$spl[1]}},$OG); ;
			$SIgenes{$spl[1]}{$OG} = $spl[0];
			$Gene2COG{$spl[0]} = $OG;
			$Gene2MGS{$spl[0]} = $spl[1];
		}
	}
	close I;
	print "Found ". scalar(keys %totalTax) ." groups with $inclGenes/$totalGenes included genes\n";
	
	#double check on low represented MGS
	my @keys = sort { $totalTax{$a} <=> $totalTax{$b} } keys(%totalTax);
	my $lcnt=0;
	print "5 lowest MGS are: \n" ; 
	foreach my $k (@keys){print "$k $totalTax{$k};\t";$lcnt++;  if ($lcnt>5){print "\n";last;}}
	
	return (\%SIgenes,\%Gene2COG,\%Gene2MGS,\%cogPrio);
}


sub checkAntiOcc{
	my ($glAR,$matHR) = @_; #(\@potDG,\%FMGm);
	my @GL  = @{$glAR}; my %mat = %{$matHR};
	my %res;
	if (@GL <= 1) {return \%res;}
	for (my $i=0;$i<@GL;$i++){
		die "can't find gene $GL[$i] in matrix\n" unless (exists($mat{$GL[$i]}));
		my @g1 = @{$mat{$GL[$i]}};
		
		my $dblOcc=0; my $singlOcc=0;
		for (my $j=($i+1);$j<@GL;$j++){
			die "can't find gene $GL[$j] in matrix\n" unless (exists($mat{$GL[$j]}));
			my @g2 = @{$mat{$GL[$j]}};
			#simple countup for non occurrence...
			for (my $k=0;$k<@g2;$k++){
				if ($g1[$k] != 0){
					if ($g2[$k] != 0){ $dblOcc++; } else {$singlOcc++;}
				} elsif ($g2[$k] != 0){
					if ($g1[$k] != 0){ $dblOcc++; } else {$singlOcc++;}
				} #else both zero.. ignore
			}
			#print "OCC = $singlOcc / $dblOcc\n";
			$res{$GL[$i]}{$GL[$j]} = $singlOcc / ($singlOcc + $dblOcc);
		}
	}
	return \%res;
}


sub readGeneIdx($){
	my ($in) = @_;
	my %ret; 
	#experimental optimization..
	keys %ret = 5e7;
	my $gCnt=0; my $dbl=0;
	open I,"<$in" or die "Can't read gene index file $in\n";
	while(my $line=<I>){
		next if ($line =~ m/#/);
		chomp $line;
		my @spl = split(/\t/,$line);
		if (exists($ret{$spl[2]})){
			print "Double entry found: $spl[2]\n$line\n"; $dbl++;
		}
		push(@{$ret{$spl[2]}}, $spl[0]);
		$gCnt++;
		#print $spl[2] ." ". $spl[0]."\n" if ($gCnt<10);
	}
	close I;
	print "Gene Index read ($gCnt)\n";
	if ($dbl>0){
		print "Double ($dbl) entries in $gCnt total entries.. \n";
	}
	return (\%ret,$gCnt);
}

sub readGeneIdxSpl($ $){
	my ($in,$splTerm) = @_;
	my %ret; 
	#experimental optimization..
	my $gCnt=0; my $dbl=0;
	open I,"<$in" or die "Can't read gene index file $in\n";
	while(my $line=<I>){
		next if ($line =~ m/#/);
		chomp $line;
		my @spl = split(/\t/,$line);
		my @spl2 = split($splTerm,$spl[2]);
		my $pre="xtraSmpls"; my $post = $spl[2];
		if (@spl2>1){
			$pre = $spl2[0];$post = $spl2[1];
		}
		unless (exists($ret{$pre})){
			keys %{$ret{$pre}} = 5e5;
		}
		if (exists($ret{$pre}{$post})){
			print "Double entry found: $spl[2]\n$line\n"; $dbl++;
		}
		push(@{$ret{$pre}{$post}}, $spl[0]);
		$gCnt++;
		#print $spl[2] ." ". $spl[0]."\n" if ($gCnt<10);
	}
	close I;
	print "Gene Index read ($gCnt)\n";
	if ($dbl>0){
		print "Double ($dbl) entries in $gCnt total entries.. \n";
	}
	return (\%ret,$gCnt);
}


sub read_matrix{
	my ($mF) = $_[0];
	my $SEP="\t";
	$SEP = $_[1] if (@_ > 1);
	my %incl; my $doIncl=0;
	if (@_ > 1){
		my $hr = $_[2]; %incl = %{$hr};
		$doIncl=1;
	}
	print "Reading matrix $mF";
	print ", with ".scalar(keys %incl)." subsets" if ($doIncl);
	my %oM;
	my $lcnt=0;
	#open I,"<$mF" or die "Can't open $mF\n";
	my ($I,$stat) = gzipopen($mF,"Matrix File",1);
	my $cnt=0;
	while (<$I>){
		chomp; $cnt++;
		my @row = split /$SEP/;
		my $ID = shift @row;
		next if ($doIncl && !exists($incl{$ID}));
		if ($cnt==1){
			$oM{header} = \@row;
		} else {
			#for (my $i=0;$i<@row;$i++){if ($row[$i] < 2){$row[$i]=0;} } 
			#for (my $i=0;$i<@row;$i++){$row[$i] = sqrt ($row[$i]);} 
			$oM{$ID} = \@row;
		}
		$lcnt++;
	}
	close $I;
	if ($lcnt <= 2){die "not enough lines in matrix $mF\n";}
	print " .. Done\n";
	return \%oM;
}

#correlation routines
sub convert_values_to_ranks {
	my $values = shift;

	# code below is slightly unintuitive, but we have to compute
	# average rank in cases where there are ties.
	my $idx = 0;
	my %sorted; # $sorted{$val} = [ $rank1, $rank2, $rank3 ];
	foreach my $val ( sort { $a <=> $b } @$values ) {
		push @{$sorted{$val}}, $idx;
		$idx++
	}

	# compute the average rank for a given value
	my %average_ranks =
		map { $_ => List::Util::sum(@{$sorted{$_}}) / scalar(@{$sorted{$_}}) }
		keys %sorted;

	# encode the values using average rank of the value in the list
	my @ranks = map { $average_ranks{$_} } @$values;

	if ( wantarray ) { return @ranks; }
	return \@ranks;
}

sub calculate_spearman_correlation {
	my $n1 = shift;
	my $n2 = shift;

	if ( scalar(@{$n1}) != scalar(@{$n2}) ) {
		die "Error: spearman correlation given two lists of unequal size!\n@{$n2}\n@{$n1}\n";
	}

	my $ranked_n1 = convert_values_to_ranks($n1);
	my $ranked_n2 = convert_values_to_ranks($n2);

	my $sum_diff_squared = 0;
	foreach my $idx ( 0 .. scalar(@$ranked_n1)-1 ) {
		my $diff = $ranked_n1->[$idx] - $ranked_n2->[$idx];
		$sum_diff_squared += $diff * $diff;
	}

	my $N   = scalar(@$ranked_n1);
	my $rho = 1 - ( 6 * $sum_diff_squared / ($N * ($N*$N-1)) );
	return $rho;
}
 
sub mean($) {
	   my ($x)=@_;
	   return if (@_ == 0);
	   my $num = scalar(@{$x}) - 1;
	   my $sum_x = '0';
	   for (my $i = 1; $i < scalar(@{$x}); ++$i){
		  $sum_x += $x->[$i];
	   }
	   my $mu_x = $sum_x / $num;
	   return($mu_x);
}
 
### ss = sum of squared deviations to the mean
sub ss {
   my ($x,$mean_x,$y,$mean_y)=@_;
   my $sum = '0';
   for (my $i=1;$i<scalar(@{$x});++$i){
     $sum += ($x->[$i]-$mean_x)*($y->[$i]-$mean_y);
   }
   return $sum;
}
 sub correl {
   my($ssxx,$ssyy,$ssxy)=@_;
   if ($ssyy==0 || $ssxy==0){return 0;}
   my $sign=$ssxy/abs($ssxy);
   my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
   return $correl;
}
sub correlation {
   my ($x,$y) = @_;
   my ($mean_x) = mean($x);
   my $mean_y = mean($y);
   my $ssxx=ss($x,$mean_x,$x,$mean_y);
   my $ssyy=ss($y,$mean_x,$y,$mean_y);
   my $ssxy=ss($x,$mean_x,$y,$mean_y);
   my $correl=correl($ssxx,$ssyy,$ssxy);
   my $xcorrel=sprintf("%.4f",$correl);
   return($xcorrel);
}
 