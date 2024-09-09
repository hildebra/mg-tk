package Mods::GenoMetaAss;
use warnings;
#use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths);

use Exporter qw(import);
our @EXPORT_OK = qw(convertMSA2NXS gzipwrite gzipopen renameFastaCnts renameFastqCnts readNCBItax   lcp prefix_find
		readMap readMapS getDirsPerAssmblGrp checkSeqTech is3rdGenSeqTech 
		renameFastHD  prefixFAhd parse_duration resolve_path
		clenSplitFastas getAssemblPath fileGZe fileGZs filsizeMB
		readClstrRev  readClstrRevGenes readClstrRevContigSubset readClstrRevSmplCtgGenSubset
		unzipFileARezip systemW is_integer 
		readGFF reverse_complement reverse_complement_IUPAC
		readFasta writeFasta readFastHD splitFastas 
		readTabByKey convertNT2AA runDiamond median mean quantile
		getRawSeqsAssmGrp getCleanSeqsAssmGrp addFileLocs2AssmGrp iniCleanSeqSetHR
		hasSuppRds
		 );#Binning Related



sub parse_duration {
	my $seconds = shift;
	my $hours = int( $seconds / (60*60) );
	my $mins = ( $seconds / 60 ) % 60;
	my $secs = $seconds % 60;
	return sprintf("00:00:%02d", $seconds) if $seconds < 60;
	return sprintf("%02d:%02d:%02d", $hours,$mins,$secs);
}

#sums up filsizes for 1 or several files; first arg can be path to files, all subsequent args need to be filenames
sub filsizeMB{
	my $totalMapSize=0;
	my $path="";
	if (-d $_[0]){$path=shift @_;}
	foreach my $fh (@_){ $totalMapSize += (-s "$path/$fh") / (1024 * 1024) if (-e "$path/$fh");}
	return $totalMapSize;
}

#check if file or file.gz exists
sub fileGZe{
	my ($fil) = @_;
	return 1 if (-e $fil);
	return 1 if (-e "$fil.gz");
	return 0;
}


#report file size, check if file or file.gz exists
sub fileGZs{
	my ($fil) = @_;
	return (-s $fil) if (-e $fil);
	return ((-s "$fil.gz")*5) if (-e "$fil.gz");
	return 0;
}

sub prefixFAhd{
	my ($hr,$nm)= @_;
	my %FNA=%{$hr};
	my %rFNA;
	foreach my $k (keys %FNA){
		my $n;#=$k;
		if ($k =~ m/^>(.*)/){
			$n = $nm.".".$1;
		} else {
			$n = $nm.".".$k;
		}
		$rFNA{$n} = $FNA{$k};
	}
	return(\%rFNA);
}

		
sub gzipwrite{
	my ($outF,$descr) = @_;
	$outF .= ".gz" if ( $outF !~ m/\.gz$/);
	open (my $O, "| gzip -c > $outF") or die "error starting gzip pipe $outF\n$!";
	#my $pigzBin = getProgPaths("piz");
	#open (my $O, "| $pigzBin -c > $outF") or die "error starting gzip pipe $outF\n$!";
	return $O;
}
sub gzipopen{
	my ($inF,$descr) = @_;
	my $dodie = 1;
	if (@_ > 2){$dodie = $_[2];}
	my $verbose=1;
	if (@_ > 3){$verbose = $_[3];}
	$inF .= ".gz" if (!-e $inF && -e $inF.".gz");
	#die "$inF";
	if ($inF =~ m/\.gz$/){
		my $inFwo = $inF; $inFwo =~ s/\.gz$//;
		$inF = $inFwo if (-e $inFwo && !-e $inF);
	} else {
		$inF .= ".gz" if (!-e $inF && -e $inF.".gz");
	}
	
	my $ISTR; my $OK = 1;
	my $msg = "Can't open $descr file $inF\n";
	#print "$dodie  $verbose  $inF\n";
	#if (!-e $inF){{if ($dodie){die $msg;} else { $OK=0;print $msg if ($verbose);}}}
	#my $pigzBin = getProgPaths("pigz");

	if($inF =~ m/\.gz$/ ){
		$msg = "Can't open a pipe to $descr file $inF\n";
		if (!-e $inF) {$OK=0; if ($dodie){die $msg;} else { print $msg if ($verbose);}
		} else {
			#if (!open($ISTR, "$pigzBin -dc $inF |")) {if ($dodie){die $msg;} else {$OK=0; print $msg if ($verbose);}}
			if (!open($ISTR, "gunzip -c $inF |")) {if ($dodie){die $msg;} else {$OK=0; print $msg if ($verbose);}}
		}
	} else{
		if (!open($ISTR, "<", "$inF") ) {if ($dodie){die $msg;} else {$OK=0; print $msg if ($verbose);}}
	}
	#print "$OK $inF\n";
	return ($ISTR,$OK);
}




sub first_index (&@) {
    my $f = shift;
    for my $i (0 .. $#_) {
	local *_ = \$_[$i];	
	return $i if $f->();
    }
    return -1;
}

#makes sure really all previous split fasta files are removed
sub clenSplitFastas($ $){
	my ($inF , $path) = @_;
	$inF =~ m/\/([^\/]+)$/;
	my $inF2 = $1;
	system "rm -f $path/$inF2.*";
}
sub splitFastas($ $ $){
	my ($inF,$num , $path) = @_;
	system "mkdir -p $path" unless (-d $path);
	my $fCnt = 0; my $curCnt=0;
	my $preNum = $num; 
	#check if a file size is given instead of splits..
	if ($preNum =~ m/[GM]$/){
		my $fileSM =  (-s $inF) / (1024 * 1024);
		$num = $preNum; $num =~ s/[MG]//; 
		$num *= 1024 if ($preNum =~ m/G$/);
		$num = int(1+$fileSM / $num);
	}
	$inF =~ m/\/([^\/]+)$/;
	my $inF2 = $1;
	my @nFiles = ("$path/$inF2.$fCnt.$num");
	#print "$nFiles[-1]\n";
	if ($num < 2){
		print "No split required!\n";
		system "rm -f $nFiles[-1];ln -s  $inF $nFiles[-1]";
		return \@nFiles;
	}
	if (-e $nFiles[-1] && -e "$path/$inF2.".($num-1).".$num" && !-e "$path/$inF2.$num.$num"){
		print "seems to exist already\n";
		for (my $i=1;$i<$num;$i++){
			push(@nFiles,"$path/$inF2.$i.$num");
		}
		return \@nFiles;
	}
	system "rm $nFiles[-1]" if (-e $nFiles[-1]);
	my $protN = `grep -c '^>' $inF`;chomp $protN;
	my $pPerFile = int($protN/$num)+10;
	open I,"<$inF"; 
	open my $out,">".$nFiles[-1];
	while (my $l = <I>){
		if ($l =~ m/^>/){
			$curCnt++;
			if ($curCnt > $pPerFile){
				$fCnt++; close $out; 
				push(@nFiles,"$path/$inF2.$fCnt.$num");
				open $out,">$nFiles[-1]";
				$curCnt=0;
			}
		}
		print $out $l;
	}
	close I; close $out;
	return \@nFiles;
	
}


sub readFasta{
	my $fils = $_[0];
	my $cutHd=0;
	my $sepChr= "\\s";
	my %subs; my $doSubs=0;
	$cutHd = $_[1] if (@_ > 1);
	$sepChr = $_[2] if (@_ > 2);
	if (@_ > 3){
		my $hr = $_[3]; %subs = %{$hr};
		$doSubs = 1;
		
	}
	my $Hseq = {}; 
	
	my @files = glob $fils;
	#if (-z $fil){ return \%Hseq;}
	foreach my $fil (@files){
		#next unless (-e $fil);
		#my $FAS; my $status=0; 
		my $doAdd = 1;
		#open($FAS,"<","$fil") || die ("Couldn't open FASTA file $fil\n");
		my ($FAS ,$status) = gzipopen($fil,"fasta file to readFasta",0);
		if (@files == 1 && $status == 0){die "Can't open fasta file $fil\n";}
		next if ($status==0);
		my $temp; my $line; 
		my $trHe =<$FAS>;  
		if (!defined $trHe){#could be empty file
			print "Empty:: $fil $status\n";
			close $FAS;return $Hseq;
		}
		$trHe = substr($trHe,1);
		if ($cutHd) {$trHe =~ s/$sepChr.*//;} chomp ($trHe); $trHe = "".$trHe;
		while($line = <$FAS>){
			chomp($line);
			if ($line =~ m/^>/){
				#check if fasta within set..
				if ($doSubs ){ 
					if ( exists($subs{$trHe})){$doAdd = 1;} else {$doAdd =0;} 
					$Hseq->{$trHe} = $temp if ($doAdd); 
				} else { #a lot faster..
					#finish old fas`
					$Hseq->{$trHe} = $temp;
				}				
				#prep new entry
				$trHe = substr($line,1);# $trHe =~ s/;size=\d+;.*//; 
				if ($cutHd) {$trHe =~ s/$sepChr.*//g;}
				chomp $trHe;
				$trHe = "".$trHe; #just to ensure it's a string
				$temp = "";
				#die $trHe."\n$sepChr\n";
				next;
			}
			$temp .= ($line);
		}
		$Hseq->{$trHe} = $temp;
		close ($FAS);
	}
	return $Hseq;
}


sub median
{
    my @vals = sort {$a <=> $b} @_;
	return 0 if (@vals == 0);
    my $len = (scalar @vals)-1;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub quantile #format: quantile(0.25,@values);
{
	my $cut = shift @_;
	if ($cut<0 || $cut > 1){die "GenoMetaAss::quantile: requested quantile <0 or >1: $cut\nAborting..";}
	if (@_ == 0){return;}
	if (@_ == 1){return ($_[0]);}
	#print "X@_\n\n" if (!defined($_[0]));
    my @vals = sort {$a <=> $b} @_;
    my $len = (scalar @vals)-1;
    return $vals[int(($len*$cut) + 0.5)];
}
sub mean{
    my @vals = @_;
    my $len = scalar @vals;
	return 0 if (@vals == 0);
	my $sum =0; foreach my $a (@vals){$sum+=$a;}
	return ($sum/$len);
}

sub convertNT2AA($){
	my ($text) = @_;
	#print $text."\n";
	 # translate a DNA 3-character codon to an amino acid. We take three letter groups from
	# the incoming string of C A T and G and translate them via a Hash.  It's listed vertically
	# so we can label each of the amino acids as we set it up.

	#$text = "aaatgaccgatcagctacgatcagctataaaaaccccggagctacgatcatcg";

	$text =~ s/[N]*$//i;
	if (length($text) % 3 != 0){
		print "Input not correct length!\n";
		my $ltr = length($text) % 3;
		my $lt = length($text);
		$text = substr $text,0,($lt -$ltr);
	}



	my %convertor = (
		'TCA' => 'S',    # Serine
		'TCC' => 'S',    # Serine
		'TCG' => 'S',    # Serine
		'TCT' => 'S',    # Serine
		'TTC' => 'F',    # Phenylalanine
		'TTT' => 'F',    # Phenylalanine
		'TTA' => 'L',    # Leucine
		'TTG' => 'L',    # Leucine
		'TAC' => 'Y',    # Tyrosine
		'TAT' => 'Y',    # Tyrosine
		'TAA' => '*',    # Stop
		'TAG' => '*',    # Stop
		'TGC' => 'C',    # Cysteine
		'TGT' => 'C',    # Cysteine
		'TGA' => '*',    # Stop
		'TGG' => 'W',    # Tryptophan
		'CTA' => 'L',    # Leucine
		'CTC' => 'L',    # Leucine
		'CTG' => 'L',    # Leucine
		'CTT' => 'L',    # Leucine
		'CCA' => 'P',    # Proline
		'CCC' => 'P',    # Proline
		'CCG' => 'P',    # Proline
		'CCT' => 'P',    # Proline
		'CAC' => 'H',    # Histidine
		'CAT' => 'H',    # Histidine
		'CAA' => 'Q',    # Glutamine
		'CAG' => 'Q',    # Glutamine
		'CGA' => 'R',    # Arginine
		'CGC' => 'R',    # Arginine
		'CGG' => 'R',    # Arginine
		'CGT' => 'R',    # Arginine
		'ATA' => 'I',    # Isoleucine
		'ATC' => 'I',    # Isoleucine
		'ATT' => 'I',    # Isoleucine
		'ATG' => 'M',    # Methionine
		'ACA' => 'T',    # Threonine
		'ACC' => 'T',    # Threonine
		'ACG' => 'T',    # Threonine
		'ACT' => 'T',    # Threonine
		'AAC' => 'N',    # Asparagine
		'AAT' => 'N',    # Asparagine
		'AAA' => 'K',    # Lysine
		'AAG' => 'K',    # Lysine
		'AGC' => 'S',    # Serine
		'AGT' => 'S',    # Serine
		'AGA' => 'R',    # Arginine
		'AGG' => 'R',    # Arginine
		'GTA' => 'V',    # Valine
		'GTC' => 'V',    # Valine
		'GTG' => 'V',    # Valine
		'GTT' => 'V',    # Valine
		'GCA' => 'A',    # Alanine
		'GCC' => 'A',    # Alanine
		'GCG' => 'A',    # Alanine
		'GCT' => 'A',    # Alanine
		'GAC' => 'D',    # Aspartic Acid
		'GAT' => 'D',    # Aspartic Acid
		'GAA' => 'E',    # Glutamic Acid
		'GAG' => 'E',    # Glutamic Acid
		'GGA' => 'G',    # Glycine
		'GGC' => 'G',    # Glycine
		'GGG' => 'G',    # Glycine
		'GGT' => 'G',    # Glycine
		);

	# We don't actually know where the groups of 3 will start in our sample piece of DNA, so we've got
	# three ways of doing the coding ... here's a loop to work out each of the possibilities in turn,
	# leaving the odd extra letters on the beginning or end.
	#for ($s=0; $s<3; $s++) {
	#last; #bs, just need on code
	#        $scrap = substr($text,0,$s);
	#        $main = substr($text,$s);
	#        $main =~ s/(...)/"$convertor{uc $1}" || "?"/eg;
	#        print "$scrap$main\n";
	#        }
	$text =~ s/(...)/"$convertor{uc $1}" || "?"/eg;
	$text =~ s/[^ACDEFGHIKLMNPQRSTVWY*]/X/g;
	#die $text;
	return $text;
}



sub lcp {#
	(join("\0", @_) =~ /^ ([^\0]*) [^\0]* (?:\0 \1 [^\0]*)* $/sx)[0];
}
sub prefix_find($){
	my ($ar) = @_;
	my @rds = @{$ar}; my @newRds = @rds;
	my $first = $rds[0]; 
	#die "$first FFF\n";
	for (my $i=0;$i<@rds; $i++){

		my $second = $rds[$i];
		my @matches;
		next unless (defined $first);
		for (my $start = 0; $start < length($first); $start++) {
			for (my $len = $start+1; $len< length($first); $len++) {
				my $substr = substr($first, $start, $len);
				push @matches, $second =~ m[($substr)]g;
			}
		}
		#print "@matches\n";
		my ($len, $longest) = 0;
		length > $len and ($longest, $len) = ($_, length) for @matches;
		$first = $longest;
		#print "$longest $i\n";
		#"$first\0$second" =~ m/^(.*)\0\1/s;
		#$first = $1;
	}
	
	#print "$first common prefix\n";
	#for (my $i=0; $i<@newRds; $i++){
	#	$newRds[$i] =~ s/^$first//;
	#	$newRds[$i] .= "/" if ($newRds[$i] !~ m/\/$/ && length($newRds[$i] ) != 0 );
	#}
	return ($first);
}

 
sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}



sub readNCBItax($){
	#reads in Jaime's ete tax file
	my ($tIn)  = @_;
	open I,"<$tIn";
	while (my $l = <I>){
		chomp $l;
		my @spl = split (/\t/,$l);
		my $txID = $spl[0];
		my $rnksLvl = $spl[2];
		my $rnks = $spl[3];
	}
	close I;
}

sub readFastHD{#only reads headers
	my $inF = $_[0];
	my $complH=0;
	if (@_ > 1){$complH = $_[1];}
	#die "$complH\n";
	my @ret;
	open I,"<$inF" or die "can t open $inF\n";
	while (my $l = <I>){
		chomp $l;
		if ($l =~ m/^>(\S+)/){
			if (!$complH){
				push @ret,$1;
			} else {
				$l =~ m/^>(.*)/;
				push @ret,$1;
			}
		}
	}
	close I;
	return \@ret;
}

sub renameFastHD($ $ $){ #set a new name for headers in fasta files
	my ($inF,$hr,$extr) = @_;
	my %COGid = %{$hr};
	my $oS = "";
	my %CidCnt;
	open I,"<$inF" or die "can t open renameFastHD  $inF\n";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			chomp $l;
			unless (exists $COGid{$l}){die "Can't find $l in COGs\n";}
			if (exists $CidCnt{$COGid{$l}}){$CidCnt{$COGid{$l}}++;
				$oS .= ">".$extr."_".$COGid{$l}."..".$CidCnt{$COGid{$l}}."\n";
			} else {$CidCnt{$COGid{$l}} = 0;
				$oS .= ">".$extr."_".$COGid{$l}."\n";
			}
			

		} else {
			$oS .= $l;
		}
	}
	close I;
	open O,">$inF" or die "Can't open rename Fasta HD out file $inF\n";
	print O $oS;
	close O;
}
sub renameFastqCnts($ $){ #set a new name for headers in fastq files, using a simple scheme of prefix and just counts afterwards
	my ($inF,$prefix) = @_;
	open I,"<$inF" or die "can t open $inF\n";
	my $cnt  = 0; my $cnt2 = 0;
	open O,">$inF.tmp"; my $plusSeen = 1;
	while (my $l = <I>){
		if ($l =~ m/^@/ && $plusSeen && $cnt2 >= 4 ){ #$l =~ m/^@/ && 
			print O "@".$prefix."_".$cnt."\n";
			$cnt ++; $plusSeen=0;$cnt2=0;
#			print L ">".$prefix."_".$cnt."\t$l";
		} else {
			print O $l;
			$plusSeen=1 if ($l =~ m/^\+\n$/);
		}
		$cnt2 ++;
	}
	close I; close O;
	system "rm $inF;mv $inF.tmp $inF";
}

sub renameFastaCnts($ $ $){ #set a new name for headers in fasta files, using a simple scheme of prefix and just counts afterwards
	my ($inF,$prefix,$logF) = @_;
	open I,"<$inF" or die "can t open $inF\n";
	my $cnt  = 0;
	open O,">$inF.tmp";
	open L,">$logF";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			print O ">".$prefix."_".$cnt."\n";
			print L ">".$prefix."_".$cnt."\t$l";
		} else {
			print O $l;
		}
		$cnt ++;
	}
	close I; close O; close L;
	system "rm $inF;mv $inF.tmp $inF";
}

sub readClstrRevSmplCtgGenSubset{
	my $inF = $_[0];
	my $subsHR = {};
	$subsHR = $_[1] if (@_ > 1);
	#my %subs = %{$subsHR};
	my $retR={}; my %retF;
	
	#my @key = keys %{$subsHR}; print "$key[0] $key[1]\n";
	print "Reading clstr $inF  .. \n";
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=""; my $cnts=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = ($arr[0]);#substr($_,0,$pos);
		
		#extract reverse: contig_gene to gene cat 
		my @tmpArr = split /,/,$arr[1];
		foreach my $gene (@tmpArr) {
			my $ctg = substr($gene,1); $ctg =~ s/_(\d+)$//; my $geneN=$1;
			#print "$ctg = $gene = $geneN\n";
			 if (exists($$subsHR{$ctg})){
				 my @spl2 = split(/__/,$ctg);
				#${$retR}{$ctg}{$geneN} = $curCl;
				${$retR}{$spl2[0]}{$spl2[1]}{$geneN} = $curCl;
				$cnts++;
			 }
		}
	}
	close $I;
	print "Found $cnts genes.\n";
	#die;
	return ($retR);
}

sub readClstrRevContigSubset{ #gets the exact assembled genes clustered in GC genes, subset by contigs (used for clusterMAGs.pl)
	my $inF = $_[0];
	my $subsHR = {};
	$subsHR = $_[1] if (@_ > 1);
	#my %subs = %{$subsHR};
	my $retR={}; my %retF;
	
	#my @key = keys %{$subsHR}; print "$key[0] $key[1]\n";
	print "Reading clstr $inF  .. ";
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=""; my $cnts=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = ($arr[0]);#substr($_,0,$pos);
		
		#extract reverse: contig_gene to gene cat 
		my @tmpArr = split /,/,$arr[1];
		foreach my $gene (@tmpArr) {
			my $ctg = substr($gene,1); $ctg =~ s/_\d+$//; 
			#print "$ctg = $gene\n";
			 if (exists($$subsHR{$ctg})){
				${$retR}{$gene} = $curCl;
				$cnts++;
			 }
		}
	}
	close $I;
	print "Found $cnts genes.\n";
	#die;
	return ($retR);
}


sub readClstrRev{ #gets the exact assembled genes clustered in GC genes #version for my shortened index file
	my $inF = $_[0];
	my $createR = 1;
	$createR = $_[1] if (@_ > 1);
	my $subsHR = {};
	$subsHR = $_[2] if (@_ > 2);
	#my %subs = %{$subsHR};
	my $doSubset=0;
	if (scalar(keys(%{$subsHR})) > 0 ){
		$doSubset = 1 ;
		if ($createR == 2){
			$doSubset = 2 ;
		}
	}
	my $retR={}; my %retF;
	print "Reading clstr $inF  .. \n";
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = ($arr[0]);#substr($_,0,$pos);
		
		if ($doSubset == 1){
			next unless (exists($$subsHR{$curCl}));
		}
		
		#my $rem = $arr[1]; #substr($_,$pos+1);
		#print $curCl."XX$rem\n";
		#foreach (split /,/,$rem){$retR{$_} = $curCl;}
		if ($createR > 0){
			my @tmpArr = split /,/,$arr[1];
			if ($doSubset == 1){
				my $gogo=0;
				foreach (@tmpArr) {
					 if (exists($$subsHR{$_})){
						 $gogo=1;
						 last;
					 }
				}
				next unless ($gogo);
			}
			@{$retR}{@tmpArr} = ($curCl) x @tmpArr;
		}
		if ($createR != 2){
			$retF{$curCl} = $arr[1];
		}
		@arr = ();

	}
	close $I;
	print "done reading ClStr\n";
	return ($retR,\%retF);
}


sub readClstrRevGenes{ #gets the exact assembled genes clustered in GC genes #version for my shortened index file
	my $inF = $_[0];
	my $subsHR = {};
	$subsHR = $_[1] if (@_ > 1);
	#my %subs = %{$subsHR};
	my $doSubset=0;
	$doSubset = 1 if (scalar(keys(%{$subsHR})) > 0 );
	my %retR; 
	print "Reading clstr2 $inF  .. \n";
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = int($arr[0]);#substr($_,0,$pos);
		
		if ($doSubset){
			next unless (exists($subsHR->{$curCl}));
		}
		my @tmpArr = split /,/,$arr[1];
		foreach (@tmpArr) {
			next unless (exists($subsHR->{$_}));
			$retR{$_} = $curCl;
		}
#			@retR{@tmpArr} = ($curCl) x @tmpArr;
	}
	close $I;
	print "done reading ClStr gene 2\n";
	return (\%retR);
}



sub unzipFileARezip($){
	my ($inFar) = @_;
	my @inFs = @{$inFar}; my $totCnt=0;
	my $bef = "gunzip "; my $aft="gzip ";
	#print "$inF.gz\n";
	foreach my $inF (@inFs){
		if (-e "$inF.gz"){
			$bef .= " $inF.gz"; $aft .= " $inF"; $totCnt++;
		}
	}
	if ($totCnt==0) {
		return ("","");
	} else {
		return $bef."\n",$aft."\n";
	}
}
sub reverse_complement_IUPAC ($) {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}


sub readGFF($){
	my ($inF) =@_;
	my %ret;my $sbcnt=1;
	my $lcnt=0;
	open I,"<$inF";
	while (<I>){
		if (m/^#/){$sbcnt=1;next;}
		chomp;
		$lcnt++;
		my @spl = split(/\t/);
		print "readGFF::too short: line $lcnt, $inF\n" unless (@spl > 7);
		$spl[8] =~ m/^ID=\d+_(\d+)/;
		my $k = ">".$spl[0]."_$1";
		#print $k;
		$ret{$k}=$_;
	}
	#die;
	close I;
	return \%ret;
}


sub getAssemblPath{
	my $cD = $_[0];
	my $newpath = "";
	$newpath = $_[1] if (@_ > 1);
	$newpath =~ s/\/\//\//g;$newpath =~ s/\/\//\//g;
	my $dieOnFail=1;
	$dieOnFail = $_[2] if (@_ > 2);
	my $firstF = "$cD/assemblies/metag/assembly.txt";
	if (!-e $firstF){
		print "GenoMetaAss::getAssemblPath::Assembly path missing: $cD\n";
		return ("") if (!$dieOnFail);
		
	}
	open my $I,"<$firstF" or die "GenoMetaAss::getAssemblPath::Assembly path missing: $cD\n";
	my $metaGD = <$I>;#= `cat $cD/assemblies/metag/assembly.txt`; 
	chomp $metaGD;
	my $metaGD2 = $metaGD;$metaGD2 =~ s/\/\//\//g;$metaGD2 =~ s/\/\//\//g;
	
	close $I;
	if ($newpath ne ""){
		if ($newpath ne $metaGD2 || !-d $metaGD  ){#needs replacement?
			print "replacing $metaGD with $newpath\n" ;
			open my $I,">$firstF"; print $I $newpath; close $I;
			$metaGD = $newpath;
		}
	}
	
	return($metaGD);
}


sub iniCleanSeqSetHR{
	my ($seqSetHR) = @_;
	my $HR = {arp1 => ${$seqSetHR}{pa1},arp2 => ${$seqSetHR}{"pa2"},singAr => ${$seqSetHR}{"pas"}, matAr => [],
					readTec => ${$seqSetHR}{seqTech}, is3rdGen => ${$seqSetHR}{is3rdGen},
					arpX1 => ${$seqSetHR}{paX1},arpX2 => ${$seqSetHR}{paX2},singArX => ${$seqSetHR}{paXs}, matArX => [],
					readTecX => ${$seqSetHR}{seqTechX}, is3rdGenX => ${$seqSetHR}{is3rdGenX}
				};
	return $HR;
}

sub addFileLocs2AssmGrp{
	my ($AsGrpsHR, $cAssGrp,$SmplName, $cleanSeqSetHR, $seqSetHR) = @_;
	${${$AsGrpsHR}{$cAssGrp}{CleanSeqs}}{$SmplName} = $cleanSeqSetHR;
	${${$AsGrpsHR}{$cAssGrp}{RawSeqs}}{$SmplName} = $seqSetHR;
	return $AsGrpsHR;
}


sub getCleanSeqsAssmGrp{
	my ($asG, $grp, $support) = @_;
	my $specSmpl = "";$specSmpl = $_[3] if (@_ >= 4); #specific sample only??
	my @smpls = ();@smpls = ($specSmpl) if ($specSmpl ne "");

	my @pa1; my @pa2; my @pas; my @readTec;
	my %raws = %{${$asG}{$grp}{CleanSeqs}};
	my @terms = ("arp1","arp2","singAr", "readTec","is3rdGen");
	@terms = ("arpX1","arpX2","singArX", "readTecX","is3rdGenX") if ($support);
	@smpls = keys %raws if (!@smpls);
	foreach my $smpl (@smpls){
	#	print "$smpl\n";
		push(@pa1, @{${$raws{$smpl}}{ $terms[0] }});
		push(@pa2, @{${$raws{$smpl}}{ $terms[1] }});
		push(@pas, @{${$raws{$smpl}}{ $terms[2] }});
		push(@readTec, ${$raws{$smpl}}{ $terms[3] });
	}
	#die "@pa1 - $grp - $support\n@pa2\n@pas\n";
	return (\@pa1, \@pa2, \@pas, \@readTec);
}

sub hasSuppRds{
	my ($asG, $grp,$smpl) = @_;
	my %raws = %{${$asG}{$grp}{RawSeqs}};
	my @smpls = keys %{${$asG}{$grp}{RawSeqs}};
	print "GRP:$grp\n@smpls\n";
#	if (  exists( ${  ${${$asG}{$grp}{RawSeqs}}{$smpl}}{"paXs"} )  ){
	if (  exists( ${  $raws{$smpl}}{"paXs"} )  ){
		return 1;
	}
	return 0;
}

sub getRawSeqsAssmGrp{
	my ($asG, $grp, $support) = @_;
	
	my $specSmpl = "";$specSmpl = $_[3] if (@_ >= 4); #specific sample only??
	my @smpls = ();@smpls = ($specSmpl) if ($specSmpl ne "");
	
	my @pa1; my @pa2; my @pas; my @libs;my @readTec;
	my %raws = %{${$asG}{$grp}{RawSeqs}};
	my @terms = ("pa1","pa2","pas", "libInfo","seqTech","is3rdGen");
	@terms = ("paX1","paX2","paXs", "libInfoX","seqTechX","is3rdGenX") if ($support);
	@smpls = keys %raws if (!@smpls);
	#die "XX @smpls YY\n$support\n@terms\n@{${$raws{$smpls[0]}}{ $terms[2] }}\n";
	foreach my $smpl (@smpls){
		#print "$smpl\n"; 
		push(@pa1, @{${$raws{$smpl}}{ $terms[0] }});
		push(@pa2, @{${$raws{$smpl}}{ $terms[1] }});
		push(@pas, @{${$raws{$smpl}}{ $terms[2] }});
		push(@libs, @{${$raws{$smpl}}{ $terms[3] }});
		push(@readTec, ${$raws{$smpl}}{ $terms[4] });
	}
	#die "@pa1\n@pa2\n@pas\n";
	return (\@pa1, \@pa2, \@pas, \@libs, \@readTec);
}


sub emptyAssGrpsObj($){
	my ($hr) = @_;
	my %AsGrps = %{$hr};
	foreach my $k (keys %AsGrps){
		$AsGrps{$k}{CntAss} = 0;
		$AsGrps{$k}{CntMap} = 0;
		#dependencies
		$AsGrps{$k}{ITSDeps} = "";
		$AsGrps{$k}{MapDeps} = "";
		$AsGrps{$k}{DiamDeps} = "";
		$AsGrps{$k}{SeqClnDeps} = "";
		#print $AsGrps{$k}{CntAimAss} ."\n";
		#remove tmp dirs
		$AsGrps{$k}{ClSeqsRm} = "";
		#copy to final folder
		@{$AsGrps{$k}{MapCopies}} = ();
		@{$AsGrps{$k}{AssCopies}} = ();
		$AsGrps{$k}{MapCopiesNoDel} = [];
		#filteredSequenceFiles for assembly
		$AsGrps{$k}{PostAssemblCmd} = "";
		$AsGrps{$k}{PostClnCmd} = "";
		$AsGrps{$k}{AssemblSmplDirs} = "";
		$AsGrps{$k}{scndMapping} = "";
		$AsGrps{$k}{ClSeqsRm} = "";
		#complex hashes replaces FilterSeq1 ..
		$AsGrps{$k}{CleanSeqs} = {};
		$AsGrps{$k}{RawSeqs} = {};
		#@{$AsGrps{$k}{FilterSeq1}} = (); @{$AsGrps{$k}{FilterSeq2}} = (); @{$AsGrps{$k}{FilterSeqS}} = ();
		#@{$AsGrps{$k}{RawSeq1}} = (); @{$AsGrps{$k}{RawSeq2}} = (); @{$AsGrps{$k}{Libs}} = ();
		
		#print $k."\n";;
	}
	$AsGrps{global}{DiamCln} = "";
	$AsGrps{global}{DiamDeps} = "";
	#die;
	return \%AsGrps;
}



sub readMapS{
	my ($inF,$folderStrClassical) = ($_[0],$_[1]);
	my $xtraCols = ""; $xtraCols = $_[2] if (@_ > 2);
	my @spl = split /,/,$inF;
	#my %ret; my %agbp;
	my @outDirs = ();
	#my $hr1 = \%ret;my $hr2 = \%agbp;  my $cnt = -1;
	my $hr1 = {};my $hr2 = {};  my $cnt = -1;
	foreach my $map (@spl){
		($hr1,$hr2) = readMap($map,$cnt,$hr1,$hr2,$folderStrClassical,$xtraCols);
		#%ret = %{$hr1};
		$cnt = $hr1->{totSmpls};
		push(@outDirs,$hr1->{opt}{outDir});
		#print $cnt."\n";
	}
	#%ret = %{$hr1};
	#print keys %ret;
	$hr1->{opt}{outDir} = join(",",@outDirs);
	#die;
	return ($hr1,$hr2);
}

#infer Assembly dirs & corrsponding bams with several Samples (compound assemblies)
sub getDirsPerAssmblGrp{
	my %AsGrps; my %map;
	my $chkDirs=0;
	if (@_ == 1){
		my ($mapF) = @_;
		my ($hrm,$asGrpObj) = readMapS($mapF);
		%AsGrps = %{$asGrpObj};		%map = %{$hrm};
		$chkDirs=1;
	} elsif (@_ == 2) {
		my ($hrm,$asGrpObj) = @_;
		%AsGrps = %{$asGrpObj};		%map = %{$hrm};
	} else {
		die "Too many/few args given in GenoMetaAss::getDirsPerAssmblGrp\n";
	}
	my @smpls = @{$map{opt}{smpl_order}};
	my $cnt=0;
	my %DOs;
	foreach my $smpl (@smpls){
		my $cntAim = $AsGrps{ $map{$smpl}{MapGroup} }{CntAimMap};
		my $dir2rd = $map{$smpl}{wrdir};
		my $cAssGrp = $map{$smpl}{AssGroup};
		my $cMapGrp = $map{$smpl}{MapGroup};
		$AsGrps{$cMapGrp}{CntMap} ++;
		#next if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} );

		my $tar = "AssmblGrp_$cAssGrp";
		if ($chkDirs && !-d "$dir2rd"){
			print "Can't read $dir2rd\nSkipping..\n";
			next;
		}
		#(-s "$curOutDir/mapping/Align_ment-smd.bam" || -s "$curOutDir/mapping/Align_ment-smd.cram")
		#my $bam = $inD."mapping/Align_ment-smd.bam";
		#print "$cAssGrp\n";
		push (@{$DOs{$cAssGrp}{wrdir}}, $dir2rd);
		push(@{$DOs{$cAssGrp}{SmplID}},$map{$smpl}{SmplID});
		$cnt++;
		#last if ($cnt > 50);
	}
	return(\%DOs,\%map);
}


#should only be used for dirs!
sub resolve_path($){
	my ($inP) = @_;
	return "" if ($inP eq "");
	$inP =~ s/\$([A-Z0-9_]*)/$ENV{$1}/g;
	#doesn't work: dir might not exist yet!
	#$inP = `realpath $inP` ; chomp $inP; $inP .= "/";
	#if (-d $inP && $inP !~ m/\/$/){$inP .= "/";}
	$inP =~ s/\/\//\//g;$inP =~ s/\/\//\//g;
	return $inP;
}

sub is3rdGenSeqTech{
	my $curReadTec = $_[0];
	my $is3rdGen = 0; #important flag for long reads (Oxford Nanopore / PacBio)
	$is3rdGen = 1 if ($curReadTec eq "ONT" || $curReadTec eq "PB");
	return $is3rdGen;
}

sub checkSeqTech{
	#ONT,PB,proto,miSeq,GAII
	my $inT = $_[0];
	my $msg = "Mapping file";
	$msg = $_[1] if (@_ > 1);
	if ($inT ne "" && $inT ne "ONT"&& $inT ne "hiSeq" && $inT ne "454" && $inT ne "SLR" && $inT ne "PB" && $inT ne "proto" && $inT ne "miSeq" && $inT ne "GAII" && $inT ne "GAII_solexa"){
		die "$msg: Can't recognize SeqTech: \"$inT\"\nHas to be one of \"ONT\",\"PB\",\"proto\",\"SLR\",\"miSeq\",\"hiSeq\",\"GAII\",\"GAII_solexa\" or \"\"\n\n";
	}
}
sub readMap{
	my $inF = $_[0];
	my $Scnt = defined $_[1] ? $_[1] : 0;
	my %ret = defined $_[2] ? %{$_[2]} : (); 
	my %agBP = defined $_[3] ? %{$_[3]} : ();
	my $folderStrClassical = defined $_[4] ? $_[4] : -1;
	my $xtraColStr = defined $_[5] ? $_[5] : "";
#die "$folderStrClassical\n";
	my @order = exists $ret{opt}{smpl_order} ? @{$ret{opt}{smpl_order}} : ();
	my %oldAssmGrps = ();
	if (@order > 0){#something already exists
		my @asgps= keys (%agBP);
		$oldAssmGrps{$_}++ for (@asgps);
	}
	my $dirCol = -1; my $SmplPrefixCol = -1; #to find location of primary reads either has to be defined.
	my $smplCol = 0; my $rLenCol = -1; 
	my $SeqTech = -1; #mate = mate pairs; SLR = synthetic long reads, e.g. 10X, TELseq..
	my $SeqTechS = -1; 
	my $AssGroupCol = -1; my $EstCovCol = -1; my $MapGroupCol = -1; my $SupRdsCol = -1;my $ExcludeAssemble = -1;
	my $cut5pR1 = -1;my $cut5pR2 = -1; my $FamGroupCol = -1;
	#some global params
	my $dir2dirs = ""; #dir on file system, where all dirs specified in map can be found (enables different indirs with different mapping files)
	my $dir2out = "";
	my $baseID = ""; my $mocatFiltPath = "";
	my $inDirSet = 0;
	my $cnt = -1;
	my $illuminaClip ="";
	my $GlbTmpD = "";
	my $NodeTmpD = "";
	my $infFoldClass = -1;
	my $xtraCol = -1;
	my $relaxSmplID = 0;
	my @dir2dirsA;
	my %trackMGs; my %trackAGs; #hashes to track the last (final) sample in each mapgroup.. important to know this to check if assembly / mapping is done 
	my %memberMGs ; my %memberAGs ;
	my %trackDirs;
	my %trackPrefixs;
	
	my $DOWARN = 1;
	my $warnDeactivateMsg = "In case you want to continue, insert \"#WARNING OFF\" underneath the header of your map file.\n";
	#print $inF."\n";
	open I,"<$inF" or die "Can't open map: $inF\n";
	#AssGrps
	while (<I>){
		#use chomp and s/\R//g; to catch also windows formated lines..
		$cnt++;s/\R//g;chomp;
		next if (length($_) ==0);
		if (m/^#/ && $cnt > 0 ){#check for ssome global parameters
			if (m/^#DirPath\s(\S+).*$/){$dir2dirs = resolve_path($1); $dir2dirs.="/" unless ($dir2dirs=~m/\/$/); push(@dir2dirsA,$dir2dirs);}
			if (m/^#OutPath\s+(\S+)/){$dir2out = resolve_path($1); $inDirSet=0;$dir2out .= "/" if ($dir2out !~ m/\/$/);}
			if (m/^#RunID\s+(\S+)/){$baseID = $1; $inDirSet=0;}
			if (m/^#mocatFiltPath\s+(\S+)/){$mocatFiltPath = resolve_path($1);}
			if (m/^#illuminaClip\s+(\S+)/){$illuminaClip = resolve_path($1);}
			if (m/^#NodeTmpDir\s+(\S+)/){$NodeTmpD = $1;}
			if (m/^#GlobalTmpDir\s+(\S+)/){$GlbTmpD = $1;}
			if (m/^#WARNING\sOFF/){$DOWARN=0;print "Warning: Deactivated Warnings while reading the map file! ..\n";}
			if (m/^#RelaxSMPLID\sTRUE/){$relaxSmplID = 1;print "Relaxed SMPLIDs\n";}
			if (!$inDirSet && $dir2out ne "" && $baseID ne ""){
				$dir2out.=$baseID unless ($dir2out =~ m/$baseID[\/]*$/); $inDirSet =1;
				$dir2out .= "/" if ($dir2out !~ m/\/$/);
			}
			$dir2out .= "/" if ($dir2out !~ m/\/$/);
			next;
		}
		my @spl = split(/\t/,$_,-1);
		if ($cnt == 0){
			#die "@spl\n";
			$smplCol = first_index { /^#SmplID$/ } @spl;
			$dirCol = first_index { /^Path$/ } @spl;
			$SmplPrefixCol = first_index { /^SmplPrefix$/ } @spl;
			$SeqTech = first_index { /^SeqTech$/ } @spl;
			$SeqTechS = first_index { /^SeqTechSingl$/ } @spl;
			$rLenCol = first_index { /^ReadLength$/ } @spl;
			$AssGroupCol = first_index { /^AssmblGrps$/ } @spl;
			$FamGroupCol = first_index { /^FamilyGrps$/ } @spl;
			$EstCovCol = first_index { /^EstCoverage$/ } @spl;
			$MapGroupCol = first_index { /^MapGrps$/ } @spl;
			$SupRdsCol = first_index { /^SupportReads$/ } @spl;
			$ExcludeAssemble = first_index { /^ExcludeAssembly$/ } @spl;
			$cut5pR2 = first_index { /^cut5PR2$/ } @spl;
			$cut5pR1 = first_index { /^cut5PR1$/ } @spl;
			if ($xtraColStr ne ""){
				$xtraCol = first_index { /$xtraColStr/ } @spl;
			}
			#die "MAP: $AssGroupCol\n";
			#die "Only \"Path\" or \"SmplPrefix\" can be defined in mapping file. Both is not supported.\n" if ($dirCol != -1 && $SmplPrefixCol != -1);
			die "Could not find \"#SmplID\" in input map\n" unless($smplCol> -1);
			die "Expected to find at least \"SmplPrefix\" or \"Path\" as column headers in .map\n" if ($dirCol == -1 && $SmplPrefixCol == -1);
			die "Either \"SmplPrefix\" or \"Path\" has to be second column in .map\n" unless ($dirCol == 1 || $SmplPrefixCol == 1);
			next;
		} #maybe later check for col labels etc
		my $smplPrefixUsed=0; my $samplePathUsed=0;
		$Scnt++;
		#die "1 $GlbTmpD 2 $NodeTmpD map\n";
		#die $spl[0]." ".$spl[1]."\n";
		die "inPath has to be set in mapping file!\n" if ($dir2dirs eq "");
		#die "$dir2out\n";
		my $curSmp = $spl[$smplCol];
		my $altCurSmp = "";
		#print $curSmp." ";
		die "\"opt\" is reserved keyword and can't be used as samples\n" if ($curSmp eq "opt");
		die "\"altNms\" is reserved keyword and can't be used as samples\n" if ($curSmp eq "altNms");
		#die "\"smpl_order\" is reserved keyword and can't be used as samples\n" if ($curSmp eq "smpl_order");
		#die "\"totSmpls\" is reserved keyword and can't be used as samples\n" if ($curSmp eq "totSmpls");
		die "Double sample ID \"$curSmp\"\n$_\n" if (exists $ret{$curSmp});
		die "Can't use character \"$\" in sampleID: $curSmp\n" if ($curSmp =~ m/\$/);
		die "Can't use character \"_\" in sampleID: $curSmp\n" if ($curSmp =~ m/_/);
		die "Can't use character \",\" in sampleID: $curSmp\n" if ($curSmp =~ m/,/);
		die "Can't use character \"-\" in sampleID: $curSmp\n" if (!$relaxSmplID && $curSmp =~ m/-/);
		die "Recommended not to use numeric as first char \"0-9\" in sampleID: $curSmp\n" if (!$relaxSmplID && $curSmp =~ m/^[0-9]/);
		die "Use alphanumeric characters (a-zA-Z0-9.) for sampleID: $curSmp\n" if (!$relaxSmplID && $curSmp =~ m/[^a-zA-Z0-9\.]/);
		
		my $cdir = ""; 
		if ($dirCol >= 0 && @spl>= $dirCol && $spl[$dirCol] ne ""){
			$cdir = $spl[$dirCol] ; 
			$samplePathUsed=1;
			if ($DOWARN && exists($trackDirs{"$dir2dirs/$cdir"}) ){ die "Warning: Found the sample path \"$cdir\" more than once. This would lead to using reads twice, aborting.\n $warnDeactivateMsg";}
			$trackDirs{"$dir2dirs/$cdir"} = 1;
		}
		my $cdir2= $cdir;
		$ret{$curSmp}{dir} = $cdir;#this one should stay without a tag
		$ret{$curSmp}{rddir} = $dir2dirs.$cdir;
		$ret{$curSmp}{clip} = $illuminaClip;
		$ret{$curSmp}{rddir} .="/" unless ($ret{$curSmp}{rddir} =~ m/\/$/);
		#die "$ret{$curSmp}{rddir} $dirCol $cdir $curSmp\n $smplCol $dirCol\n";
		if ($SmplPrefixCol>=0 && @spl>= $SmplPrefixCol && $spl[$SmplPrefixCol] ne ""){
			$cdir2 = $spl[$SmplPrefixCol];
			$ret{$curSmp}{prefix} = $cdir2;
			if ($DOWARN && exists($trackPrefixs{"$dir2dirs/$cdir2"}) ){ die "Warning: Found the sample path \"$cdir2\" more than once. This would lead to using reads twice, aborting.\n $warnDeactivateMsg";}
			$trackPrefixs{"$dir2dirs/$cdir2"} = 1;
			$smplPrefixUsed =1;
		} else {
			$ret{$curSmp}{prefix} = "";
		}
		$cdir2.="/" unless ($cdir2 =~ m/\/$/);
		
		#basic check that no twice usage..
		if ($DOWARN && $smplPrefixUsed && $samplePathUsed){die"Warning: in mapping file both \"SmplPrefix\" and \"Path\" are set for sample $curSmp!\nThis is not supported\n$warnDeactivateMsg";}

		
		
		if ($folderStrClassical== -1){
			if (-d  $dir2out.$cdir2 && !-d $dir2out.$curSmp){
				if ($infFoldClass==0){die "readMap: Inferring old/new folder structure failed, as both folders seem to be valid\n";}
				$ret{$curSmp}{wrdir} = $dir2out.$cdir2;
				$infFoldClass = 1;
			} else {
				if ($infFoldClass==1){die "readMap: Inferring old/new folder structure failed, as both folders seem to be valid\n";}
				$ret{$curSmp}{wrdir} = $dir2out.$curSmp."/";
				$infFoldClass = 0;
			}
		}elsif ($folderStrClassical == 1){
			$ret{$curSmp}{wrdir} = $dir2out.$cdir2;
		} else {
			$ret{$curSmp}{wrdir} = $dir2out.$curSmp."/";
		}
		#die "$ret{$curSmp}{wrdir}\n";
		
		$ret{$curSmp}{SmplID} = $curSmp;
		$ret{$curSmp}{mapFinSmpl} = $curSmp;
		$ret{$curSmp}{assFinSmpl} = $curSmp;
		#ONT,PB,proto,miSeq,GAII
		if ($SeqTech >= 0) { my $RT=$spl[$SeqTech];checkSeqTech($RT);$ret{$curSmp}{SeqTech} = $RT; } else {$ret{$curSmp}{SeqTech} = "";}
		if ($SeqTechS >= 0) { my $RT=$spl[$SeqTechS];checkSeqTech($RT);$ret{$curSmp}{SeqTechSingl} = $RT; } else {$ret{$curSmp}{SeqTechSingl} = "";}
		
		
		#die "$SupRdsCol\t@spl\n$spl[$SupRdsCol]\n";
		if ($SupRdsCol >= 0 && $SupRdsCol < @spl) { 
			if(length($spl[$SupRdsCol]) > 0 && $spl[$SupRdsCol] !~ m/\/$/ && -d $spl[$SupRdsCol]) {
				$spl[$SupRdsCol].="/";
			}
			#resolve paths.. split on , for multiple files/cases
			my @spl2 = split /;/,$spl[$SupRdsCol];
			my @spl4 = ();
			foreach my $case (@spl2){
				my @spl3 = split /:/,$case;
				if (@spl3>=2){
					my $out = shift @spl3;$out .=  ":";
					foreach (@spl3){$out .=resolve_path($_);}
					push (@spl4,$out);
				} else {
					push (@spl4,$case);
				}
			}
			$ret{$curSmp}{SupportReads} = join(",",@spl4); 
			
			#print "\n\n$ret{$curSmp}{SupportReads} \n\n";
			
		} else {
			$ret{$curSmp}{SupportReads} = "";
		}
		if ($rLenCol >= 0){$ret{$curSmp}{readLength} = $spl[$rLenCol];} else {$ret{$curSmp}{readLength} = 0;}
		if ($ExcludeAssemble >= 0){$ret{$curSmp}{ExcludeAssem} = $spl[$ExcludeAssemble];} else {$ret{$curSmp}{ExcludeAssem} = 0;}
		$ret{$curSmp}{cut5pR2} = 0;
		if ($cut5pR2 >= 0){
			if (defined ($spl[$cut5pR2]) && length($spl[$cut5pR2]) > 0){
				$ret{$curSmp}{cut5pR2} = $spl[$cut5pR2];
			}
		} 
		$ret{$curSmp}{cut5pR1} = 0;
		if ($cut5pR1 >= 0){
			if (defined ($spl[$cut5pR1]) && length($spl[$cut5pR1]) > 0){
				$ret{$curSmp}{cut5pR1} = $spl[$cut5pR1];
			}
		} 
		if ($xtraCol != -1){$ret{$curSmp}{$xtraColStr} = $spl[$xtraCol];}
		if ($EstCovCol >= 0){$ret{$curSmp}{DoEstCoverage} = $spl[$EstCovCol];} else {$ret{$curSmp}{DoEstCoverage} = 0;}
		if ($AssGroupCol >= 0 && $spl[$AssGroupCol] ne ""){
			my $curAG = $spl[$AssGroupCol];
			$ret{$curSmp}{AssGroup} = $curAG ;
			if (!exists($agBP{$curAG}{CntAimAss})){$agBP{$curAG}{CntAimAss}=0;}
			$agBP{$curAG}{CntAimAss}++;
			$altCurSmp = $curSmp."M".$agBP{$curAG}{CntAimAss};
			$memberAGs{$curAG} = [] unless (exists $memberAGs{$curAG});
			$trackAGs{$curAG} = $curSmp;  push(@{$memberAGs{$curAG}},$curSmp);
			$agBP{$curAG}{prodRun} = "";
			if (exists($oldAssmGrps{$curAG})){
				die "Warning: assembly group \"$curAG\" seems to exist in one of the previous maps. This is currently not supported for MATAFILER, please make sure all assembly groups are unique to a single map!\n";
			}
			#print $agBP{$spl[$AssGroupCol]}{CntAimAss}. " :$spl[$AssGroupCol]\n" ;
		} else {$ret{$curSmp}{AssGroup} = $Scnt; $agBP{$Scnt}{CntAimAss}=0;$agBP{$Scnt}{prodRun} = "";}
		if ($MapGroupCol >= 0 && $spl[$MapGroupCol] ne ""){
			my $curM = "M_".$spl[$MapGroupCol];
			$ret{$curSmp}{MapGroup} = $curM;
			if (!exists($agBP{$curM}{CntAimMap})){$agBP{$curM}{CntAimMap}=0;}
			$agBP{$curM}{CntAimMap}++;
			$memberMGs{$curM} = [] unless (exists $memberMGs{$curM});
			$trackMGs{$curM} = $curSmp; push(@{$memberMGs{$curM}},$curSmp);
			#print $agBP{$spl[$MapGroupCol]}{CntAimMap}. " :$spl[$MapGroupCol]\n" ;
		} else {$ret{$curSmp}{MapGroup} = $Scnt; $agBP{$Scnt}{CntAimMap}=0;}
		
		if ($FamGroupCol >= 0 && $spl[$FamGroupCol] ne ""){
			my $curF = $spl[$FamGroupCol];
			$ret{$curSmp}{FamGroup} = $curF;
			#if (!exists($agBP{$curF}{CntAimFam})){$agBP{$curF}{CntAimFam}=0;}
			$agBP{$curF}{CntAimFam}++;
			$memberMGs{$curF} = [] unless (exists $memberMGs{$curF});
			$trackMGs{$curF} = $curSmp; push(@{$memberMGs{$curF}},$curSmp);
		} else {$ret{$curSmp}{FamGroup} = ""; $agBP{$Scnt}{CntAimFam}=0;}
		
		if ($altCurSmp ne ""){$ret{altNms}{$altCurSmp} = $curSmp;} 
		push(@order,$curSmp);
		#print $spl[0]."\n";
	}
	
	#insert final sample destination for all AssGroups and MapGroups
	foreach my $k (keys %memberMGs){
		foreach (@{$memberMGs{$k}}){  $ret{$_}{mapFinSmpl} = $trackMGs{$k};  }
	}
	foreach my $k (keys %memberAGs){
		foreach (@{$memberAGs{$k}}){  
			$ret{$_}{assFinSmpl} = $trackAGs{$k};  
			#insert members of each AG back into object..
			$ret{$_}{AG_members} = $memberAGs{$k};
		}
	}
	
	#my @forbiddenSmplIDs = qw(opt totSmpls smpl_order inDir outDir baseID mocatFiltPath);
	
	#die();
	close I;
	$ret{opt}{folderStruct} = $infFoldClass;
	$ret{opt}{inDir} = join(",",@dir2dirsA) ;#if ($dir2dirs ne "");
	$ret{opt}{outDir} = $dir2out ;#if ($dir2out ne "");
	$ret{opt}{baseID} = $baseID ;#if ($baseID ne "");
	$ret{opt}{mocatFiltPath} = $mocatFiltPath ;#if ($mocatFiltPath ne "");
	$ret{opt}{GlbTmpD} =$GlbTmpD;
	$ret{opt}{NodeTmpD} =	$NodeTmpD ;
	$ret{opt}{totSmpls} = $Scnt;
	$ret{opt}{smpl_order} = \@order;
	#make this redundant at later point...
	#$ret{inDir} = join(",",@dir2dirsA) if ($dir2dirs ne "");
	#$ret{outDir} = $dir2out if ($dir2out ne "");
	#$ret{baseID} = $baseID if ($baseID ne "");
	#$ret{mocatFiltPath} = $mocatFiltPath if ($mocatFiltPath ne "");
	#@order = keys %agBP;die "@order\n";
	my $asGrpHr = emptyAssGrpsObj(\%agBP);

	return (\%ret,$asGrpHr);
}



sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub systemW{
	my ($cmddd) = $_[0];
	if ($cmddd eq ""){return;}
	my $killOnDead = 1;
	$killOnDead = $_[1] if (@_ > 1);
	$cmddd = "set -e\nulimit -c 0;\n$cmddd"; #make sure only single line excecuted.. (and not a huge core dump)
	my $stat= system $cmddd;
	#can use $? also instead for status
	if ($stat){
		print "system call \n$cmddd\nfailed with code $stat.\n";
		exit($?) if ($killOnDead);
	}
	#chomp $stat;
	return $stat;
}

sub readTabByKey{
	my ($inF) = @_;
	my %ret;
	my ($I,$OK) = gzipopen($inF,"tab file");
	my $maxTabs = 0;
	while (my $l = <$I>){
		chomp $l; my @tmp = split /\t/,$l;
		$ret {$tmp[0]} = $tmp[1];
		if (@tmp > $maxTabs){$maxTabs = @tmp;}
	}
	close $I;
	if ($maxTabs > 2){print "Warning in Tab reader: more than 2 columns were present ($maxTabs)\n";}
	return %ret;
}

sub writeFasta{
	my ($hr,$of) = ($_[0],$_[1]);
	my $maxFs = -1;
	$maxFs = $_[2] if (@_ > 2);
	my %FA = %{$hr};
	
	if ($maxFs <0){	$maxFs =  scalar(keys %FA); } #$maxFs+=1000;}
	#die "$maxFs\n";
	my $cnt=0;
	open O,">$of" or die "can't open out fasta $of\n";
	foreach my $k (keys %FA){
		$cnt++; 
		if ($k =~ m/^>/){
			print O $k."\n".$FA{$k}."\n";
		} else {
			print O ">".$k."\n".$FA{$k}."\n";
		}
		last if ( $cnt > $maxFs);
	}
	close O;
	#die $of;
}

sub convertMSA2NXS{
	my $filename = $_[0];
	my $outF = "";
	$outF = $_[1] if (@_>1);
	my $hr = readFasta($filename);
	my %FNAs = %{$hr};
	my @kk = keys %FNAs;
	my $ostr="";
	my $numtaxa = scalar(@kk);
	my $maxlength = length($FNAs{$kk[0]}); #all seqs should be same length in MSA format Format datatype=dna missing=? gap=-;
	$ostr = "#NEXUS\nBegin data;\nDimensions ntax=$numtaxa nchar=$maxlength;\nFormat datatype=dna missing=? gap=-;\nMatrix\n";

	foreach my $k (@kk) {
		my $len=length$FNAs{$k};
		die "Error nexus format conversion: $len != $maxlength in $k\n" if ($len != $maxlength);
		#if ($len<$maxlength) { my $add=$maxlength-$len; for (my $j=0; $j<$add; $j++) {$seqs[$i]=$seqs[$i].'-';}}
		$ostr.= "\n$k\t$FNAs{$k}";
	}

	$ostr .= "\n;\nend;";
	if ($outF ne ""){
		open O,">$outF" or die "Can't open nxs out file $outF";
		print O $ostr;
		close O;
	}
	return $ostr;
}


 