#!/usr/bin/env perl
use strict;
use warnings;


sub LCA_rework_tmpLines;
sub getSLVtaxo;


my ($search, $index, $sample) = ("","","");
my (@samples, @ambiglines);
my (%map, %idx, %idx2, %count, %reads, %ambig, %ambiHist);
@ambiHist{(0..200)} = 0;
#ini values for LCA
my @idThr = (97,93,93,91,88,78); #id cutoffs were a hit is still considered
my $LCAfraction = 0.9; my $lengthTolerance = 0.8;

my $id_file="/g/bork/sunagawa/shared/Falk/LCA/16S.OTU.SILVA.reference.sequences.ids";
my $map_file="/g/bork/sunagawa/shared/Falk/LCA/read2SampleMap.tsv";
my $LCAtaxf = "/g/bork3/home/hildebra/dev/lotus//DB//SLV_123_SSU.tax";
my $result_file=$ARGV[0];

#load taxonomy
my $SLVthr = getSLVtaxo($LCAtaxf); 

# create array of SILVA ids and a hash of array indices
my @array = do {
    open (ID, "<$id_file") || die "could not open $id_file: $!";
    <ID>;
};
close ID;
chomp @array;
@idx{@array} = (0..$#array);

# create 1) hash that maps read ids to sample ids, 2) array with unique sample ids, 3) hash that maps sample ids to read numbers, 4) hash of sample array indices
#open (M, "<$map_file");
#while (<M>){
#    my @F = split "\t", $_;    chomp @F;
#    $map{$F[0]}=$F[1];    push @samples, $F[1];
#} close M;

my @unique_samples = do { grep { !$reads{$_}++ } @samples };
@idx2{@unique_samples} = (0..$#unique_samples);

# create hash of samples with arrays filled with zeros as values
foreach my $k (@unique_samples){
    my @arr = ("0") x scalar @array;
    $count{$k} = [ @arr ];
    $ambig{$k} = 0;
}
my @arr = ("0") x scalar @array;

# parse usearch result
open (IN, "<$result_file");
my $previous="";
my $flag=0;
my @F;
while (<IN>){
	@F = split "\t", $_;
	#print $F[0]." ";
    if($F[0] ne $previous){
		if ($flag == 1){
			$ambig{$sample}++;
			$ambiHist{ scalar(@ambiglines) } ++;
			$flag = 0;
	# new edits
			# parse ambiglines and define LCA for this read
			#print "LCA\n";
			my $taxAR = LCA_rework_tmpLines(\@ambiglines,$SLVthr,7);
			my @tmp = @{$taxAR};
			print "$previous  @tmp\n";
			# here $count{$sample}[$index]++; where $index in the index of LCA array
			@ambiglines=();
			#die "" if ($ambig{$sample} == 2);
		}else{
			$count{$sample}[$index]++ unless $. == 1;
			$flag = 0;
		}
		
		#TODO remove after debug
		if (0){
		if (defined $map{$F[0]}){
			$sample = $map{$F[0]};
		}else{
			die "no sample id found for read: $F[0]\n";
		}
		$search = $F[1];
		if (defined $idx{$search}){
			$index = $idx{$search};
		}else{
			die "no index found for SILVA id: $search\n";
		}
		}
    }else{
		$flag=1;
	# new edits
		push (@ambiglines, [ @F ] );   
    }
    $previous=$F[0];
}
if ($flag == 1){
    $ambig{$sample}++;
}else{
    $count{$sample}[$index]++;
}
close IN;
system "touch tmp.tmp";
die "Blast file processed\n";

# create table
my @rows = ();
my @transposed = ();

# add colnames
my @col = @array;
unshift (@col, "ambiguous");
unshift (@col, "OTUrep");
push (@rows, [ @col ]);

# add ambiguous
#my @ambig = @array;
#foreach my $i (@array){
#    my $index = $idx{$i};
#    @{ $count{"ambiguous"} }[$index] = $ambig{$i};
#    $ambig[$index] = $ambig{$i};
#}
#unshift(@{ $count{"ambiguous"} }, "ambiguous");
#unshift (@ambig, "ambiguous");
#push (@rows, [ @ambig ]);

foreach my $k ( keys %count ) {
    unshift (@{ $count{$k} }, $ambig{$k});
    unshift (@{ $count{$k} }, $k);
    push(@rows, [@{ $count{$k} }]);
}

# transpose and print
for my $row (@rows) {
    for my $column (0 .. $#{$row}) {
	push(@{$transposed[$column]}, $row->[$column]);
    }
}
for my $new_row (@transposed) {
    print join "\t", @{ $new_row };
    print "\n";
}

exit;





#added by Falk
#load tax database
sub getSLVtaxo($){
	my ($ggTax) = @_;
	open TT,"<",$ggTax or die "Can't open taxonomy file $ggTax\n";
	#my @taxLvls = ("domain","phylum","class","order","family","genus");
	my %ret;
	while (my $line = <TT>){
		chomp $line;
		my @spl = split("\t",$line);
		my $tmp =  $spl[1];
		if (@spl < 2){		die("Taxfile line missing tab separation:\n".$line."\n");}
		$tmp =~ s/__;/__\?;/g;
		$tmp =~ s/__unidentified;/__\?;/g;
		$tmp =~ s/s__$/s__\?/g;
		$tmp =~ s/\s*[kpcofgs]__//g;

		$tmp =~ s/\"//g;
		my @sp2 = split(";",$tmp);
		foreach (@sp2){s/\]\s*$//; s/^\s*\[//; chomp;}
		my $taxv = join("\t",@sp2);
		#die($taxv."\n");
		$ret{$spl[0]} = \@sp2;
	}
	close TT;
	#printL "Read $refDBname taxonomy\n",0;
	return \%ret;
}

#LCA helper
sub maxTax($){
	my ($in) = @_;
	my @spl22 = @{$in};
	my $cnt=0;
	foreach (@spl22){
		last if ($_ eq "?");
		$cnt++;
	}
	return $cnt;
}
#LCA helper
sub correctTaxString($ $){
	my ($sTax2,$sMaxTaxX) = @_;
	my @ta = @{$sTax2};
	my @ta2 = @ta;
	#die "@ta\n".$ta[$sMaxTaxX]." ".$sMaxTaxX."\n";
	for (my $i=$sMaxTaxX; $i<@ta2; $i++){
		$ta2[$i] = "?";
	}
	return \@ta2;
}
#LCA helper
sub add2Tree($ $ $){
	my ($r1,$r2,$mNum) = @_;
	my %refT = %{$r1};
	my @cT = @{$r2};
	my $k="";
	#print $cT[3]."\n";
	#my $tmp = join("-",@cT); print $tmp." ".$mNum."\n";
	for (my $i =0; $i<$mNum; $i++){
		last if $cT[0] eq "?";
		if ($i == 0){$k=$cT[0];
		} else {$k .= ";".$cT[$i];}
		if (exists($refT{$i}{$k})){ $refT{$i}{$k} ++; 
		} else { $refT{$i}{$k} = 1; }
	}
	return \%refT;
}

#main LCA function
sub LCA($ $ $){
	my ($ar1,$ar2,$maxGGdep) = @_;
	my @sTax = @{$ar1};
	my @sMaxTaxNum = @{$ar2};
	if (scalar(@sTax) == 1){
		#print"early";
		my @tmpX = @{$sTax[0]};
		my @tmp = ();
		for (my $i=0;$i<$sMaxTaxNum[0];$i++){
			push (@tmp,$tmpX[$i]);
		}
		for (my $re=scalar(@tmp);$re<$maxGGdep;$re++){push(@tmp,"?");}
		return(\@tmp);
	}
	my $r1 = {};
	for (my $i=0; $i<scalar(@sTax); $i++){
		#my @temp = split($sTax[$i]);
		#print @{$sTax[$i]}[0]."  $sMaxTaxNum[$i] \n";
		#next if ($sTax[$i] =~ m/uncultured /);
		$r1 = add2Tree($r1,$sTax[$i],$sMaxTaxNum[$i]);
		
	}
	my %refT = %{$r1};
	my $fini = 0;
	my $latestHit = "";
	#determine which taxa has the highest number of hits
	my $dk;
	my $numHits = int(@sTax) + 1;
	foreach $dk (sort {$a<=>$b} (keys %refT)){
		#print $dk." ";
		my @curTaxs = keys %{$refT{$dk}};
		foreach my $tk (@curTaxs){
			#if ($dk == 2){print int($LCAfraction*$numHits). " ". $refT{$dk}{$tk}.":";}
			#if ($refT{$dk}{$tk} < $numHits){#need to get active
			if ($refT{$dk}{$tk} >= int($LCAfraction*$numHits)){
				$latestHit = $tk;
				$fini=0;
				$numHits = $refT{$dk}{$latestHit};
				last;
			#} #else {#$fini = 1;#last;}
			} else {
				$fini = 1;
				#$latestHit = $tk;
			}
		}
		
		if ($fini){last;}
	}
	#die;
	#my $winT = join("\t",@refT);
	#print "LAT ".$latestHit."\n";
	my @ret = split(";",$latestHit);
	for (my $re=scalar(@ret);$re<$maxGGdep;$re++){push(@ret,"?");}
	
	return(\@ret);
}

#get a group of usearch lines that are considered valid hits and apply LCA to these
sub LCA_rework_tmpLines(){
	my ($tmpLinesAR,$GGhr,$maxGGdep) = @_;
	my @tmpLines= @{$tmpLinesAR};
	my $debug_flag = 0;
	my ($sID,$sLength) = (0,0);
	#find best vals
	foreach my $lin2 (@tmpLines){
		$sID = ${$lin2}[2] if (${$lin2}[2] > $sID);
		$sLength = ${$lin2}[3] if (${$lin2}[3] > $sLength);
	}
	if (@tmpLines == 0 || $sID == 0){#prob no entry passed inclusion criteria
		return ([]);
	}
	my %GG = %{$GGhr};
	my  $tolerance = 3;
	my @sTax=(); my @sMaxTaxNum = ();
	foreach my $lin2 (@tmpLines){#just compare if the tax gets any better
		my @spl2 = @{$lin2};#split("\t",$lin2);
		#TODO  : recheck uparse params
		if ($spl2[2] < ($sID - $tolerance)) {next;}
		if ($spl2[3] < ($sLength * $lengthTolerance)){next;}
		my $sMax2 = 0;
		foreach (@idThr) {if($spl2[2] < $_){$sMax2++}};
		$sMax2 = 7 - $sMax2;
		unless (exists $GG{$spl2[1]} ){die "Can't find GG entry for $spl2[1]\n";}
		my $tTax = $GG{$spl2[1]} ;
		#my @tmp = @{$tTax}; print "@tmp\n";
		my $sMax3 = maxTax($tTax);
		if ($sMax3 <= $sMax2){$sMax2 = $sMax3;} 
		push(@sTax,$tTax); #push the tax string
		push(@sMaxTaxNum,$sMax2); #push the actual numbers of entries with tax in taxstring
	}
	#entry for last OTU with best results etc..
	#die "sTax not defined: LC=".@tmpLines."\n@{$tmpLines[0]}\n@{$tmpLines[1]}\n@{$tmpLines[2]}\n@{$tmpLines[3]}\n" unless ( @sTax > 0);
	my($sTaxX) = LCA(\@sTax,\@sMaxTaxNum,$maxGGdep);
	my @tmp = @{$sTaxX}; print "$sID @tmp\n";
	return ($sTaxX);#,\%retD);

}


sub splitBlastTax($ $){
	my ($blf,$num) = @_;
	my $blLines = `wc -l $blf | cut -f1 -d ' ' `;
	my $endL = int($blLines / $num); 
	if ($endL < 3000) {$endL = 3000;}
	my $subLcnt = $endL;
	my @subf;my $fcnt = 0; my $totCnt=0;
	open I,"<$blf" or die "Can't open to split $blf\n";
	my $lstHit = "";my $OO;
	while (my $l = <I>){
		$l =~ m/^(\S+)\s/;
		#my $hit = $1;
		if ($1 ne $lstHit && $subLcnt >= $endL){
			#open new file
			#print "$lstHit  $1 $subLcnt\n";
			$subLcnt = 0; $lstHit = $1;
			close $OO if (defined $OO);
			open $OO,">$blf.$fcnt"; 
			push (@subf,"$blf.$fcnt");
			$fcnt ++;
		} else {
			$lstHit = $1; #has to continue until next time a change occurs..
		}
		print $OO $l; $subLcnt++; $totCnt++;
		
	}
	close $OO; close I;
	#die $blLines." $blf\n";
	#die "@subf\n$totCnt\n";
	return @subf;
}
