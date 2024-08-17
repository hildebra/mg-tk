#use strict;
#use warnings;
use Mods::TamocFunc qw(sortgzblast uniq);
#use List::MoreUtils qw(uniq);
use strict; 
my $ardbfile = "/g/bork1/forslund/tara_resistome/ardb.tabs.parsed";
my $mapfile = "/g/bork1/forslund/tara_resistome/ardb_and_reforg_mapping";
my $besthitfile = "/g/bork1/forslund/tara_resistome/ardb_vs_reforg9f.overlap90shortest_famthres_or_symbol.sorted.besthit";
my $outputfile = $ARGV [1]; # "/g/bork1/forslund/tara_resistome/test.gz";
my $outputfilecats = $ARGV [2]; # "/g/bork1/forslund/tara_resistome/testCats.gz";
my $inputfile = $ARGV [0]; # "/g/scb/bork/hildebra/Tamoc/FinSoil/Sample_AV110_4/diamond/dia.ABR.blast.gz";



$inputfile = sortgzblast($inputfile);
 
# die $inputfile."\n";
 
 
$outputfilecats =~ m/^(.*\/)[^\/]+$/;
my $outD = $1;
system "mkdir -p $outD" unless (-d $outD);
 
#reads ardb tabs
my %ssym = (); #cat 1
my %scat = (); #cat 2
my %sthres = (); #cat 3

#read CAT DB
open (FH, $ardbfile);
while (<FH>) {

    my $aLine = $_;
    chomp ($aLine);

    my @words = split (/\t/, $aLine);

	if (uc ($words [3]) eq "BACA") { $words [7] = "80"; } # assume/fix typo in ardb file
    $ssym {$words [0]}{uc ($words [2])} = "";
    $scat {$words [0]}{uc ($words [3])} = "";
    $sthres {$words [0]}{uc ($words [7])} = "";
}
close (FH);

open (FH, $mapfile);
my %sym2drug = (); #specific drug resistance
while (<FH>) {
    my $aLine = $_;
    chomp ($aLine);
    my @words = split (/\t/, $aLine);
    $sym2drug {$words [1]} = $words [3];
}
close (FH);


#read in ID cutoffs
open (FH, $besthitfile); # besthit file
while (<FH>) {
	my $aLine = $_;
	chomp ($aLine);
	my @words = split (/\t/, $aLine);
	foreach my $asym (keys % {$ssym {$words [1]}}) {
		$ssym {$words [0]}{$asym} = "";
	}
	foreach my $acat (keys % {$scat {$words [1]}}) {
		$scat {$words [0]}{$acat} = "";
	}
	foreach my $athres (keys % {$sthres {$words [1]}}) {
		$sthres {$words [0]}{$athres} = "";
	}
}

close (FH);

open (FH2, " > $outputfile");
open (FH3, " > $outputfilecats");
open (FH, "zcat $inputfile |") or die "Can't opne input file $inputfile\n";
my $quOld = "";
my %wordv1; my %wordv2;my ($okhit,$retstr,$jnLine) ;
while (<FH>) {
	my $aLine = $_;
	chomp ($aLine);

	my @words = split (/\t/, $aLine);
	my $query = $words [0];
	$query =~ s/\/\d$//;
	
	if ($quOld eq "" ){$quOld = $query;
	} elsif ($quOld ne $query){  $quOld = $query;
		($okhit,$retstr,$jnLine) = workwords(\%wordv1,\%wordv2);
		
		if ($okhit){
			print FH2 "$jnLine\n";
			print FH3 "$retstr";
		}
		
		undef %wordv2; undef %wordv1;

	}
	
	#fill in array
	if (@words > 10 && $words[0] =~ m/2$/){
		$wordv2{$words [1]} = \@words;
	} elsif (@words > 10) {$wordv1{ $words [1]} = \@words;}
}

($okhit,$retstr,$jnLine) = workwords(\%wordv1,\%wordv2);

if ($okhit){
	print FH2 "$jnLine\n";
	print FH3 "$retstr";
}
close (FH);
close (FH2);
close (FH3);


system "touch $inputfile.stone";




sub combineBlasts($ $){
	my ($wh1,$wh2) = @_;
	my %bl1 = %{$wh1}; my %bl2 = %{$wh2};
	my %ret;
	
	my @allKs = uniq ( keys %bl1, keys %bl2); #
	#die "@allKs\n";
	foreach my $k (@allKs){
		my $ex1 = exists ($bl1{$k});
		unless ($ex1 && exists ($bl2{$k}) ){
			if ($ex1){$ret{$k} = $bl1{$k};
			} else { $ret{$k} = $bl2{$k};}
			next;
		}
		#pair
		$k =~ s/\/\d$/\/12/;
		$ret{$k} = $bl1{$k};
		my @hit1 = @{$bl1{$k}};
		my @hit2 = @{$bl2{$k}};
		
		my @sbss1 = sort{ $a <=> $b }(($hit1[8],$hit1[9]));
		my @sbss2 = sort{ $a <=> $b }($hit2[8],$hit2[9]);
		#sort 
		#print "$sbss1[0] > $sbss2[0]\n";
		if ($sbss1[0] > $sbss2[0]){
		#print"X";
			my @tmp = @hit1; @hit1 = @hit2; @hit2 = @tmp;
			@tmp = @sbss1; @sbss1 = @sbss2; @sbss2 = @tmp;
		}
		my @quss1 = sort { $a <=> $b }($hit1[6],$hit1[7]);my @quss2 = sort { $a <=> $b }($hit2[6],$hit2[7]);
		#overlap?
		my $overlap = 0;
		if ($sbss1[1] > $sbss2[0]){ $overlap= ( $sbss1[1] - $sbss2[0]); }#die "${$ret{$k}}[3] = $hit1[3] + $hit2[3] - ( $sbss1[1] - $sbss2[0])\n";}
		${$ret{$k}}[3] = $hit1[3] + $hit2[3] - $overlap; #ALlength
		#print "$sbss1[1] > $sbss2[0] $sbss1[0]  ${$ret{$k}}[3] $overlap\n";
		${$ret{$k}}[2] = ($hit1[2] + $hit2[2] ) /2;#%id
		${$ret{$k}}[11] = ($hit1[11] + $hit2[11]) * (1- $overlap/($hit1[3] + $hit2[3] ) );#bitscore
		${$ret{$k}}[10] = ($hit1[10], $hit2[10])[$hit1[10] > $hit2[10]];  # min($hit1[10] + $hit2[10]);
		#die "@{$ret{$k}}\n";
		#print "$k \n@sbss1 @sbss2\n@hit1\n";
	}
	
	return \%ret;
}

sub bestBlHit($){
	my ($hr) = @_;
	#my $emode=0; my $scomode = 1;
	my %blasts = %{$hr}; my $bestBit=0; my $bestk=""; my $bestScore=0;my $bestID=0; my $bestLen=0; my $bestIDever=0;
	my $k = "";
	foreach $k (keys %blasts){
	#print $k."\n";
		#my ($Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc) = ${$blRes{$k}}[11];
		#print $eval."\n";
		#sort by eval #changed from bestE -> bestScore
		#if ( ( ($emode && $bestE > $eval) || ($scomode && $bestScore< ($id * $AlLen)) ) 
		#			#&& ($eval <= $minBLE || $bitSc >= $minScore)
		#			&& ($quCovFrac == 0 || $AlLen > $DBlen{$Subject}*$quCovFrac) 
		#			&& (($fndCat || $noHardCatCheck) || exists $c2CAT{$Subject}) ) {
		#			#print "Y";
		
		#my $curScore = ${$blasts{$k}}[2] * ${$blasts{$k}}[3];
		#if ($bestScore < $curScore){$curScore = $bestScore;$bestk = $k;	}
		
		#	$bestAlLen=$AlLen;$bestE = $eval;#$bestQuery = $Query;
		
		#just sort by bitscore
		#print "@{$blasts{$k}}\n";
		#if ($bestBit < ${$blasts{$k}}[11]){
		#	$bestBit = ${$blasts{$k}}[11];  $bestk = $k;
		#}
		my $cID = ${$blasts{$k}}[2]; my $cLe = ${$blasts{$k}}[3]; my $cSc=${$blasts{$k}}[11];
		if (($bestID -5)< $cID){
			if ( ( $cID >= ($bestID *0.97) && $cID >= $bestIDever*0.9   && ( $cLe >= $bestLen * 1.15 ) )  ||  #length is just better (15%+)
					(   ($cID >= $bestID *1.03) && ( $cLe >= $bestLen * 0.9) ) ||#id is just better, while length is not too much off
					$cSc > $bestBit*0.8){   #convincing score
				$bestID = $cID;  $bestk = $k; $bestLen = $cLe; $bestBit = $cSc;
				if ($bestID > $bestIDever){$bestIDever=$bestID;}
			}
		}
	}
	if ($bestk eq ""){die "something went wrong with ABR blast:\n$k:$blasts{$k}\n";}
	return $blasts{$bestk};
}

sub workwords(){
	my ($wh1,$wh2) = @_;
	if (keys %{$wh2} != 0 && keys %{$wh1}==0){my $tmp = $wh1; $wh2 = $wh1; $wh1 = $tmp;}
	if (keys %{$wh1} == 0){return (0,"");}
	#1st combine scores
	my $whX = $wh1;
	if (keys %{$wh2} != 0){$whX = combineBlasts($wh1,$wh2);}
	#2nd: sort these, find best hit
	#print %{$whX}."\n";
	my $arBhit = bestBlHit($whX);
	#print "found\n";
	my @words = @{$arBhit};
	#contains already combined blast scores
	my $query = $words [0];
	my $subject = $words [1];
#	$qlen = $words [3];
#	$slen = $words [4];
	my $id = $words [2];
#	$bitscore = $words [5];
#	$all = $words [6];
	my $evalue = $words [10];

##	$overlap = $all / $llen;
	my @thress = sort {$a <=> $b} keys % {$sthres {$subject}};
	my $thres = $thress [0];
	my $okhit=0;my $retstr2 = ""; my $jnLine = "";
	if ($evalue < 0.00001 && $id >= $thres) {
		my %ldrugs = ();
		foreach my $sym (keys % {$ssym {$subject}}) {
			foreach my $drug (split (/\,/, $sym2drug {$sym})) {
				$ldrugs {$drug} = "";
			}
		}
		$okhit = 1;
		$jnLine = join ("\t",@words);
	#print "@words\n";
		$retstr2 = "$query\t$subject\t$id\t".join (",", keys % {$ssym {$subject}})."\t".join (",", keys % {$scat {$subject}})."\t".join (",", keys % ldrugs)."\n";
		#print FH2 "$aLine\n";
		#print FH3 "$query\t$subject\t$id\t".join (",", keys % {$ssym {$subject}})."\t".join (",", keys % {$scat {$subject}})."\t".join (",", keys % ldrugs)."\n";
	}
	return($okhit,$retstr2,$jnLine);
}
