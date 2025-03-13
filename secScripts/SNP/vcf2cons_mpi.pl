#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw( min max );

sub v2q_post_process;
sub sumCL;
my $indelWin=5; 
#v0.1: 23.10.24: added support for merged .vcf
#v.011: 12.3.25: fix for multiple inputs
my $vcf2fastCons = 0.11; 

#my $line1tmp = <>;

die "vcf2cons_mpi.pl:: not enough args\n" if (@ARGV <3);
my $statfile="";
$statfile = $ARGV[0];
my $minDepthPar = $ARGV[1];#X <= is replaced by 'N' char..
my $minCallQual = $ARGV[2];

my $startT = time;

my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
             GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K');
my ($last_chr,  $qual, $last_pos, @gaps, @spos, @sfreq);
my $seq = "";
$last_pos=0;$last_chr="";
my $_Q=20; my $_d=1; my $_D=30000; 
my $bcnts=0; my $ccnts =0;
my $lbcnts=0; my $lccnts =0;
my $lcnt=0;
my $chromL=0;
my $disagreeCall=0;my $overridePrevCall=0; my $disagreeCallLocal=0; #stats
my $bpAdded =0; my $bpIsN =0;my $lbpIsN =0; my $entryNum=0; #stats
my $prevLine="";
my %list_of_Chrs; # to track what was already processed..
#my $minDepthPar = 0; 

print STDERR "vcf2cons_mpi v $vcf2fastCons\nParameters: minDepth=$minDepthPar minCallQual=$minCallQual\n";

my $places = 2; my $factor = 10**$places; #for rounding frequencies
my $addMode=0;#in case of merged vcf: same entry from different vcfs
my $prevDep=0; my $prevDidAlt=0; my $prevFreq=0; my $prevB=""; my $prevQual=0;
my $prevAltF=0; my $prevRcnt=0;

#depth related stats
my %depthStat; my %altFreq; my %allF;
while (my $line = <STDIN>) {
	$lcnt++;
	next if ($line =~ /^#/);
	$addMode=0;
	my @t = split/\t/,$line;
	if (@t<5){print STDERR "too short: @t\n";next;} #something wrong with file
	if ($t[7] =~ /INDEL/){next;} #print STDERR "INDEL\n";
	if ($last_chr eq ""){$last_chr = $t[0];$last_chr =~ m/L=(\d+)=/;  $chromL = $1; }#die "$chromL\n";}
	
 	if ( $last_chr ne $t[0]) {
		#print "$last_pos > ".length ($seq)."\n";
		 if ($last_chr ne ""){
			$last_pos = $chromL+1 if ($chromL > $last_pos);
			#die "$last_pos $chromL \n";
			if ($last_pos > length ($seq)){
				my $ext = $last_pos - length($seq) - 1;
				$seq .= 'n' x ($ext);
				$lbpIsN += $ext; 

			}
			&v2q_post_process($last_chr);
		}
		($last_chr, $last_pos) = ($t[0], 0);
		$seq = $qual = '';@spos=(); @sfreq = (); %allF = ();
		@gaps = ();
		
		$bcnts+=$lbcnts; $lbcnts=0; 
		$ccnts+=$lccnts;$lccnts =0; 
		$bpIsN += $lbpIsN; $lbpIsN =0;
		$prevDep=0; $prevDidAlt=0;
		$last_chr =~ m/L=(\d+)=/;  $chromL = $1;
		#ensure fasta entry is only reported once
		if (exists($list_of_Chrs{$last_chr})){
			die "$last_chr was processed already previously.\nOn line $lcnt : $line\n";
		}
		$list_of_Chrs{$last_chr}=1;
	}
	#print STDERR "$last_pos pos \n";
	if ($t[1] - $last_pos > 1) {
		#die @t;
		#print STDERR "add $t[1] - $last_pos\n";
		$seq .= 'n' x ($t[1] - $last_pos - 1);
		$lbpIsN += ($t[1] - $last_pos - 1); 
		#$qual .= '!' x ($t[1] - $last_pos - 1);
    }
	if ($t[1] == ($last_pos) ) { #
		$addMode=1;
		#die "$line\n$prevLine\n$last_pos $t[1]";
		#next;
	}
	#die("Fatal vcf2cons_mpi.pl::@t"."\nlast_pos == $t[1] $last_pos\n") if ($last_pos == $t[1] );
    if ($t[1] - $last_pos < 0){
		die("Fatal vcf2cons_mpi.pl: [vcf2cons] unsorted input\n$t[1] - $last_pos\non line $lcnt\n@t\n") 
	}

	my ($ref, $alt) = ($t[3], $t[4]);
	if (length($ref) >= 1 && $t[7] !~ /INDEL/){# && $t[4] =~ /^([A-Za-z.])(,[A-Za-z])*$/) { # a SNP or reference
		#print "$ref, $alt";
		my ($b, $q);
		$q=$t[5];  $q=0 if ($t[5] eq '.');
		my $didAlt=0;
		$b = $ref;
		
		my @spl = split /:/,$t[9];
		my $splL = (scalar @spl)-1;
		my $dep= sumCL($spl[$splL]);
		if ($addMode){
			$dep += $prevDep;
		} else {
			$prevDep = $dep;
			$prevDidAlt=0;
			$prevQual = $q;
			$prevAltF =0; $prevRcnt=0;
		}
		
		#print "$dep\n";
		#if ($t[7] =~ /DP=(\d+)/){$dep= $1;}
		  
		if ($dep < $minDepthPar){$b = 'N';$lbpIsN++}
		my @altV = ();
		if ($alt =~ m/,/ || length($alt)>1){
			@altV = split /,/,$alt;
			$alt = $altV[0];
			#die "@t\n";
		}
#		if ($t[7] =~ /;AF=([^;]+);AN=(\d)/){ #freebayes
		if ($splL == 6){#$spl[1] =~ m/,/ && $t[7] =~ /;AN=(\d)/){#$t[7] =~ /;AC=([^;]+);AN=(\d)/){
			my $NumAlleles=$1;
			my @dps =split/,/,$spl[$splL];
			my $Rcnt = shift @dps;
			my $altF = $dps[0];
			 #1/2:31:0,11,20:0:0:11,20:415,764:-92.9498,-60.2943,-56.983,-33.0199,0,-26.9993 
			 if (@dps > 2){
				my $mxI = max @dps;
				$altF = $dps[$mxI]; $alt = $altV[$mxI];
			 } 
			my $freq = $altF/($altF+$Rcnt);
			$freq = int($freq * $factor) / $factor; 
			$allF{int((100*($freq)))} ++;
			if ($addMode){$freq = ($prevAltF + $altF)/( $prevAltF + $prevRcnt+ $altF+$Rcnt);}
			#print STDERR "$spl[$splL]   $freq q = $altF/($altF+$Rcnt);  $ref, $alt  $last_pos  $spl[$splL]\n";
			  #print STDERR ("@t\n$q : $freq\n") if ($freq > 0.5 && $freq < 1);
			if ($dep>=$minDepthPar && $freq  > 0.501 && $alt ne '.' && $q>= $minCallQual ){
				#this is an alternate allele
				$b = $alt; $lccnts++;
				push(@spos,$t[1]);push(@sfreq,$freq);
				$depthStat{$dep}{alt}++;
				$depthStat{$dep}{altF} += $freq;
				$altFreq{int((100*($freq)))} ++;
				$didAlt=1;
				if (!$addMode){
					$prevDidAlt = 1 ;
					$prevFreq = $freq;
					$prevAltF = $altF; $prevRcnt = $Rcnt;
				} else {
					$prevAltF += $altF; $prevRcnt += $Rcnt;
					#$prevFreq = ($prevFreq + $freq)/2;
				}

			}
		}
		
		$depthStat{$dep}{norm}++ if (!$didAlt);
		#die "$q $b\n";
		#$b = lc($b);
		#$b = uc($b) if (($t[7] =~ /QA=(\d+)/ && $1 >= $_Q && $1 >= $_d && $1 <= $_D) );
		if ($addMode ){
			#replace previous nucl?
			if ($prevB ne $b && $b ne 'N' ){
				if ($prevB ne 'N'){ #otherwise no reason to count..
					$disagreeCall++;$disagreeCallLocal++;
				}
				if ($q> $prevQual){
					if ($prevB ne 'N'){
						$overridePrevCall++ ;
						#DEBUG
						#print STDERR "$prevLine\n$line\n $prevB  $b \n";
					}
					
					substr($seq,-1,1,$b) ;
					$prevB = $b;
					$prevQual=$q;
				} else {
					#print STDERR "$prevLine\n$line\n $prevB  $b \n";
				}
			}
			#clear some stats..
			if ($prevDidAlt){
				$depthStat{$prevDep}{alt}--;
				$depthStat{$prevDep}{altF} -= $prevFreq;
				$altFreq{int((100*($prevFreq)))} --;

				$prevDidAlt=0;$prevDep=0; 
			} else {
				$depthStat{$prevDep}{norm}--;
			}
		} else {
			$seq .= $b;
			$prevB = $b;
		}
		$prevDep = $dep;#in case more than two lines here..
		#$bpAdded++; #$bpIsN++ if ($b eq 'N' || $b eq 'n');
		
		#print STDERR "$seq\n";
		#print STDERR "$b\n";
		#die "$seq\n" if (length($seq) > 10);
		#substr($seq,$t[1],length($b)) = $b;
		#die "yes";
		#$q = int($q + 33 + .499);
		#$q = chr($q <= 126? $q : 126);
		#$qual .= $q;
		$lbcnts++; 
    } elsif ($t[4] ne '.') { # an INDEL
		die "Fatal vcf2cons_mpi.pl:: @t\nne .\n";
		push(@gaps, [$t[1], length($t[3])]);
    }
	#print "$last_pos\n ";
	#print STDERR "@t\n" if ($alt ne ".");
    $last_pos = $t[1];#length($seq);#$t[1];
	if (length($seq) != $t[1]){
		die "Sequence has been malformed:\n$seq\n $line\n$prevLine\n";
	}
	$prevLine = $line;

  }
#die;

$last_pos = $chromL+1 if ($chromL > $last_pos);
if ($last_pos > length ($seq)){
	my $ext = $last_pos - length($seq) - 1;
	$seq .= 'n' x ($ext);
}

$bcnts+=$lbcnts; $ccnts+=$lccnts;$bpIsN += $lbpIsN; 
&v2q_post_process($last_chr) if ($lcnt>0);


#print STDERR "$bcnts $ccnts\n";
#print depth stat instead
my $dsStr=""; my @cumSum;
my @deps = sort {$a <=> $b} keys %depthStat;
for (my $i=0; $i<$#deps;$i++){
	if (exists($depthStat{$i})){
		my $lval = $#cumSum;
		$lval += $depthStat{$i}{norm} if (exists($depthStat{$i}{norm}));
		$lval += $depthStat{$i}{alt} if (exists($depthStat{$i}{alt}));
		push (@cumSum,$lval);
	} else {
		push(@cumSum,$#cumSum);
	}
}
#die "@cumSum\n";
$dsStr .= "Coverage\tRef_cov\tAlt_cov\tAvg_freq\n";#header
foreach my $dep (@deps){
	
	if (exists($depthStat{$dep}{alt})){ #only print depth if alt exists..
		$dsStr .= $dep."\t";
		if (exists($depthStat{$dep}{norm})){
			$dsStr.=$depthStat{$dep}{norm}."\n";
		} else{$dsStr .= "0\n";}
		$dsStr.=$depthStat{$dep}{alt}."\t";
		if ($depthStat{$dep}{alt} > 0){
		$dsStr .=  $depthStat{$dep}{altF}/$depthStat{$dep}{alt}."\t";
		} else {$dsStr .= "0\t";}
	} else{
		#$dsStr .= "0\t0\t";
	}

	
}

#write depth stats
open OD,">$statfile" or die "Fatal vcf2cons_mpi.pl:: Could not open outfile $statfile!!\n";
print OD $dsStr."\n";
print OD "Alternate allele freqs:Alt_freq\tOccurrence\n";
foreach my $i (sort {$a <=> $b} keys %altFreq){
	if (exists($altFreq{$i})){
		print OD ((($i)) / 100)  ."\t$altFreq{$i}\n";
	}
}

close OD;

print STDERR "Finished vcf2cons_mpi.pl.\nElapsed time: " . (time - $startT ) . "s\n";
print STDERR "Total SNPs detected: $ccnts\n";
print STDERR "Total bp written: $bcnts ($bpIsN not resolved) on $entryNum entries\n";
print STDERR "Conflicting calls: $disagreeCall Resolved with second line: $overridePrevCall\n";




#finished
exit;

sub sumCL($){#sums comma sep list
	my ($in) = @_;
	my @spl = split(/,/,$in);
	my $sum=0; foreach (@spl){$sum+=$_;}
	return($sum);
}
  
sub v2q_post_process {
  my ($chr) = @_;
 # print $chr." ".length($seq)." ".@gaps."\n";
  #my %allF = %{$hrF};
  
  my $aFs = "";
  foreach my $k (0 .. 100){
	#DEBUG
	#print $k."\n";
	if (exists($allF{$k})){
		$aFs .= $allF{$k}.",";
	} else {
		$aFs .= ",";
	}
  }
  my $lengthS = length($seq);
  foreach my $g (@gaps) {
	#print "@{$g} ".length($seq)."\n";
    my $beg = $g->[0] > $lengthS ? $g->[0] - $lengthS : 0;
    my $end = $g->[0] + $g->[1] + $lengthS;
    $end = length($seq) if ($end > length($seq));
    substr($seq, $beg, $end - $beg) = lc(substr($seq, $beg, $end - $beg));
  }
  #" . join(",",@gaps) . "
  print ">$chr COV=$lbcnts REPL=$lccnts POS=".join(",",@spos)." FR=".join(",",@sfreq)." FREQT=$aFs CONFL=$disagreeCallLocal\n$seq\n"; 
  $entryNum++;
  $disagreeCallLocal=0;
  #&v2q_print_str($seq);
  #print "+\n"; &v2q_print_str($qual);
}