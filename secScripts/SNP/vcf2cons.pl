#!/usr/bin/perl -w
use strict;
use warnings;

sub v2q_post_process;
sub sumCL;
my $indelWin=5; 

#my $line1tmp = <>;
#my $statfile="";
#if ($line1tmp =~ m/#depthStat (\S+)/){	$statfile = $1;}

my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
             GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K');
my ($last_chr,  $qual, $last_pos, @gaps, @spos);
my $seq = "";
$last_pos=0;$last_chr="";
my $_Q=20; my $_d=1; my $_D=30000; 
my $bcnts=0; my $ccnts =0;
my $lbcnts=0; my $lccnts =0;
my $lcnt=0;
my $chromL=0; 
my %chrmLs;my $headIni=0; my @chrOrd; my $chrIdx = 0;
#depth related stats
my %depthStat; my %altFreq;

my $regF = $ARGV[0]; #use fai to get length and order..
open I,"<$regF" or die "fai file $regF not found!\n";while(<I>){chomp; my@spl=split/\t/;push(@chrOrd ,$spl[0]);$chrmLs{$spl[0]}=$spl[1]; }close I; $headIni=3;
#die "@chrOrd\n";

while (<STDIN>) {
	$lcnt++;
	if (/^#/){
		next if ($headIni>=2);
		if (/^##contig=<ID=([^,]+),length=(\d+)>/){
			$chrmLs{$1} = $2;  
			if ($headIni<2){
				push (@chrOrd, $1) ;
				$headIni=1;
			}
		} elsif($headIni==1) {
			@chrOrd = sort(@chrOrd);
			#die "@chrOrd\n";
			$headIni=2;
		}
		next;
	}
	#die "@chrOrd\n";
	my @t = split;
	if ($last_chr eq ""){$last_chr = $t[0];$chromL = $chrmLs{$last_chr};}#$last_chr =~ m/L=(\d+)=/;  $chromL = $1; }#die "$chromL\n";}
	die @t."\nlast_pos == t[1]\n" if ($last_pos == $t[1] );
	
 	if ( $last_chr ne $t[0]) {
		print STDERR "$last_chr ne $t[0] : $chrOrd[$chrIdx]\n";
		#print "$last_pos > ".length ($seq)."\n";
		
		 if ($last_chr ne ""){
			#insert missing chromosomes
			while ($chrIdx< @chrOrd && $last_chr ne $chrOrd[$chrIdx]){
				my $seq2 = 'n' x ($chrmLs{$chrOrd[$chrIdx]} );
				#die "$last_chr ne $chrOrd[$chrIdx]\n$chrmLs{$chrOrd[$chrIdx]}\n";
				&v2q_post_process($chrOrd[$chrIdx], \$seq2, \$qual, \@gaps, $indelWin,0,0,[]);
				$chrIdx++;
			}
			$last_pos = $chromL+1 if ($chromL > $last_pos);
			#die "$last_pos $chromL \n";
			if ($last_pos > length ($seq)){
				my $ext = $last_pos - length($seq) - 1;
				$seq .= 'n' x ($ext);
			}
			&v2q_post_process($last_chr, \$seq, \$qual, \@gaps, $indelWin,$lbcnts,$lccnts,\@spos);
			$chrIdx++;
			print STDERR "post $chrOrd[$chrIdx]\n";
		}
		($last_chr, $last_pos) = ($t[0], 0);
		$seq = $qual = '';@spos=();
		@gaps = ();$lbcnts=0; $lccnts =0;
		if (exists($chrmLs{$last_chr})){
			$chromL = $chrmLs{$last_chr};
		} else {
			$last_chr =~ m/L=(\d+)=/;  $chromL = $1;
		}
	}
	#print STDERR "$last_pos pos \n";
	if ($t[1] - $last_pos > 1) {
		#die @t;
		#print STDERR "add $t[1] - $last_pos\n";
		$seq .= 'n' x ($t[1] - $last_pos - 1);
		#$qual .= '!' x ($t[1] - $last_pos - 1);
    }
    die("[vcf2cons] unsorted input\n$t[1] - $last_pos\non line $lcnt\n") if ($t[1] - $last_pos < 0);
	if (length($t[3]) >= 1 && $t[7] !~ /INDEL/){# && $t[4] =~ /^([A-Za-z.])(,[A-Za-z])*$/) { # a SNP or reference
		
		my ($ref, $alt) = ($t[3], $t[4]);
		my ($b, $q);
		$q=0;
		$b = $ref;

		
		my @spl = split /:/,$t[9];
		my $dep= sumCL($spl[2]);
		
		#if ($t[7] =~ /DP=(\d+)/){$dep= $1;}
		  
		if ($dep<=1){$b = "N";}
		my @altV = ();
		if ($alt =~ m/,/ || length($alt)>1){
			@altV = split /,/,$alt;
			$alt = $altV[0];
			#die "@t\n";
		}
		if ($t[7] =~ /;AF=([^;]+);AN=(\d)/){ #freebayes
			$q = $1; my $NumAlleles=$2;
			
			 #1/2:31:0,11,20:0:0:11,20:415,764:-92.9498,-60.2943,-56.983,-33.0199,0,-26.9993 
			 my $freq=0;
			 if ($spl[5] =~ m/(\d+),(\d+)/){
				my $altF;
				if ($1>$2){$altF = $1; } else {$altF = $2; $alt = $altV[1];}
				$freq = $altF/($spl[3]+$altF);
			 } else {
				$freq = $spl[5]/($spl[3]+$spl[5]);
			}
			$q = $freq;
			  #print STDERR ("@t\n$q : $freq\n") if ($freq > 0.5 && $freq < 1);
			if ($q > 0.501 && $alt ne '.' ){
				if ( $dep>1 ) { 
					#this is an alternate allele
					$b = $alt;$ccnts++; $lccnts++;push(@spos,$t[1]);
					$depthStat{$dep}{alt}++;$depthStat{$dep}{norm}--;
					$depthStat{$dep}{altF}+= $freq;
					#print STDERR "$q\n"; 
					#$q=0.6;
					$altFreq{int((100*($freq)))} ++;
				}
			}
		}
		$depthStat{$dep}{norm}++;
		#die "$q $b\n";
		#$b = lc($b);
		#$b = uc($b) if (($t[7] =~ /QA=(\d+)/ && $1 >= $_Q && $1 >= $_d && $1 <= $_D) );
		$seq .= $b;
		#print STDERR "$seq\n";
		#print STDERR "$b\n";
		#die "$seq\n" if (length($seq) > 10);
		#substr($seq,$t[1],length($b)) = $b;
		#die "yes";
		#$q = int($q + 33 + .499);
		#$q = chr($q <= 126? $q : 126);
		#$qual .= $q;
		$bcnts++;$lbcnts++; 
    } elsif ($t[4] ne '.') { # an INDEL
		die "@t\nne .\n";
		push(@gaps, [$t[1], length($t[3])]);
    }
	#print "$last_pos\n ";
    $last_pos = length($seq);#$t[1];
}
#die;

while ($chrIdx < @chrOrd && $last_chr ne $chrOrd[$chrIdx] ){
	#print STDERR "$last_chr != $chrOrd[$chrIdx]  $chrIdx\n";
	my $seq2 = 'n' x ($chrmLs{$chrOrd[$chrIdx]} );
	
	&v2q_post_process($chrOrd[$chrIdx], \$seq2, \$qual, \@gaps, $indelWin,0,0,[]);
	$chrIdx++;
}
$last_pos = $chromL+1 if ($chromL > $last_pos);
if ($last_pos > length ($seq)){
	my $ext = $last_pos - length($seq) - 1;
	$seq .= 'n' x ($ext);
}
&v2q_post_process($last_chr, \$seq, \$qual, \@gaps, $indelWin,$lbcnts,$lccnts,\@spos);
$chrIdx++;

#check if unfinished chroms are present..
for (my $j=$chrIdx;$j<@chrOrd;$j++){
	#print STDERR "post   $chrOrd[$j]\n";next;
	my $chrN = $chrOrd[$j];
	my $seq2 = 'n' x ($chrmLs{$chrN} );
	&v2q_post_process($chrN, \$seq2, \$qual, \@gaps, $indelWin,0,0,[]);
}


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
foreach my $dep (@deps){
	$dsStr .= $dep."\t";
	if (exists($depthStat{$dep}{alt})){
		$dsStr.=$depthStat{$dep}{alt}."\t";
		$dsStr .= $depthStat{$dep}{altF}/$depthStat{$dep}{alt}."\t";
	} else{
		$dsStr .= "0\t0\t";
	}
	if (exists($depthStat{$dep}{norm})){
		$dsStr.=$depthStat{$dep}{norm}."\n";
	} else{
		$dsStr .= "0\n";
	}
	
}
print STDERR $dsStr."\n";
print STDERR "Alternate allele freqs:\n";
foreach my $i (sort {$a <=> $b} keys %altFreq){
	if (exists($altFreq{$i})){
		print STDERR (($i)) ."\t$altFreq{$i}\n";
	}
}



exit;

sub sumCL($){#sums comma sep list
	my ($in) = @_;
	my @spl = split(/,/,$in);
	my $sum=0; foreach (@spl){$sum+=$_;}
	return($sum);
}
  
sub v2q_post_process {
  my ($chr, $seq, $qual, $gaps, $l,$reports,$replaces, $ARpos) = @_;
 # print $chr." ".length($$seq)." ".@gaps."\n";
  my @pos = @{$ARpos};
  for my $g (@$gaps) {
	#print "@{$g} ".length($$seq)."\n";
    my $beg = $g->[0] > $l? $g->[0] - $l : 0;
    my $end = $g->[0] + $g->[1] + $l;
    $end = length($$seq) if ($end > length($$seq));
    substr($$seq, $beg, $end - $beg) = lc(substr($$seq, $beg, $end - $beg));
  }
  
  print ">$chr COV=$reports REPL=$replaces POS=".join(",",@pos)."\n$$seq\n"; #&v2q_print_str($seq);
  #print "+\n"; &v2q_print_str($qual);
}