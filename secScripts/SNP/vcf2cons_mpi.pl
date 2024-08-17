#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw( min max );


#/ei/projects/1/115b210e-aa45-4469-9eb1-d0d85d879fc0/data/results2022/IBD/HMP2I473/LOGandSUB/SNP

sub v2q_post_process;
sub sumCL;
my $indelWin=5; 

#my $line1tmp = <>;

die "vcf2cons_mpi.pl:: not enough args\n" if (@ARGV <1);
my $statfile="";
$statfile = $ARGV[0];



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
my $minDepthPar = 0; #X <= is replaced by 'N' char..

my $places = 2; my $factor = 10**$places; #for rounding frequencies

#depth related stats
my %depthStat; my %altFreq; my %allF;
while (<STDIN>) {
	$lcnt++;
	next if (/^#/);
	my @t = split;
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
			}
			&v2q_post_process($last_chr, \$seq, \$qual, \@gaps, $indelWin,$lbcnts,$lccnts,\@spos,\@sfreq,\%allF);
		}
		($last_chr, $last_pos) = ($t[0], 0);
		$seq = $qual = '';@spos=(); @sfreq = (); %allF = ();
		@gaps = ();$lbcnts=0; $lccnts =0;
		$last_chr =~ m/L=(\d+)=/;  $chromL = $1;
	}
	#print STDERR "$last_pos pos \n";
	if ($t[1] - $last_pos > 1) {
		#die @t;
		#print STDERR "add $t[1] - $last_pos\n";
		$seq .= 'n' x ($t[1] - $last_pos - 1);
		#$qual .= '!' x ($t[1] - $last_pos - 1);
    }
	die("Fatal vcf2cons_mpi.pl::@t"."\nlast_pos == $t[1] $last_pos\n") if ($last_pos == $t[1] );
    if ($t[1] - $last_pos < 0){
		die("Fatal vcf2cons_mpi.pl: [vcf2cons] unsorted input\n$t[1] - $last_pos\non line $lcnt\n@t\n") 
	}

	my ($ref, $alt) = ($t[3], $t[4]);
	if (length($ref) >= 1 && $t[7] !~ /INDEL/){# && $t[4] =~ /^([A-Za-z.])(,[A-Za-z])*$/) { # a SNP or reference
		
		#print "$ref, $alt";
		my ($b, $q);
		$q=0;
		$b = $ref;
		
		my @spl = split /:/,$t[9];
		my $splL = (scalar @spl)-1;
		my $dep= sumCL($spl[$splL]);
		
		#print "$dep\n";
		#if ($t[7] =~ /DP=(\d+)/){$dep= $1;}
		  
		if ($dep <= $minDepthPar){$b = "N";}
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
			#print STDERR "$spl[$splL]   $freq q = $altF/($altF+$Rcnt);  $ref, $alt  $last_pos  $spl[$splL]\n";
			  #print STDERR ("@t\n$q : $freq\n") if ($freq > 0.5 && $freq < 1);
			if ($dep>=1 && $freq  > 0.501 && $alt ne '.' ){
				#this is an alternate allele
				$b = $alt;$ccnts++; $lccnts++;
				push(@spos,$t[1]);push(@sfreq,$freq);
				$depthStat{$dep}{alt}++;$depthStat{$dep}{norm}--;
				$depthStat{$dep}{altF}+= $freq;
				$altFreq{int((100*($freq)))} ++;
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
		die "Fatal vcf2cons_mpi.pl:: @t\nne .\n";
		push(@gaps, [$t[1], length($t[3])]);
    }
	#print "$last_pos\n ";
	#print STDERR "@t\n" if ($alt ne ".");
    $last_pos = length($seq);#$t[1];
  }
#die;

$last_pos = $chromL+1 if ($chromL > $last_pos);
if ($last_pos > length ($seq)){
	my $ext = $last_pos - length($seq) - 1;
	$seq .= 'n' x ($ext);
}
&v2q_post_process($last_chr, \$seq, \$qual, \@gaps, $indelWin,$lbcnts,$lccnts,\@spos,\@sfreq,\%allF) if ($lcnt>0);
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
	
	if (exists($depthStat{$dep}{alt})){ #only print depth if alt exists..
		$dsStr .= $dep."\t";
		$dsStr.=$depthStat{$dep}{alt}."\t";
		$dsStr .= $depthStat{$dep}{altF}/$depthStat{$dep}{alt}."\t";
		if (exists($depthStat{$dep}{norm})){
			$dsStr.=$depthStat{$dep}{norm}."\n";
		} else{
			$dsStr .= "0\n";
		}
	} else{
		#$dsStr .= "0\t0\t";
	}

	
}

#write depth stats
open OD,">$statfile" or die "Fatal vcf2cons_mpi.pl:: Could not open outfile $statfile!!\n";

print OD $dsStr."\n";
print OD "Alternate allele freqs:\n";
foreach my $i (sort {$a <=> $b} keys %altFreq){
	if (exists($altFreq{$i})){
		print OD (($i)) ."\t$altFreq{$i}\n";
	}
}

close OD;


#finished
exit;

sub sumCL($){#sums comma sep list
	my ($in) = @_;
	my @spl = split(/,/,$in);
	my $sum=0; foreach (@spl){$sum+=$_;}
	return($sum);
}
  
sub v2q_post_process {
  my ($chr, $seq, $qual, $gaps, $l,$reports,$replaces, $ARpos,$ARfreq,$hrF) = @_;
 # print $chr." ".length($$seq)." ".@gaps."\n";
  my %allF = %{$hrF};
  
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
  my @pos = @{$ARpos};
  my @feq = @{$ARfreq};
  for my $g (@$gaps) {
	#print "@{$g} ".length($$seq)."\n";
    my $beg = $g->[0] > $l? $g->[0] - $l : 0;
    my $end = $g->[0] + $g->[1] + $l;
    $end = length($$seq) if ($end > length($$seq));
    substr($$seq, $beg, $end - $beg) = lc(substr($$seq, $beg, $end - $beg));
  }
  
  print ">$chr COV=$reports REPL=$replaces POS=".join(",",@pos)." FR=".join(",",@feq)." FREQT=$aFs\n$$seq\n"; 
  #&v2q_print_str($seq);
  #print "+\n"; &v2q_print_str($qual);
}