#!/usr/bin/perl
#effectively reduces bwt hits of incomplete genes into .clstr clustering of complete genes
#further, adds btw hits not matching criteria into the remainder fasta files
#./mergeCls.pl g/bork1/hildebra/SNP/GNMassSimu/GeneCatalog/B0/ 95
use strict;
use warnings;

sub readSam; sub readCls;sub readFasta;
sub mergeClsSam; sub secondaryCls;
#my $inD = "/g/bork1/hildebra/SNP/GNMassSimu/GeneCatalog/B0/";my $idP = "95";
my $idP = $ARGV[1];
my $inD = $ARGV[0];
mergeClsSam($inD,$idP);
secondaryCls($inD,$idP);

sub mergeClsSam(){
	my ($inD,$idP) = @_;
	my $samF = $inD."incompl.$idP.align.sam";
	my $samF2 = $inD."35compl.$idP.align.sam";
	my $clFile = $inD."compl.$idP.fna.clstr";
	my $outF = $inD."compl.incompl.$idP.fna.clstr";
	my $logf = "$inD/log/Cluster.log";
	system("mkdir -p $inD/log/");
	open LOG,">$logf";
	my ($hr1,$hr2,$totN,$totM) = readCls($clFile);
	print LOG "ComplGeneClus	$totN\nComplGeneClusMember	$totM\n";
	my %lnk  = %{$hr2};
	my %cls = %{$hr1};
	my %Hitin; my $remSeqs;
	my $totCnt = 0;
	for (my $k=0;$k<2;$k++){
		if ($k==0){
			#my $fref = readFasta("$inD/incompl.fna");
			open my $OO,">>","$inD/incompl.NAl.$idP.fna" or die "Can't open $inD/incompl.NAl.$idP.fna\n";
			($hr1,$remSeqs) = readSam($samF,$OO);
		#	print OO $remSeqs;
			close $OO;
		} elsif ($k==1){
			#my $fref = readFasta("$inD/3Pcompl.fna");
			open my $OO,">>","$inD/P35compl.NAl.$idP.fna" or die "Can't open $inD/P35compl.NAl.$idP.fna\n";
			($hr1,$remSeqs) = readSam($samF2,$OO);
		#	print OO $remSeqs;
			close $OO;
		}
		%Hitin = %{$hr1};
		foreach my $hit (keys %Hitin){
			unless (exists($lnk{$hit})){die "Can't find link $hit\n";}
			my $clCnt = scalar(split(/\n/,$cls{$lnk{$hit}}));
			$cls{$lnk{$hit}} .= "\n$clCnt\t".$Hitin{$hit};
			$totCnt++;
			#die "" if ($totCnt > 10);
		}
	}

	close LOG;

	open O,">$outF";
	foreach my $cl (keys %cls){
		my $ostr = $cl."\n".$cls{$cl};
		print O $ostr."\n";
	}
	close O;


	print "Results $outF\n";
}

sub secondaryCls(){
	my ($bdir,$cdhID) = @_;
	my $cmd="";
	#$cmd .= $cdhitBin."-est -i $bdir/P35compl.NAl.$cdhID.fna -o $bdir/P35compl.$cdhID.fna -n 9 -G 0 -M 5000 -aL 0.5 -aS 0.95 $defaultsCDH\n";
	$cmd .= "cat $bdir/P35compl.NAl.$cdhID.fna >> $bdir/incompl.NAl.$cdhID.fna\n";
	$cmd .= $cdhitBin."-est -i $bdir/incompl.NAl.$cdhID.fna -o $bdir/incompl.$cdhID.fna -n 9 -G 0 -M 5000 -aL 0.3 -aS 0.8 $defaultsCDH\n";
	$cmd .= "rm -f $bdir/incompl.NAl.$cdhID.fna\n";
	$cmd .= "rm -f $bdir/P35compl.NAl.$cdhID.fna $bdir/35compl.$cdhID.align.sam $bdir/incompl.$cdhID.align.sam\n";
}

sub readCls(){
	my ($iF) = @_;
	my %retCls; my %retRepSeq;
	open I,"<$iF";
	my $clName = "";
	my $clNum=0; my $totMem=0;
	while (my $line = <I>){
		chomp $line;
		if ($line =~ m/^>/){#open new cluster
			$clName = $line; $clNum++; $totMem++; next;
		}
		$totMem++;
		if (exists($retCls{$clName})){
			$retCls{$clName} .= "\n".$line;
		} else {
			$retCls{$clName} = $line;
		}
		if ($line =~ m/\*$/){#cluster seed
			#my @tmp = split(/\s*/,$line);
			$line =~ m/>(.*)\.\.\.\s*/;
			#print $1."\n";;
			$retRepSeq{$1} = $clName;
		}
	}

	return(\%retCls,\%retRepSeq,$clNum, $totMem);
}

sub readSam($){
	my ($iF,$OO) = @_;#$fref,
	#my %fas = %{$fref};
	open I,"<$iF" or die "Can't open $iF\n";
	my $cnt = 0; my $totLines = 0;
	my %ret; #my $fasStr = "";
	while (my $line = <I>){
		chomp $line; $totLines++;
		my @sam = split(/\t/,$line);
		my $qu = $sam[0]; my $ref = $sam[2];
		my $xtrField = join("\t",@sam[11..$#sam]);
		#die $xtrField;
		$xtrField =~ m/XM:i:(\d+)\s.*XO:i:(\d+)\s.*XG:i:(\d+)/;
		my $refL = length($sam[9]);
		#print "$1 $2 $3 $refL ".$1/$refL." ".($2+$3)/$refL."\n";
		#95% id || 90% seq length
		my $pid = $1/($refL-$2-$3);
		if ( ($2+$3)/$refL > 0.1 || $pid > 0.05){
			#TODO: mark down / read
			my $nqu = ">".$qu;
			#die "Can't find $qu in ref fasta\n" unless(exists($fas{$nqu}));
			#experimental TODO
			print $OO $nqu."\n".$sam[9]."\n";#$fas{$nqu}."\n";
			#die $nqu."\n".$sam[9]."\n";
			next;
		}
		my $mism = $1; my $gaps = $2+$3;
		$cnt++;
		$ret{$ref} = $refL-$2-$3."nt, >".$qu."... at +\/". int((1-$pid)*100) .".00%";
		#print $ret{$ref}."\n";
		#die if ($cnt == 10);
	}
	close I;
	print LOG $cnt." hits (of ".$totLines.") in $iF\n";
	return (\%ret);
}


sub readFasta($){
  my ($fil) = @_;
  my %Hseq;
  if (-z $fil){ return \%Hseq;}
  open(FAS,"<","$fil") || die("Couldn't open FASTA file $fil.");
    
     my $temp; 
     my $line; my $hea=<FAS>; chomp ($hea);
      my $trHe = ($hea);
      #my @tmp = split(" ",$trHe);
      #$trHe = substr($tmp[0],1);
      # get sequence
    while($line = <FAS>)
    {
      #next if($line =~ m/^;/);
      if ($line =~ m/^>/){
        chomp($line);
        $Hseq{$trHe} = $temp;
        $trHe = ($line);
       # @tmp = split(" ",$trHe);
		#$trHe = substr($tmp[0],1);
		#$trHe =~ s/\|//g;
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