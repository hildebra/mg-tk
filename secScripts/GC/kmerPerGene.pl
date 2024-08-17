#!/usr/bin/perl
#calcs the kmer freq per gene, taking the average of gene occurrence in +- 5 genes
# perl kmerPerGene.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/
use warnings;
use strict;

use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta getAssemblPath gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
use Mods::math qw( roundF);

sub readKmers_Smpl;

my $GCd = $ARGV[0];
my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;
my ($hr1,$hr2) = readMapS($mapF,-1);
my %map = %{$hr1}; my %AsGrps = %{$hr2};
my $kmerOut = "$GCd/compl.incompl.95.fna.kmer";
system "rm -f $kmerOut" if (-e $kmerOut);

#read gene clusters jsut to get an idea of which genes are present
($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",0); $hr1 = {};
my @totGenes = sort  { $a <=> $b } keys %{$hr2};
print @totGenes."\n";
#empty mem
$hr2 = (); 
my %finalKmer; my %finalKcnt;

my $binSize = 7e6;
my @cntGR = (0,0,0,0,0,0);
for (my $subs = 0; $subs < int(scalar(@totGenes)/$binSize)+1; $subs++){
	print "\nAt genes ". $binSize*$subs . " - " . $binSize*($subs+1) . "\n\n";
	my %subHs;
	my $numGenesT = scalar @totGenes;
	for (my $i = $binSize*$subs; $i < $binSize*($subs+1); $i++){
		$subHs{$totGenes[$i]} = 1;
		last if ($i+1 >= $numGenesT);
	}
	($hr1,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx",0,\%subHs);my %cl2gene = %{$hr2}; $hr1 = {}; undef %subHs;
	print scalar(keys %cl2gene)."\n";
	my %cl2gene2; my $numGenesRep = 0;my $numKmerReps=0;
	foreach my $gene (keys %cl2gene){
		my $geneStr = $cl2gene{$gene};
		$geneStr =~ s/>//g;
		my @genegenes = split /,/,$geneStr;
		#sort contigs by length and decide on length cutoff
		my %srtH; my $gr10k=0; my $gr5k=0; my $gr2k=0; my $gr1k=0;
		foreach my $sg (@genegenes){
			$sg =~ m/_L=(\d+)=/;
			$srtH{$sg} = $1;
			$gr10k++ if ($1 > 10000);
			$gr5k++ if ($1 > 5000);
			$gr2k++ if ($1 > 2000);
			$gr1k++ if ($1 > 1000);
		}
		my $uses=0;
		if ($gr10k>2){
			$uses = $gr10k ;
			$cntGR[0]++;
		} elsif ($gr5k>2){
			$uses = $gr5k ;
			$cntGR[1]++;
		} elsif ($gr2k>2){
			$uses = $gr2k ;
			$cntGR[2]++;
		} elsif ($gr1k>2){
			$uses = $gr1k ;
			$cntGR[3]++;
		} elsif ($gr1k){
			$uses = $gr1k ;
			$cntGR[4]++;
		} else {
			$uses = scalar @genegenes;
			$cntGR[5]++;
		}
		my @srtKs = sort {$srtH{$b} <=> $srtH{$a}} keys %srtH;
		#if (@genegenes > 10){print "@srtKs"."\n\n$uses\n";}
		my $genesFound=0;
		foreach my $sg (@srtKs){
#			next if ($sg =~ m/_1$/); #very simple filter to remove genes occurring on small contigs (less info)
#			$sg =~ m/_L=(\d+)=/; next if ($1 < 2000);
			$sg =~ m/^(\S+)__(\S+)$/;
			$cl2gene2{$1}{$2} =  $gene;
			#print "$1 __ $2 ";
			$genesFound++;
			last if ($genesFound >= $uses);
		}
		if ($genesFound <= 0){#even accept genes on bad contigs..
			foreach my $sg (@genegenes){ 
				#my @spls = split /__/,$sg;
				$sg =~ m/^(\S+)__(\S+)$/;
				$cl2gene2{$1}{$2} =  $gene;
				$genesFound++;
			}
		}
		$numGenesRep++ if ($genesFound > 0);
		$numKmerReps += $genesFound;
	}
	undef %cl2gene;
	print "\n\nFound $numGenesRep unique genes, represented by $numKmerReps entries - @cntGR\n";
	print "Reading in from ".scalar( keys (%cl2gene2 )) . " samples\n"; 
	#go over each sample to read kmers...
	my $smplCnt=0;
	my $genesProcessed=0;
	foreach my $smpl (sort keys %cl2gene2){
		my $geneRd = readKmers_Smpl($smpl,$cl2gene2{$smpl},$smpl);
		print "$smplCnt:$smpl:$geneRd:".scalar(keys %{$cl2gene2{$smpl}}) ."\t";
		print "\n" if ($smplCnt % 10 == 0);
		$smplCnt++;
		undef %{$cl2gene2{$smpl}};
		$genesProcessed += $geneRd;
	}
	print "Read $smplCnt samples' and $genesProcessed kmers\n";
	#write out kmers..
	my @kKeys = sort keys %finalKmer;
	print "Writing ".scalar @kKeys . " kmers\n";
	open O,">>$kmerOut";
	foreach my $gen (@kKeys){
		print O "$gen"; my $div = 0;#$finalKcnt{$gen};
		my @curA = @{$finalKmer{$gen}};
		for (my $i=0;$i<@curA;$i++){
			$div += $curA[$i];
		}
		$div /= 100;
		for (my $i=0;$i<@curA;$i++){
#			print O "\t".round($curA[$i]/$div,2);
			print O "\t".roundF($curA[$i]/$div,1000);
		}
		print O "\n";
	}
	close O;
	%finalKmer = ();%finalKcnt = ();
}



print "CntGrps - @cntGR\n$kmerOut\nDone\n";
exit(0);





sub readKmers_Smpl{
	my $sd = $_[0]; #this is current sample
	my $hr = $_[1];
	my $smpl = $_[2];
	my %genes = %{$hr};
	my $sd2 = $sd;
	#print "$sd ";
	if (exists(  $map{altNms}{$sd}  )){		$sd2 = $map{altNms}{$sd}; 	}
	unless (exists ($map{$sd2}) ) {
		print "Can't find map entry for $sd\n"; die;
	}
	
	my $cD = $map{$sd2}{wrdir}."/";
	my $metaGD = getAssemblPath($cD);

	
	my $inK = "$metaGD/ContigStats/scaff.pergene.4kmer.pm5.gz";
	#print "$inK\n";
	my $geneCnt=0;
	my ($I,$ok) = gzipopen($inK,"K-mer per gene",1);
	die "can't open input kmer\n" if (!$ok);
	my $genesFound=0;
	while (my $l = <$I>){
		next if ($l =~ m/^Contig/);
		$l =~ m/^\S+__(\S+)\s/;
		my $locGene = $1;
		#my @spl = split /__/,$locGene1 ;
		#my $locGene =  $spl[1];
		#$locGene = $smpl."__".$locGene;
		unless (exists($genes{$locGene})){
			#print "Could not find gene $locGene\n";
			next;
		}
		#$genesFound++;
		chomp $l; my @s = split /\t/,$l; shift @s;
		my $curGCgene= $genes{$locGene};
		#print "$curGCgene  $s[0]\n";
		if (exists($finalKmer{$curGCgene})){#add to existing gene
			#print $s[0]."\t". ${$finalKmer{$curGCgene}}[0]."\n";
			for (my $i=0;$i<@s;$i++){
				${$finalKmer{$curGCgene}}[$i] += $s[$i];
			}
			#die $s[0]."\t". ${$finalKmer{$curGCgene}}[0]."\n";
			#$finalKcnt{$curGCgene}++;
		} else { #create new entry
			$finalKmer{$curGCgene} = \@s;
			#$finalKcnt{$curGCgene}=1;
		}
		$geneCnt ++;
	}
	#print "$genesFound
	close $I;
	return $geneCnt;
}


























