#!/usr/bin/perl -w
use strict;
use warnings;


use Mods::GenoMetaAss qw( readFasta );

my $refG = shift @ARGV;
my @allFNAs = @ARGV;

my $hr = readFasta($refG); my %FR = %{$hr};
my %coPos;my %coPosRA;my %coPosOCC;

foreach my $eGf (@allFNAs){
	$eGf =~ m/\/([^\/]+)$/; my $genoN = $1; 
	#print $genoN;
	#die $genoN . "$eGf\n";
	$hr = readFasta($eGf,0); my %curG = %{$hr};
	foreach my $kk( keys %curG){
		$kk =~ m/^(\S+)/;
		my $shrtHd = $1;
		die "can't find $shrtHd fasta entry\n" unless (exists($FR{$shrtHd}));
		my $seq = $curG{$kk};
		#now find deviant positions and find NT change
		#>MM4TEC2__C1_L=46754= COV=995 REPL=5 POS=31797,31985,32191,32192,32193
		next if ($kk =~ m/ REPL=0 /);
		$kk =~ m/ POS=(.*)$/;
		#check each positions for deviations
		my @pos = split /,/,$1;
		foreach my $pp (@pos){
			my $N = uc substr($seq,$pp-1,1);
			my $O = uc substr($FR{$shrtHd},$pp-1,1);
			#print "$N $O $pp @pos\n";
			$coPos{$shrtHd}{$pp}++;
			if (exists($coPosRA{$shrtHd}{$pp})){
				#print "$coPosRA{$shrtHd}{$pp-1}  X  $O\t$N\n";
				if ($coPosRA{$shrtHd}{$pp} ne "$O\t$N") {print STDERR "ohoh\n";}
			}
			$coPosRA{$shrtHd}{$pp} = "$O\t$N";
			$coPosOCC{$shrtHd}{$pp} .= "$genoN\t";
		}
	}
}
#now rework data into something sensible
my $vcfl; my $SNPcnt=0;
#pseudo header
$vcfl = "##fileformat=VCFv4.2
##fileDate=20170116
##source=freeBayes v1.1.0-dirty
##reference=$refG
##phasing=none\n";
foreach my $chr (keys %coPos){
	#  ##contig=<ID=MM4TEC2__C1_L=46754=,length=46754>";
	$chr =~ m/=(\d+)=/;
	$vcfl .= "##contig=<ID=$chr,length=$1>\n";
}
foreach my $chr (keys %coPos){
	my @ks = sort { $coPos{$chr}{$b} <=> $coPos{$chr}{$a} } keys %{$coPos{$chr}};
	foreach my $pos (@ks){
		if ($coPos{$chr}{$pos} > 1 ){print STDERR "multi: $coPos{$chr}{$pos} $chr $pos  $coPosOCC{$chr}{$pos}\n";}
		#create the VCF line
		#MM4TEC2__C5_L=178675=	178066	.	T	G	0	.	AB=0;ABP=0;AC=0;AF=0;AN=7;AO=38;CIGAR=1X;DP=1752;DPB=1752;DPRA=57.0667;EPP=5.06748;EPPR=48.2918;GTI=0;LEN=1;MEANALT=1.33333;MQM=41;MQMR=41.9463;NS=7;NUMALT=1;ODDS=300.759;PAIRED=0.552632;PAIREDR=0.908932;PAO=0;PQA=0;PQR=0;PRO=0;QA=1316;QR=63007;RO=1713;RPL=38;RPP=85.5263;RPPR=41.8319;RPR=0;RUN=1;SAF=16;SAP=5.06748;SAR=22;SRF=665;SRP=188.96;SRR=1048;TYPE=snp;technology.ILLUMINA=1	GT:DP:AD:RO:QR:AO:QA:GL	.	0:15:15,0:15:555:0:0:0,-48.2921	.	.	.	0:4:4,0:4:156:0:0:0,-13.7396	0:14:14,0:14:524:0:0:0,-45.6663	.	.	.	.	0:7:7,0:7:254:0:0:0,-22.4531	.	0:16:14,2:14:541:2:74:0,-40.0465	.	.	0:1069:1054,14:1054:38669:14:486:0,-3316.8	0:627:605,22:605:22308:22:756:0,-1868
		$vcfl .= "$chr	$pos	.	$coPosRA{$chr}{$pos}	0.1	.	ConsensusEvidence\n";#AB=0;ABP=0;AC=0;AF=0;AN=7;AO=38;CIGAR=1X;DP=1752;DPB=1752;DPRA=57.0667;EPP=5.06748;EPPR=48.2918;GTI=0;LEN=1;MEANALT=1.33333;MQM=41;MQMR=41.9463;NS=7;NUMALT=1;ODDS=300.759;PAIRED=0.552632;PAIREDR=0.908932;PAO=0;PQA=0;PQR=0;PRO=0;QA=1316;QR=63007;RO=1713;RPL=38;RPP=85.5263;RPPR=41.8319;RPR=0;RUN=1;SAF=16;SAP=5.06748;SAR=22;SRF=665;SRP=188.96;SRR=1048;TYPE=snp;technology.ILLUMINA=1	GT:DP:AD:RO:QR:AO:QA:GL	.	0:15:15,0:15:555:0:0:0,-48.2921	.	.	.	0:4:4,0:4:156:0:0:0,-13.7396	0:14:14,0:14:524:0:0:0,-45.6663	.	.	.	.	0:7:7,0:7:254:0:0:0,-22.4531	.	0:16:14,2:14:541:2:74:0,-40.0465	.	.	0:1069:1054,14:1054:38669:14:486:0,-3316.8	0:627:605,22:605:22308:22:756:0,-1868\n";
		#die "$vcfl\n";
		$SNPcnt++;
	}
}

print STDERR "Foung $SNPcnt SNPs in consensus\n";

print $vcfl;
















