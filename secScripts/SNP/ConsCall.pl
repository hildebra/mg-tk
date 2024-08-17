#!/usr/bin/env perl
#builds from a contig file several stats to seperate contigs into single species

#./SNPcalls.pl /g/bork3/home/hildebra/data/AnnaPry/161013_M00758_0551_000000000-ATWLT/Bsubtillis168.fasta /g/scb/bork/hildebra/Tamoc/ANNA/GlbMap/BS168 BS168
use warnings;
use strict;
use Mods::GenoMetaAss qw(readMap median);
use Mods::Subm qw(qsubSystem emptyQsubOpt);

sub merge_vcf;
sub createConsensus;

my $QSBoptHR = emptyQsubOpt(1,"");

my $vcfcnsScr = "perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/SNP/vcf2cons.pl ";

my $frDir = "/g/bork3/home/hildebra/bin/freebayes/bin/";#required for some python scripts in that dir
my $frbBin = "$frDir/freebayes";
my $vtBin = "/g/bork3/home/hildebra/bin/vt/vt";
my $bcBin = "/g/bork3/home/hildebra/bin/bcftools-1.3.1/./bcftools"; 


my $par =0; #trigger to start parallizing freebayes
my $totRds = 10; #here you need to figure out, how many reads are in your bam
if ($totRds < 1e6){$myPar=0;
} else {$myPar=1;}

#prob some for loop around your bams..
my @vcfCons; my @vcfConsDep;
for (my $i=0;$i<10;$i++){
	my $tar = "something.bam";
	my $oVcfCons = "/sad/out.vcf";
	#arguments: inbam, 
	push(@vcfCons,$oVcfCons);
	push(@vcfConsDep,createConsensus($tar,$oVcfCons,$i,0));
}

#merge all the vcfs per sample, that were created before
merge_vcf($conFastaDir,\@vcfConsDep,$refFA,$outVCF);


sub merge_vcf{
	my ($vcfDir,$arDeps,$ref,$finalVCF) = @_;
	#die "@$arVCFs\n";
	my $mrgCmd = "$bcBin merge $vcfDir/*.pre.vcf.gz > $finalVCF.1 \n";
	$mrgCmd .= "$vtBin normalize -r $ref $finalVCF.1 >  $finalVCF";
	$mrgCmd .= "\nbgzip $finalVCF; tabix -p vcf $finalVCF.gz;\n";
	$mrgCmd .= "rm $finalVCF.1";
	$mrgCmd .= "rm $vcfDir/*.pre.vcf.gz";

	die "$mrgCmd\n";
	my ($dep,$qcmd) = qsubSystem($qsubDir."mergeVCF.sh",$mrgCmd,1,"40G","ConsMrg",join(";",@$arDeps),"",1,[],$QSBoptHR);
	
}
sub getRegionsBam(){
	my @curReg;
	my $regionFile = "$odir/regions_par.txt";
	system "samtools faidx $refFA" unless (-e "$refFA.fai");
	if (system "python $frDir/fasta_generate_regions.py $refFA.fai $splitFAsize > $regionFile\n"){
		print "python $frDir/fasta_generate_regions.py $refFA.fai $splitFAsize > $regionFile\n";
		die "Can't gemerate region file:\n";
	}
	open my $handle, '<', $regionFile;
	chomp(@curReg = <$handle>);
	close $handle;
	return (\@curReg);
}

sub createConsensus($ $ $ $ ){
	my ($tar,$oVcfCons,$x,$overwrite)  = @_;
	my $cmd = ""; my $rdep="";
	my $hereCtgs = 0;
	if (-e $oVcfCons){
		#only works for fastas
		#my $tmps = `grep -c '>' $oVcfCons`; $tmps=~m/^(\d+)/;$hereCtgs=$1;
	}
	#print "$tarCtgs != $hereCtgs\n";
	#die;
	#if ($tarCtgs != $hereCtgs){$overwrite=1;}
	if (!-e $oVcfCons && !-e "$oVcfCons.gz"){$overwrite=1;}

	#	my $cmd = "$smtBin bam2fq $tar | $seqtkBin seq -A >$ofasCons";#only gives reads back
#	$cmd = "$bam2cnsBin --bam $tar  --coverage 10 --qual-weighted --qv-offset 28 --prefix $ofasCons\n";#--ignore-weak-reads=20 
#	$cmd = "$smtBin mpileup -B  -C 20 -d 5000 -Q 28 -v -u -f $refFA  $tar | $vcfcnsScr  >$ofasCons\n";
	#bcftools consensus -f /g/bork3/home/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa test
	$cmd = "$frbBin -f $refFA  -u -i -m 30 -q 30 -C 1 -F 0.4 -k --pooled-continuous --report-monomorphic  --min-repeat-entropy 1 ";
	$cmd .= "--use-best-n-alleles 2 --haplotype-length 0 -G 1 "; #|  $vcfTools/vcfbreakmulti
	my @allDeps2;
	#die $cmd;
#	system $cmd;
	#implement in parallel as too slow in single core mode :/
	my $tmpOut = "$odir/${name[$idx]}.$x.cons.vcf";
	my @curReg = ("1");
	if ($myPar){
		my $refAR = getRegionsBam();
		@curReg = @{$refAR};
	}
	my $qsubDir = $odir."/qsubs/";
	system "mkdir -p $qsubDir" unless (-d $qsubDir);
	for (my $i=0;$i<@curReg;$i++){
		my $cmd2 = $cmd;
		if ($myPar){
			$cmd2 .= " --region '$curReg[$i]' $tar > $tmpOut.$i \n";
		} else {
			$cmd2 .= " $tar > $oVcfCons\n";
			#two options: 1) create vcf
			$cmd2 .= "bgzip $oVcfCons; tabix -p vcf $oVcfCons.gz\n";
			#option 2) create fna via perl script that I sent before ($vcfcnsScr)
			#| $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat\n";
		}
		#if (-s $tmpOut){[print " tmp exists ";next;}
		#die $cmd2;
	
		if ($overwrite ){
			die "overwr\n";
			my ($dep,$qcmd) = qsubSystem($qsubDir."FB_Cons$x.$i.sh",$cmd2,1,"5G","FBC$x.$i","","",1,[],$QSBoptHR);
			push (@allDeps2,$dep);
			$rdep =$dep;
		}
	}
	
	my $postcmd ="";

	if ($myPar && $overwrite ){
		my $vcf1stHd = "/g/bork3/x86_64/bin/python /g/bork3/home/hildebra/bin/vcflib/bin/vcffirstheader";
		my $vcfStrSrt = "/g/bork3/home/hildebra/bin/vcflib/bin/./vcfstreamsort";
		$postcmd .= "cat $tmpOut.* | $vcf1stHd | $vcfStrSrt -a > $oVcfCons\n";
		$postcmd .= "bgzip $oVcfCons; tabix -p vcf $oVcfCons.gz\n";
		#| $vcfcnsScr  >$ofasCons 2> $ofasCons.depStat\n\n";
		$postcmd .= "rm $tmpOut.*\n";
		#die "$postcmd\n";
		my ($dep,$qcmd) = qsubSystem($qsubDir."postCmd.sh",$postcmd,1,"40G","Cons$x",join(";",@allDeps2),"",1,[],$QSBoptHR);
		$rdep =$dep;
		#die;
	}
	if(-e $oVcfCons && !-e "$oVcfCons.gz"){
		my $fix .= "bgzip $oVcfCons; tabix -p vcf $oVcfCons.gz\n";
		system "$fix\n";
	}
	return $rdep;
}