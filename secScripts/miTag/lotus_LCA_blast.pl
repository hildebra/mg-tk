#!/usr/bin/env perl
#assigns LCA tax to LSU / SSU / ITS seq fragments (extracted from metag)
#ex ./lotus_LCA_blast.pl /g/scb/bork/hildebra/Tamoc/FinSoil/Sample_G2677_115/ribos FS31 20
#  ./lotus_LCA_blast.pl /g/scb/bork/hildebra/SNP/SimuG/simulated_FungiSA_5/ribos S5 20 /tmp/hildebra/simu/S5/
use strict;
use warnings;
use threads;
use Mods::GenoMetaAss qw(reverse_complement_IUPAC readFasta);
use Mods::IO_Tamoc_progs qw(getProgPaths );



sub getTaxForOTUfromRefBlast;
sub splitBlastTax;
sub writeBlastHiera; sub runBlastLCA;
sub merge;
sub fastq2fna; sub flashed;
sub rework_tmpLines;

#binaries
my $flashBin = getProgPaths("flash");#"/g/bork3/home/hildebra/bin/FLASH-1.2.10/flash";
my $lambdaIdxBin = getProgPaths("lambdaIdx");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda_indexer";
my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";


#some more pars for LCA	
my @idThr = (97,95,93,91,88,85);
my $lengthTolerance = 0.85;
my $maxHitOnly = 0;
my $DoPar = 0;

if (@ARGV ==0 ){die"no Input arguments given\n";}

my $inD = $ARGV[0]; #expects TAMOC output files in this dir
$inD .= "/" unless ($inD =~ m/\/$/);
#dir to store LCA
my $outdir = $inD."ltsLCA/";
#dir to store merged reads in
my $mergeD = $inD."flashMerge/";
my $SmplID = $ARGV[1];#Name of the Sample
my $BlastCores = $ARGV[2]; #parrallel execution
my $tmpD = ""; #tmp dir for faster I/O
if (@ARGV > 4){
	$tmpD = $ARGV[4] ;
	system "rm -rf $tmpD;mkdir -p $tmpD";
}
my $DBdir = $ARGV[3];

#datbases
my $LSUdbFA = "$DBdir/SLV_128_LSU.fasta";
my $LSUtax = "$DBdir/SLV_128_LSU.tax";
my $SSUdbFA= "$DBdir/SLV_128_SSU.fasta";
my $SSUtax = "$DBdir/SLV_128_SSU.tax";
#my $ITSdbFA = "$DBdir/sh_refs_qiime_ver7_99_02.03.2015.fasta";my $ITStax = "$DBdir/sh_taxonomy_qiime_ver7_99_02.03.2015.txt";
my $ITSdbFA = "$DBdir/ITS_comb.fa";my $ITStax = "$DBdir/ITS_comb.tax";
my $PR2dbFA = "$DBdir/gb203_pr2_all_10_28_99p.fasta";
my $PR2tax = "$DBdir/PR2_taxonomy.txt";

#some cleanups (includes prev runs)
system "rm -f $inD/reads_SSU.fq $inD/reads_LSU.fq $inD/reads_ITS.fq $inD/*.blast $outdir/*riboRun_bl";

system "mkdir -p $outdir" unless (-e $outdir);

my $blMode = 2; #2=lambda,1=Blast,3=sortmeRNA
my $inputOK=1; #flag
my $curStone = ""; #flag file for each subpart (SSU,LSU,ITS)
my @allInFa = ($inD."reads_LSU.r1.fq",$inD."reads_LSU.r2.fq",$inD."reads_ITS.r1.fq",$inD."reads_ITS.r2.fq",$inD."reads_SSU.r1.fq",$inD."reads_SSU.r2.fq");

#pretty circular safety catch.. maybe remove later?
if (!-z $allInFa[2] && (-z $allInFa[0] ||  -z $allInFa[1])){#LSU obviously wrong
	system "rm -f $inD/LSU_pull.sto $outdir/LSU_ass.sto";
}
if (!-z $allInFa[2] && (-z $allInFa[4] || -z $allInFa[5])){#SSU obviously wrong
	system "rm -f $inD/SSU_pull.sto $outdir/SSU_ass.sto";
}

#all are empty (also wrong)
if (-z $inD."reads_LSU.r1.fq" && -z $inD."reads_LSU.r2.fq" && -z $inD."reads_ITS.r1.fq" && -z $inD."reads_ITS.r2.fq"
			&& -z $inD."reads_SSU.r1.fq" && -z $inD."reads_SSU.r2.fq" &&
			-z $inD."reads_LSU.r1.fq.gz" && -z $inD."reads_LSU.r2.fq.gz" && -z $inD."reads_ITS.r1.fq.gz" && -z $inD."reads_ITS.r2.fq.gz"
			&& -z $inD."reads_SSU.r1.fq.gz" && -z $inD."reads_SSU.r2.fq.gz"){
	print "input seems completely wrong.. run has to be repeated\n";
	system "rm -f $mergeD/*_pull.sto";
	exit(11);
}

#all done already
if (-e "$outdir/Assigned.sto" && -e "$outdir/LSU_ass.sto" && -e "$outdir/ITS_ass.sto" && -e "$outdir/SSU_ass.sto"){
	print "All assigned already\n";
	exit (0);
}



#------------- LSU ------------------
$curStone = "$outdir/LSU_ass.sto";
unless (-e $curStone && -e "$outdir/LSUriboRun_bl.hiera.txt"){
	my @dbfa = ($LSUdbFA); my @dbtax = ($LSUtax);
	#merge
	my $go = flashed($inD."reads_LSU.r1.fq",$inD."reads_LSU.r2.fq",$mergeD,"LSU",$inD); 
	if ($go==1){
		#sim search & LCA
		runBlastLCA($mergeD."LSU",\@dbfa,\@dbtax,"LSUriboRun_bl",$blMode,$SmplID,"LSU",$go);
		system "touch $curStone";
	} elsif($go==2) {$inputOK=0;print"Problem with LSU primary input; run needs to be repeated";}
}

#die();


#------------- SSU ------------------
$curStone = "$outdir/SSU_ass.sto";
unless (-e $curStone && -e "$outdir/SSUriboRun_bl.hiera.txt"){
	my @dbfa = ($PR2dbFA,$SSUdbFA); my @dbtax = ($PR2tax,$SSUtax);
	my $go = flashed($inD."reads_SSU.r1.fq",$inD."reads_SSU.r2.fq",$mergeD,"SSU",$inD);
	if ($go==1){
		runBlastLCA($mergeD."SSU",\@dbfa,\@dbtax,"SSUriboRun_bl",$blMode,$SmplID,"SSU",$go);
		system "touch $curStone";
	} elsif($go==2) {$inputOK=0; print"Problem with SSU primary input; run needs to be repeated";}
}


#------------- ITS ------------------ 
$curStone = "$outdir/ITS_ass.sto"; 
unless (-e $curStone && -e "$outdir/ITSriboRun_bl.hiera.txt"){
# $go checks for empty files (no assigned ITS in metag)
	my @dbfa = ($ITSdbFA); my @dbtax = ($ITStax);
	my $go = flashed($inD."reads_ITS.r1.fq",$inD."reads_ITS.r2.fq",$mergeD,"ITS",$inD);
	if ($go==1){
		runBlastLCA($mergeD."ITS",\@dbfa,\@dbtax,"ITSriboRun_bl",$blMode,$SmplID,"ITS",$go);
		system "touch $curStone";
	} elsif($go==2) {$inputOK=0;print"Problem with ITS primary input; run needs to be repeated";}
}
#die();
#cleanup
system "rm -r $tmpD" if ($tmpD ne "");
if ($inputOK){
	system "touch $outdir/Assigned.sto";
	system "gzip ".join (" ",@allInFa);
	exit(1);
}
print "Finished\n";
exit(0);


#flash merge & some file checks
sub flashed($ $ $ $ $){
	my ($r1,$r2,$outD,$outT,$primD) = @_;
	if ( -e $r2.".gz"){system "rm -f $r2; gunzip $r2.gz";}
	if ( -e $r1.".gz"){system "rm -f $r1; gunzip $r1.gz";}
	if (    !-e $r2 || !-e $r1 || (-z $r1 && !-z $r2) || (!-z $r1 && -z $r2) 
	||  (`wc -l $r1 | cut -f1 -d' '` != `wc -l $r2 | cut -f1 -d' '` ) ){
		print "Empty input $r1\n";
		system "rm -f $primD/$outT"."_pull.sto";
		#exit 33;
		return 2;
	}
	if (-z $r1 && -z $r2){return 0;}
	print "running flash..\n";
	my $mergCmd = "$flashBin -M 200 -o $outT -d $outD -t $BlastCores $r1 $r2";
	if (system $mergCmd){
		print "\n$mergCmd\nfailed\n";
		system "rm -f $primD/$outT"."_pull.sto";
		return 2;
	}# or 
	return 1;
}

sub fastq2fna($){
	my ($in) = @_;
	print $in."\n";
	my $out = $in;
	$out =~ s/fastq$/fna/;
	system "cat $in | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\\t' '\\n' > $out";
	#die $out;
	return $out;
}

#find reads that could not match any LCA and redo with different DB
sub findUnassigned($ $ $ ){
	my ($BTr,$Fr,$outF) = @_;
	my %BT = %{$BTr}; my %Fas = %{$Fr};
	#my @t =keys %Fas; print "$t[0]\n";
	my $cnt =0; my $dcn =0 ;
	my @kk = keys %Fas;
	my $totFas = @kk;
	
	if ($outF eq ""){
		foreach my $k (keys %BT){
			my @curT = @{$BT{$k}};
			if ( @curT==0 || $curT[0] eq "?" ){$dcn++;} #|| $curT[1] eq "?"
			$cnt++;
		}
		if (@kk == 0){
			print "Total of ". ($cnt-$dcn)." / $cnt reads have LCA assignments\n";
		} else {
			print "$dcn / $cnt reads failed LCA assignments, checked $totFas reads.\n";
		}
		return;
	}
	foreach my $k (keys %BT){
		my @curT = @{$BT{$k}};
		#print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n" 
		if ( @curT==0 || $curT[0] eq "?" ){#|| $curT[1] eq "?"){
			delete $BT{$k};
			$dcn++;
			
			#print ">".$k."\n".$Fas{$k}."\n";
		} #else {print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n" ;}
		else {
			die "Can't find fasta entry for $k\n" unless (exists $Fas{$k});
			delete $Fas{$k};
		}
		$cnt ++;
		#die if ($cnt ==100);
	}
	@kk = keys %Fas;
	print "$dcn / $cnt reads failed LCA assignments\nWriting ".@kk." of previous $totFas reads for next iteration.\n";
	open O,">$outF" or die "can;t open unassigned fasta file $outF\n";
	foreach my $k(@kk){
		print O ">".$k."\n".$Fas{$k}."\n";
	}
	close O;
	return ($outF,\%BT,$dcn,\%Fas);
}

#main routine that does sim search & starts the LCA
sub runBlastLCA(){
	my ($queryO,$DBar,$DBtaxar,$id,$doblast,$SmplID,$MKname,$go) =@_;
	my $taxblastf = $outdir."$id";
	if ($tmpD ne ""){
		$taxblastf = $tmpD."$id";
	}
	my $doInter=1; my $doQuery=1; #fine control of what subparts to do..
	my $r1 = $queryO; my $r2 = $queryO;
	my $interLeaveO = "";
	my $hof1 = "$outdir/$id.hiera.txt"; my $hof2 = "$outdir/$id"."_inter.hiera.txt";
	#check for empty input
	if ($go == 2){ return ;}
	if (!$go){
		print "$id has empty input files\n";
		system "touch $hof1";
		return;
	}

	unless($queryO =~ m/\.fna$/){#merged fastq that need to be processed
		$queryO .= ".extendedFrags.fastq";
		$r1 .= ".notCombined_1.fastq";
		$r2 .= ".notCombined_2.fastq";
		$interLeaveO = $outdir."inter$id.fna";
		#print $r1." V\n";
		$r1 = fastq2fna($r1); $r2 = fastq2fna($r2); $queryO = fastq2fna($queryO);

		merge($r1,$r2,$interLeaveO) if ($doInter);
	} else {
		print "no interLeaveO\n";
		$doInter=0;
	}
	if (-z $interLeaveO){print "No interleaved files $id\n";$doInter=0;}
	if (-z $queryO){print "No merged files $id\n";$doQuery=0;}
	my $taxblastf2=$taxblastf.".i.blast";
	#my $BlastCores = 20;

	my $fasrA = readFasta( $queryO,1); 
	#my %fas=%{$fasrA}; my @t = keys %fas; die "$queryO\n@t\n$t[0]\n";
	my $fasrI = readFasta($interLeaveO,1); 
	#print $DBar."\n";
	my @DBa = @{$DBar}; my @DBtaxa = @{$DBtaxar};
	my $query= $queryO; my $interLeave = $interLeaveO;
	my $BlastTaxRi = {}; my $BlastTaxR = {};
	my $fullBlastTaxRi = {}; my $fullBlastTaxR = {};
	my $simName = ""; my $leftover = 111;
	my $contInter = 1; my $contQuery = 1;

	for (my $DBi=0;$DBi<@DBa; $DBi++){
		my $DB = $DBa[$DBi]; my $DBtax = $DBtaxa[$DBi];
		print "Running sim search $doblast..\n";
		if (-z $interLeave){$doInter=0;} else {$doInter=1;}
		if (-z $query){$doQuery=0;} else {$doQuery=1;}

		if ($doblast==1){
			$taxblastf.=".blast"; $simName="blast";
			print "Running Blast\n";
			my $mkBldbBin = "/g/bork5/hildebra/dev/lotus//bin//ncbi-blast-2.2.29+/bin/makeblastdb";
			my $blastBin = "/g/bork5/hildebra/dev/lotus//bin//ncbi-blast-2.2.29+/bin/blastn";
			my $cmd = "$mkBldbBin -in $DB -dbtype 'nucl'\n";
			unless (-f $DB.".nhr"){	system($cmd);}
			my $strand = "both";
			#-perc_identity 75
			$cmd = "";
			if ($doQuery && $contQuery){
				$cmd .= "$blastBin -query $query -db $DB -out $taxblastf -outfmt 6 -max_target_seqs 50 -evalue 0.1 -num_threads $BlastCores -strand $strand \n"; #-strand plus both minus
			}
			unless ( -e $taxblastf){	print "Running blast on combined files\n";
				system($cmd);
			} else {	print "Blast output $taxblastf does exist\n";}
			if ($doInter){
				$cmd .= "$blastBin -query $interLeave -db $DB -out $taxblastf2 -outfmt 6 -max_target_seqs 50 -evalue 0.1 -num_threads $BlastCores -strand $strand \n"; #-strand plus both minus
				unless ( -e $taxblastf2){	print "Running blast on interleaved files\n";system($cmd);
				} else {	print "Blast output $taxblastf2 does exist\n";}
			}
		} elsif ($doblast==2){
			$simName = "lambda";
			if (!-f $DB.".dna5.fm.sa.val"  ) {
				print "Building LAMBDA index anew (may take up to an hour)..\n";
				my $cmdIdx = "$lambdaIdxBin -p blastn -t $BlastCores -d $DB";
				#system "touch $DB.dna5.fm.lf.drv.wtc.24";
				if (system($cmdIdx)){die ("Lamdba ref DB build failed\n$cmdIdx\n");}
			}
			print "Starting LAMBDA similarity search..\n";
			my $tmptaxblastf = "$outdir/tax.m8";
			$tmptaxblastf = "$tmpD/tax.m8" unless ($tmpD eq "");
			my $cmd = "";
			$cmd = "";
			if ($doQuery && $contQuery){
				$cmd .= "$lambdaBin -t $BlastCores -id 75 -nm 200 -p blastn -e 1e-5 -so 7 -sl 16 -sd 1 -b 5 -pd on -q $query -d $DB -o $tmptaxblastf\n";
			}
			$cmd .= "\nmv $tmptaxblastf $taxblastf\n";
			if ($doInter && $contInter){
				$cmd .= "$lambdaBin -t $BlastCores -id 75 -nm 200 -p blastn -e 1e-5 -so 7 -sl 16 -sd 1 -b 5 -pd on -q $interLeave -d $DB -o $tmptaxblastf\n";
				$cmd .= "\nmv $tmptaxblastf $taxblastf2\n";
				#unless ( -e $taxblastf2){	
				print "Running blast on interleaved files\n";
			}
			#die $cmd;
			system($cmd);
			#} else {	print "Blast output $taxblastf2 does exist\n";}}
			#system $cmd ;#or die "\n$cmd\n failed\n";
		} elsif ($doblast==3) {
			$simName = "smRNA";
			my $smrnaBin = "/g/bork5/hildebra/bin/sortmerna-2.0/sortmerna";
			my $smrnaMkDB = "/g/bork5/hildebra/bin/sortmerna-2.0/indexdb_rna";
			my $refDB = "$DB,$DB.sidx";
			unless (-e "$DB.sidx.stats"){	print "making DB..";system "$smrnaMkDB --ref $refDB";	print " done\n";}
			#die $outdir."inter$id.fna\n";
			print "Running sortmerna merge\n";
			my $cmd = "$smrnaBin --best 50 --reads $interLeave ";
			$cmd .= "--blast 1 -a 20 -e 0.1 -m 10000 --paired_in --fastx --aligned $taxblastf.i ";
			$cmd .= "--ref $refDB\n";
			print "$taxblastf.i\n";
			system $cmd;
			print "Running sortmerna\n";
			$cmd = "";
			if ($doQuery && $contQuery){
				$cmd .= "$smrnaBin --best 50 --reads $query ";
				$cmd .= "--blast 1 -a 20 -e 0.1 -m 10000 --aligned $taxblastf ";
				$cmd .= "--ref $refDB\n";
			}
			system $cmd;
			$taxblastf.=".blast";
			print "done sortmeRNA assignment\n";
		}

		#return; #for testing purpose

		#rewriteBlastSmplNm($taxblastf,$SmplID,$MKname); #not required, as moved to catchLSUSSU.pl script.. better to do on fastq's
		#system "sort $taxblastf > $taxblastf.tmp;rm $taxblastf;mv $taxblastf.tmp $taxblastf";
		print "Running LCA..\n";
		print $DBtax."\n";
		my %GG = getGGtaxo($DBtax);
		print "Read TaxDB\n";
		my $BlastTaxR = {}; 
		if ($doQuery && $contQuery){
			if ($DoPar){
				my @subf = splitBlastTax($taxblastf,$BlastCores);
				#paralellize getTaxForOTUfromRefBlast
				my @thrs;
				for (my $i=0;$i<@subf;$i++){
					$thrs[$i] = threads->create(\&getTaxForOTUfromRefBlast,$subf[$i],\%GG,0,1);
				}
				#combine tax object
				for (my $i=0;$i<@subf;$i++){my $ret = $thrs[$i]->join();$BlastTaxR = {%$BlastTaxR,%$ret};}
			} else {
				$BlastTaxR = getTaxForOTUfromRefBlast($taxblastf,\%GG,0,$BlastCores);
			}
			#check for unassigned reads
			
			$query = $tmpD."Query__U$DBi".".fna" if ($tmpD ne "");
			
			if ($DBi < (@DBa-1)){
				($query,$BlastTaxR,$leftover,$fasrA)= findUnassigned( $BlastTaxR,$fasrA,$query );
			} else {findUnassigned($BlastTaxR,$fasrA,"");}
			$contQuery=0 if ($leftover ==0 );
		}
		#save good hits..
		$fullBlastTaxR = {%$fullBlastTaxR,%$BlastTaxR};

		print "Interleaved LCA..\n";
		if ($doInter && $contInter){
			if ($DoPar){
				my @subf = splitBlastTax($taxblastf2,$BlastCores);
				#paralellize getTaxForOTUfromRefBlast
				my @thrs;
				for (my $i=0;$i<@subf;$i++){$thrs[$i] = threads->create(\&getTaxForOTUfromRefBlast,$subf[$i],\%GG,0,1);}
				#combine tax object
				for (my $i=0;$i<@subf;$i++){my $ret = $thrs[$i]->join();$BlastTaxRi = {%$BlastTaxRi,%$ret};}
			} else {
				$BlastTaxRi = getTaxForOTUfromRefBlast($taxblastf2,\%GG,1,$BlastCores);
			}
			$interLeave = $tmpD."Inter__U$DBi".".fna" if ($tmpD ne "");
			if ($DBi < (@DBa-1)){
				($interLeave,$BlastTaxRi,$leftover,$fasrI)= findUnassigned( $BlastTaxRi,$fasrI,$interLeave );
			} else {findUnassigned($BlastTaxRi,$fasrI,"");}
			$contInter=0 if ($leftover ==0 );

		}
		$fullBlastTaxRi = {%$fullBlastTaxRi,%$BlastTaxRi};

		undef %GG ; system "rm -f $taxblastf $taxblastf2" ; #die();
	}
	findUnassigned($fullBlastTaxR,{},"");
	findUnassigned($fullBlastTaxRi,{},"");

	system "rm -f $queryO"."__U*"." $interLeaveO"."__U*";
	$hof1 = writeBlastHiera($fullBlastTaxR,$id,$simName);  undef $BlastTaxR ;
	#merging of interleave results
	if ($doInter){
		$hof2 = writeBlastHiera($fullBlastTaxRi,$id."_inter",$simName); $BlastTaxRi = {};
		#system "cat $taxblastf2 >> $taxblastf; ";
		system "cat $hof2 >> $hof1; rm $hof2;";
	}
	
	system "rm -f $queryO"."__U* $interLeaveO"."__U*";
}

#file operations on paired end files
sub merge($ $ $){
	my ($r1,$r2,$interleave) = @_;
	print "Interleaving 2 files..\n";
	open I1,"<$r1" or die "1   $r1\n".$!;open I2,"<$r2"or die "2   $r2\n".$!;open O,">$interleave"or die "interl   ".$!;
	my $line1=<I1>;my $line2=<I2>; $line2="";#read header & get rid of it..
	while(1){
		my $tmp="";
		while ($tmp !~ m/^>/){chomp $tmp;$line1.=$tmp;$tmp=<I1>; last unless defined $tmp;}
		#print $line1."\n";
		print O $line1; if (defined $tmp ){chomp $tmp; $line1 = $tmp."\n"; }
		$tmp=""; 
		while ($tmp !~ m/^>/){chomp $tmp; $line2.=$tmp;$tmp=<I2>;last unless defined $tmp;}
		#get one more line to get rid of header:
		print O reverse_complement_IUPAC($line2)."\n"; $line2 = "";#$tmp;
		last unless defined $tmp;
	}
	close O; close I1; close I2;
#	my $mergeScript = "/g/bork5/hildebra/bin/sortmerna-2.0/scripts/merge-paired-reads.sh";
#	system "bash $mergeScript $r1 $r2 $interLeave";

}

#load tax database
sub getGGtaxo($){
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
	return %ret;
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
	my $LCAfraction = 0.9;
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

sub empty_proc(){return {};}
#parallize LCA helper
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

#function that goes line by line through Blast m8 table
sub getTaxForOTUfromRefBlast($ $ $ $){
	my ($blastout,$GGref,$interLMode,$BlastCores) = @_;
	#sp,ge,fa,or,cl,ph
	my %GG = %{$GGref};
	my @ggk = keys(%GG);
	my $maxGGdep = scalar(@{$GG{$ggk[0]}});
	@ggk = ();
	open B,"<",$blastout or die "Could not read $blastout\n";
	my $sotu = "";my $sID=0 ; my $sLength=0; 
	#my @sTax=(); my @sMaxTaxNum = ();
	my %retRef ;# my $retDRef ={};
	my $cnt=0;
	my @tmpLines = (); #stores Blast lines
	my @spl; #temp line delim
	my $minBit = 120; my $minEval = 1e-14;
	my %prevQueries = ();
	my $line = "";
	my @thrs; my $thrsCnt = 0;
	#for (my $i=0;$i<$BlastCores;$i++){$thrs[$i] = threads->create(\&empty_proc);}
	while ($line = <B>){
		$cnt++; chomp $line; #$line2 = $line;
		my @spl = split("\t",$line);
		my $totu = $spl[0]; #line otu
		$totu =~ s/^>//;
		if ($cnt == 1) {$sotu = $totu;}
		#print $line." XX $spl[11] $spl[10]\n"; die () if ($cnt > 50);
		
		#check if this is a 2 read-hit type of match (interleaved mode) & merge subsequently
		if ($interLMode){
			if (@tmpLines>0 && exists($prevQueries{ $spl[1]}) ){
				my @prevHit = @{$prevQueries{$spl[1]}};
				die "something went wrong with the inter matching: $prevHit[1] - $spl[1]\n" unless ( $prevHit[1] eq $spl[1] );
				$prevHit[11] += $spl[11];#bit score
				$prevHit[3] += $spl[3];$prevHit[5] += $spl[5];$prevHit[4] += $spl[4];#alignment length,mistmatches,gap openings
				$prevHit[2] = ($prevHit[2] + $spl[2]) / 2;
				#$tmpLines[-1] = \@prevHit;
				@spl = @prevHit;
				#next;
			} else {$prevQueries{ $spl[1] } = \@spl;}
		}
		if (($spl[11] < $minBit) || ($spl[10] > $minEval) ){ #just filter out..
			#print "ss\n"; 
			#next;
			$spl[2] =0; #simply deactivate this way...
		}
		
		if ($sotu eq $totu){
			push(@tmpLines,\@spl);
			if ($spl[2] > $sID && $spl[3] > $sLength){$sID = $spl[2]; $sLength = $spl[3];}
			if ($spl[3] > ($sLength*1.4) && $spl[2] > ($sID*0.9)) {$sID = $spl[2]; $sLength = $spl[3];} #longer alignment is worth it..
			#print $sID."\n";
		} else {
			#print "DD";
			#print "Maybe\n";
			
#			if (0 && $sotu ne ""){				if ($thrsCnt>= $BlastCores){$thrsCnt=0;}				#get result of old job
#				my $ret = $thrs[$thrsCnt]->join();				foreach (keys %{$ret}){$retRef{$_} = ${$ret}{$_};}				$thrs[$thrsCnt] = threads->create(\&rework_tmpLines,\@tmpLines,$sotu,$sID,$sLength,\%GG,$maxGGdep,$sMaxTax) ;
#				$thrsCnt++; print $thrsCnt." ";			}
			if ($sotu ne ""){
				my ($AR) = rework_tmpLines(\@tmpLines,$sID,$sLength,\%GG,$maxGGdep,); 
				$retRef{$sotu} = $AR;
			}
			$sotu = $totu; undef @tmpLines ; undef %prevQueries ;
			push(@tmpLines,\@spl);$prevQueries{ $spl[1] } = \@spl;
			$sID =  $spl[2]; $sLength = $spl[3];
		}
	}
	#last OTU in extra request
	if ($sotu ne ""){
		#my @spl = split("\t",$line);
		
		my ($AR) = rework_tmpLines(\@tmpLines,$sID,$sLength,\%GG,$maxGGdep);
		$retRef{$sotu} = $AR;
#		my($sTaxX,$sMaxTaxX) = LCA(\@sTax,\@sMaxTaxNum,$maxGGdep);
#		$ret{$sotu} = $sTaxX;
#		$retD{$sotu} = $sMaxTaxX;
	}
	close B;
	#debug 
	#my %ret = %{$retRef};	my @tmp = @{$ret{$sotu}};print "\n@tmp  $sotu\n";
	undef  @tmpLines;undef %prevQueries ;
	#print "Assigned $refDBname Taxonomy to OTU's\n",0;
	
	return (\%retRef);
}

#get a group of Blast lines that are considered "valid" hits and apply LCA to these
sub rework_tmpLines(){
	my ($tmpLinesAR,$sID,$sLength,$GGhr,$maxGGdep) = @_;
	my @tmpLines= @{$tmpLinesAR};
	my $debug_flag = 0;
	if (@tmpLines == 0 || $sID == 0){#prob no entry passed inclusion criteria
		return ([]);
	}
	
	#if ($sotu eq "S1__ITS_26"){print "\nYAAYYA\n\n\n";}
	my %GG = %{$GGhr};
	#my @ggk = keys %GG; print @ggk."NN\n";
	#my @spl = @{$splAR};
	my @sTax=(); my @sMaxTaxNum = ();
	my $tolerance = 1.5;
	#extend tolerance upon lesser hits (to include more spurious, off-target hits)
	if ($maxHitOnly==1){$tolerance=0;}	elsif ($sID == 100){$tolerance=0.1;}
	elsif ($sID >= 99.5){$tolerance=0.25;}	elsif ($sID >= 99){$tolerance=0.5;}
	elsif ($sID >= 98){$tolerance=1;}	elsif ($sID >= 97){$tolerance=1.25;}
	#print "XX $sID $tolerance $sLength\n";
	foreach my $lin2 (@tmpLines){#just compare if the tax gets any better
		#only hits within 1.5% range of best (first) hit considered
		my @spl2 = @{$lin2};#split("\t",$lin2);
		#print "$spl2[2] < ($sID - $tolerance\n";
		if ($spl2[2] < ($sID - $tolerance)) {next;}
		if ($spl2[3] < ($sLength * $lengthTolerance)){next;}
		my $sMax2 = 0;
		foreach (@idThr) {if($spl2[2] < $_){$sMax2++}};
		$sMax2 = 7 - $sMax2;
		unless (exists $GG{$spl2[1]} ){die "Can't find GG entry for $spl2[1]\n";}
		my $tTax = $GG{$spl2[1]} ;
		#print $tTax." JJ\n";
		my $sMax3 = maxTax($tTax);
		#die "$tTax  $sMax3\n";
		if ($sMax3 <= $sMax2){$sMax2 = $sMax3;
		} 
		push(@sTax,$tTax);
		push(@sMaxTaxNum,$sMax2);
	}
	#print "@sTax\n";
	#print $sID."\n";
	#entry for last OTU with best results etc..
	die "sTax not defined: LC=".@tmpLines."\n@{$tmpLines[0]}\n@{$tmpLines[1]}\n@{$tmpLines[2]}\n@{$tmpLines[3]}\n" unless ( @sTax > 0);
	my($sTaxX) = LCA(\@sTax,\@sMaxTaxNum,$maxGGdep);
	return ($sTaxX);#,\%retD);

}


#purpose is to rewrite the sample names in the blast file to smplID and a count to make adding up later easier..
sub rewriteBlastSmplNm($ $ $){
	my ($blastfile,$smplID,$MKtype) = @_;
	my $cnt = 0;
	my $blastfile2 = $blastfile.".tmp";
	open I,"<$blastfile";  open O,">$blastfile2"; open L,">$blastfile.transcribe.log";
	while (my $l = <I>){
		chomp $l;
		my @spl = 
		my $rdNm = "$smplID"."__$MKtype$cnt";
		print O $l."\n";
		$cnt++;
	}
	close I; close O; close L;
	system "rm $blastfile;mv $blastfile2 $blastfile";
}

#write Blast hierachy helper
sub blastFmt($ $){
	my ($aref,$tD) = @_;
	$aref = correctTaxString($aref,$tD);
	my @in  = @{$aref};
	return join("\t",@{$aref});
}

#write out result of LCA, while sorting in non hits etc, do some final formatting for consistency etc
sub writeBlastHiera($ $ $ ){
	my ($taxr,$id,$simTerm) = @_;
	my @avOTUs = keys %{$taxr};
	my $cnt = 0;
	my $outfile = "$outdir/$id.hiera.txt";
	open H,">",$outfile;#_$simTerm
	print H "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOTU\n";
	#
	foreach (@avOTUs){
		my $tdep = 0;
		$tdep = maxTax(${$taxr}{$_}) if ( exists( ${$taxr}{$_}) &&  ${$taxr}{$_} );
		if (  $tdep > 1 ){
			print H blastFmt(${$taxr}{$_},$tdep)  ."\t".$_."\n";
		} else {
			print H "?\t?\t?\t?\t?\t?\t?\t".$_."\n";
		}
		$cnt ++;
	}
	close H;
	print "$outfile\n";
	return $outfile;
}
