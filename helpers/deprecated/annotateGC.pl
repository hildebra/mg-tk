#!/usr/bin/perl
#usage
#./annotateGC.pl /g/scb/bork/hildebra/SNP/GNMass2_singl/
#./annotateGC.pl /g/scb/bork/hildebra/SNP/GNMassSimuC/

my $inD = $ARGV[0];
my $GCd = "$inD/GeneCatalog/";


#KEGG KOs
my $hmmBin = "/g/bork5/hildebra/bin/hmmer-3.0/hmm30/bin/hmmsearch";
my $N = 30;
my $FOAMhmm = "/g/bork5/hildebra/DB/FOAM/FOAM-hmm_rel1.hmm";
my $query = "$GCd/compl.incompl.95.prot.faa";
system "mkdir -p $GCd/tmp" unless (-d "$GCd/tmp");
my $tmpOut = "$GCd/tmp/FOAM.hmm.dom";
my $outF = "$GCd/KOassignment.txt";
my $cmd = "$hmmBin --cpu $N -E 1e-05 --noali --domtblout $tmpOut $FOAMhmm $query > 1& /dev/null\n";
$cmd .= "sort $tmpOut > $tmpOut.sort\n";
$cmd .= "python /g/bork3/home/hildebra/DB/FOAM/scripts/bmn-HMMerBestHit.py $tmpOut.sort > $tmpOut.sort.BH\n";
$cmd .= "awk '{print $1,$4}' $tmpOut.sort.BH > $outF\n";
$cmd .= "rm -f $tmpOut $tmpOut.sort";
die $cmd."\n";
#sort Sample1_hmmsearch.out > Sample1_hmmsearch.out.sort
#python bmn-HMMerBestHit.py Sample1_hmmsearch.out.sort > Sample1_hmmsearch.BH
#awk '{print $4}' Sample1_hmmsearch.BH > Sample1_hmmsearch.BH.tmp1