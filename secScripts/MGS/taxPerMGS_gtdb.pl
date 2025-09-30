#!/usr/bin/perl
#uses MGS to cluster genes and kraken tax assignments
#./taxPerMGS.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//MB2.clusters.ext.can.Rhcl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/
use warnings; use strict;

use Mods::Binning qw(readMGSrev );
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::GenoMetaAss qw( systemW gzipopen);


#my $COND = getProgPaths("CONDA"); #source conda..
#my $py3activate = getProgPaths("py3activate"); #source conda.. 
my $GTDBtkBin = getProgPaths("GTDBtk");
my $GTDBtkDB = getProgPaths("GTDBtk_DB");
my $GTDBtkMash = getProgPaths("GTDBtk_mash",0);

#die "$GTDBtkDB\n$GTDBtkMash\n";
my $refMGd = $ARGV[0];
my $ncore = $ARGV[1];
my $tmpD = $ARGV[2];
my $Bdir = $ARGV[3];

my $pplacer_cores = $ncore;
$pplacer_cores = 2 if ($ncore > 2);

$tmpD.="/GTDB/";
my $oDir = "$tmpD/GTDBTK/";
system "mkdir -p $oDir" unless (-d $oDir);
system "mkdir -p $tmpD" unless (-d $tmpD);

my $cmd = "";
#$cmd .= "$COND\n$py3activate\n";
#--scratch_dir $tmpD 
#$cmd .= "export GTDBTK_DATA_PATH=$GTDBtkDB\n";
if ($GTDBtkBin =~ m/ activate /){
	$GTDBtkBin =~ s/activate (\S+)/activate $1\nexport GTDBTK_DATA_PATH=$GTDBtkDB\n/;
}

my $mashArg="" ;
if ($GTDBtkMash ne ""){ #for newer GTDBtk versions not supported
	$mashArg = "--mash_db $GTDBtkMash/\$MVERSION/";
	$cmd .= "MVERSION=`mash --version`"
}

$cmd .= "$GTDBtkBin classify_wf -x fna $mashArg --cpus $ncore --pplacer_cpus $pplacer_cores  --genome_dir $refMGd --out_dir $oDir"; #--scratch_dir $tmpD/GTtmp/ --genes
print "\n\n".$cmd."\n\n";
#die;
systemW $cmd;

systemW "rm -f $Bdir/gtdbtk.summary.tsv";
system "cat $oDir/gtdbtk.bac*.summary.tsv >> $Bdir/gtdbtk.summary.tsv\n";# if (-e "$oDir/gtdbtk.bac120.summary.tsv");
system "cat $oDir/gtdbtk.ar*.summary.tsv >> $Bdir/gtdbtk.summary.tsv\n";# if (-e "$oDir/gtdbtk.ar122.summary.tsv");
systemW "cut -f1,2 $Bdir/gtdbtk.summary.tsv | sed 's/.__//g' > $Bdir/GTDBTK.tax";
systemW "tar -zcvf $refMGd/GTDBtk.tar.gz $oDir";
systemW "rm -r $oDir";
systemW "rm -r $tmpD";

#transfer files
