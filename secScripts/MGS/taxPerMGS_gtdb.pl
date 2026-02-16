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



if ($GTDBtkBin =~ m/ activate /){
	$GTDBtkBin =~ s/activate (\S+)/activate $1\nexport GTDBTK_DATA_PATH=$GTDBtkDB\n/;
}


#get GTDBtk version
my $verSt = `$GTDBtkBin --version`;
$verSt =~ m/version (\S+) /;
my $GTDBver = $1+0.0;

#print "$verSt\n\n$1\n";

print "Using GTDBtk ver $GTDBver\n";
#print "$GTDBtkBin --version\n";
#die;

my $cmd = "";
#$cmd .= "$COND\n$py3activate\n";
#--scratch_dir $tmpD 
#$cmd .= "export GTDBTK_DATA_PATH=$GTDBtkDB\n";

my $mashArg="" ;my $hook = "";
if ($GTDBtkMash ne "" && $GTDBver >= 2.1){ #for newer GTDBtk versions not supported
	$mashArg = "--mash_db $GTDBtkMash/\$MVERSION/";
	$hook = "MVERSION=`mash --version`;\n"
}

$cmd .= "$GTDBtkBin classify_wf -x fna $mashArg --cpus $ncore --pplacer_cpus $pplacer_cores  --genome_dir $refMGd --out_dir $oDir"; #--scratch_dir $tmpD/GTtmp/ --genes

#get hook inserted before command
$cmd =~ s/gtdbtk classify/${hook}gtdbtk classify/;

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
