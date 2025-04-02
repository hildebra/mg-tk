#!/usr/bin/env perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw(systemW);
use Mods::Binning qw (runMetaBat runCheckM runSemiBin runMetaDecoder );
use Mods::IO_Tamoc_progs qw(getProgPaths jgi_depth_cmd);

#v0.1: 2.3.24: ini version FH
#v0.11: 23.4.24: added support for hybrid assemblies (i.e. 2 crams / sample)
my $version = 0.11;


my $DoMetaBat2 = "";
my $BinDir = "";
my $smplIDs1 = "";
my $nodeSpTmpD2 = "";
my $cAssGrp = "";
my $metaGassembly = "";
my $MB2coresL = 1;
my $pathsPre = "";
my $seqTec = "ill";
my $BinnerName = "";
my $giveSBenv = ""; #human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global


#"$BinnerScr -binner $DoMetaBat2 -binD $Bindir -smplID $smplIDs1 -tmpD $nodeSpTmpD2 -assmbl $metaGassembly -assmblGrp $cAssGrp -cores $MB2coresL -smplDirs " . join(",",@paths) . " -seqTec $seqTec 

GetOptions(
#Directories/files
	"binner=i"  => \$DoMetaBat2, #1: metabat2 2: semibin 3: metadecoder
	#"tmp=s" => \$tmpDir, 
	"binD=s" => \$BinDir, #dir where bin result will be stored
	"smplID=s" => \$smplIDs1, #id of current (assembled) sample
	"tmpD=s" => \$nodeSpTmpD2, #node specific (local) tmp dir
	"assmbl=s" => \$metaGassembly, #file path to assembly, fna
	"assmblGrp=s" => \$cAssGrp, #assembly group of current sample
	"cores=i" => \$MB2coresL , #num cores to be used (locally in this case)
	"smplDirs=s" => \$pathsPre, #dirs to MF output dirs, where the map etc will be looked up
	"seqTec=s" => \$seqTec  #PB, ONT, ill etc
);




my @paths = split /,/,$pathsPre;

print "======================================================================\n";
print "   runBinners.pl v $version\n";
print "     found " . scalar @paths . " sample dirs for sample \"$smplIDs1\"\n";
print "     using assembly $metaGassembly\n";
print "     using $MB2coresL cores, binner \"$DoMetaBat2\" to outdir $BinDir\n";
print "======================================================================\n";



my $MBcmd = "";#"mkdir -p $nodeSpTmpD2";
system "rm -rf $nodeSpTmpD2; mkdir -p $nodeSpTmpD2" unless (-d $nodeSpTmpD2);
system "rm -rf $BinDir; mkdir -p $BinDir;";

#actual work happens here..
if ($DoMetaBat2 == 1){
	$MBcmd .= jgi_depth_cmd(\@paths,$nodeSpTmpD2."/depth$cAssGrp",95,$MB2coresL,$metaGassembly);# unless (-e );
	my $tmp = runMetaBat("$nodeSpTmpD2/depth$cAssGrp.jgi.depth.txt",$BinDir,$smplIDs1,$metaGassembly,$MB2coresL);
	if ($tmp eq "" ){$MBcmd="";
	} else {$MBcmd .= $tmp;}
	$BinnerName = "MB2";
} elsif ($DoMetaBat2 == 2){#SemiBin
	#could work with jgidepth, but only for single sample.. not for co-assembly :(
	#"$nodeSpTmpD2/depth$cAssGrp.jgi.depth.txt"
	$MBcmd .= runSemiBin("",$BinDir,$nodeSpTmpD2, $smplIDs1,$metaGassembly,$MB2coresL,\@paths, $seqTec, $giveSBenv);
	$BinnerName = "SB";
}elsif ($DoMetaBat2 == 3){#MetaDecoder
	#could work with jgidepth, but only for single sample.. not for co-assembly :(
	#"$nodeSpTmpD2/depth$cAssGrp.jgi.depth.txt"
	$MBcmd .= runMetaDecoder("",$BinDir,$nodeSpTmpD2, $smplIDs1,$metaGassembly,$MB2coresL,\@paths);
	$BinnerName = "MD";
}

my $stone = "$BinDir/Binning.stone";
$MBcmd .= "\ntouch $stone\n";


print "running: $MBcmd\n";
#die ;
systemW $MBcmd;


print "Done executing binner $BinnerName\n";