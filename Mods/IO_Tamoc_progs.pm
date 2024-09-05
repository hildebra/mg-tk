package Mods::IO_Tamoc_progs;
use warnings;
use Cwd 'abs_path';
use strict;

use vars qw($CONFIG_FILE @CONFIG_TEXT %CONFIG_HASH);
$CONFIG_FILE="";
@CONFIG_TEXT = ();
%CONFIG_HASH = ();
sub setConfigFile;

#TAMOC programs related to IO to other programs, program paths .. not real subroutines that do anything

use Exporter qw(import);
our @EXPORT_OK = qw(getProgPaths 
					inputFmtSpades inputFmtMegahit jgi_depth_cmd createGapFillopt setConfigFile 
					buildMapperIdx mapperDBbuilt decideMapper checkMapsDoneSH greaterComputeSpace convert2Gb);


sub decideMapper($ $){
	my ($MapperProg,$readTec) = @_;
	if ($MapperProg == -1 ) {
		if ($readTec eq "PB" || $readTec eq "ONT"){  #default to minimap2 for 3rd gen seq
			$MapperProg = 3;
		} else {#otherwise use bowtie2
			$MapperProg = 1;
		}
	}
	if ($MapperProg<1 || $MapperProg>5){
		die "IO_Tamoc_progs.pm::decideMapper:: unknown MapperProg provided: $MapperProg\n! Aborting..\n";
	}
	return $MapperProg;
}

sub checkMapsDoneSH{
	my ($inAR) = @_;
	my @dirSS = @{$inAR};
	my $ctrlStr = "";
	foreach my $DDI (@dirSS){
		my $iBAM = $DDI;
		if ( $DDI =~ m/\/$/  ){
			$ctrlStr .= "if [ ! -e $DDI/mapping/done.sto ]; then echo \"Can't find $DDI/mapping/done.sto !! Aborting .. \"; exit 1; fi \n"; 
		} elsif ($DDI !~ m/\.bam$/ || $DDI !~ m/\.cram$/)  {
			$ctrlStr .= "if [ ! -e $DDI ]; then echo \"Can't find $DDI !! Aborting ..\"; exit 1; fi \n"; 
		}
	}
	return $ctrlStr;
}

#computes which string in the style of "120G" or "120" is greatest and returns this (in "G" format)
sub greaterComputeSpace{
	my @inA  = @_;
	my $ret = 0; #unit is G (gigabyte)
	foreach my $in (@inA){
		if ($in =~ m/\D/){
			#print "NON\n";
			if ($in =~ /^\d+G$/){
				$in =~ s/G//;$in = $in + 0;
			}elsif ($in =~ /^\d+M$/){
				$in =~ s/M//; $in = 0+$in; $in /= 1024;
			}elsif ($in =~ /^\d+T$/){
				$in =~ s/T//; $in = 0+$in; $in *= 1024;
			} else {
				die "Unrecognized compute space: $in\n";
			}
		} else{
			$in = 0+$in; #convert to num
		}
		$ret = $in if ($in > $ret);
	}
	#$ret .= "G";
	#die "$ret\n";
	return $ret;
}

#converts string of style "120G" or "120" to GB. If no "G/M/T" given, assumes already in Gb ("G")
sub convert2Gb($){
	my ($tmpSpace) = @_;
	if ($tmpSpace !~ m/\D/ && ($tmpSpace eq "0" || ($tmpSpace+0) == 0) ){$tmpSpace = 0 ;
	}elsif ($tmpSpace =~ s/G$//){$tmpSpace = int($tmpSpace+0.5);
	}elsif ($tmpSpace =~ s/T$//){$tmpSpace *= 1024 ;
	}elsif ($tmpSpace =~ s/M$//){$tmpSpace /= 1024 ;
	} else {$tmpSpace = int($tmpSpace+0.5);}
	return $tmpSpace;
}

sub setConfigFile{
	my @var = @_;
	my $customCfg = 0;
	if (@var == 1 && $var[0] eq "internal"){
		my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"};
		$modDir =~ s/IO_Tamoc_progs.pm//;
		$CONFIG_FILE = "$modDir/../Mods/config_internal.txt";
	} elsif (@var == 1 && $var[0] eq "DBconfig"){
		my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"};
		$modDir =~ s/IO_Tamoc_progs.pm//;
		$CONFIG_FILE = "$modDir/../Mods/config_DBs.txt";
	} elsif (@var == 1 && $var[0] ne ""){
		$CONFIG_FILE = $var[0];
		$customCfg = 1;
	} else {#default value
		my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"};
		$modDir =~ s/IO_Tamoc_progs.pm//;
		$CONFIG_FILE = "$modDir/MATAFILERcfg.txt";
	}
	die "Can't find MATAFILER config file: $CONFIG_FILE\nConsider changing path to config file via \"-config\" argument.\n Aborting..\n" unless (-e $CONFIG_FILE);
	print "Using config file : $CONFIG_FILE\n" if ($customCfg);
}

sub truePath{
	my ($TMCpath) = @_;
	if ($TMCpath =~ m/^\$/){
	$TMCpath =~ s/^\$//; 
	$TMCpath = $ENV{$TMCpath};
	}
	return $TMCpath;

}

sub loadConfigs{
	#loads once in every program run the entire config file(s) into hash %CONFIG_HASH
	if (scalar @CONFIG_TEXT == 0){
		print "READING config files..";
		setConfigFile() if ($CONFIG_FILE eq "");
		open I,"<$CONFIG_FILE" or die "Can't open $CONFIG_FILE\n";
		chomp(@CONFIG_TEXT = <I>);
		close I;
		setConfigFile("internal") ;
		open I,"<$CONFIG_FILE" or die "Can't open internal $CONFIG_FILE\n";
		my @INTtmp;
		chomp(@INTtmp = <I>);
		close I;
		push(@CONFIG_TEXT,@INTtmp);
		#DB config read..
		setConfigFile("DBconfig") ;
		open I,"<$CONFIG_FILE" or die "Can't open DBconfig $CONFIG_FILE\n";
		@INTtmp=();
		chomp(@INTtmp = <I>);
		close I;
		push(@CONFIG_TEXT,@INTtmp);
	}
	#my $condaA = getProgPaths("CONDA");
	#die "@CONFIG_TEXT\n";
	print "converting config files.. ";
	my $TMCpath = "";my $Tset=0; my $BINpath = "";my $Bset=0; 
	my $DBpath = "";my $DBset=0; my $SINGcmd = ""; my $SINGset=0; 
	my $CONDset = 0; my $CONDcmd="";my $CONDset2 = 0; my $CONDA="";my $CONDset3 = 0; my $CONDAbaseEnv="MFF";
	my $PY3cmd = ""; my $PY3set=0; my $Rscriptcmd = ""; my $Rscriptset=0;
	foreach my $l (@CONFIG_TEXT){
		next if ($l =~ m/^#/ || length($l) == 0);
		if (!$Tset && $l =~ m/^MFLRDir\t([^#]+)/){
			$Tset=1;$TMCpath = truePath($1);
#			next;
		} elsif (!$Bset && $l =~ m/^BINDir\t([^#]+)/){
			$BINpath = truePath($1); $Bset=1;
		} elsif (!$DBset && $l =~ m/^DBDir\t([^#]+)/){
			$DBpath = truePath($1); $DBset=1;
		} elsif (!$SINGset && $l =~ m/^SINGcmd\t([^#]+)/){
			$SINGcmd = truePath($1); $SINGset=1;
		} elsif (!$CONDset && $l =~ m/^CONDcmd\t([^#]+)/){
			$CONDcmd = truePath($1); $CONDset=1;
		} elsif (!$CONDset2 && $l =~ m/^CONDA\t([^#]+)/){
			$CONDA = truePath($1); $CONDset2=1;
		} elsif (!$CONDset3 && $l =~ m/^CONDAbaseEnv\t([^#]+)/){
			$CONDAbaseEnv = truePath($1); $CONDset3=1;
			
		} elsif (!$PY3set && $l =~ m/^PY3cmd\t([^#]+)/){
			$PY3cmd = $1; $PY3set=1;
			if ($SINGset){$PY3cmd =~ s/\[SINGcmd\]/$SINGcmd/;}
		} elsif (!$Rscriptset && $l =~ m/^Rscript\t([^#]+)/){
			$Rscriptcmd = $1; $Rscriptset=1;
			#print "Rscript set!! : $Rscriptcmd  \n\n\n";
		
		
		} else { #bit strange way of doing this.. but compliant with old style
			my $XVar = "";
			#print "$l\n";
			my @spl = split (/\t/,$l);
			$XVar = $spl[0];
			if (@spl == 1) {$CONFIG_HASH{$XVar} = "";next;}
			
			$l =~ m/^$XVar\t([^#^\t]+)/;
			my $reV = $1;
			
			#die "$reV  $XVar  $l\n";
			
			$reV =~s/\[MFLRDir\]/$TMCpath/;
			$reV =~ s/\[BINDir\]/$BINpath/;
			$reV =~ s/\[DBDir\]/$DBpath/;
			$reV =~ s/\[SINGcmd\]/$SINGcmd/;
			$reV =~ s/\[PY3\]/$PY3cmd/;
			$reV =~ s/\[Rscript\]/$Rscriptcmd/;
			if ($l =~ m/env:([^#^\t]+)/){
				
				$reV = "$CONDA;$CONDcmd activate $1\n$reV";
			}
			
			#return $reV;
			$CONFIG_HASH{$XVar} = $reV;
		}
	}
	#some check about basic params being set..
	if (!$Rscriptset){die "Could not find \"Rscript\" correctly configured in config file, please check your local config! Aborting..\n";}
	if (!$CONDset){die "Could not find \"CONDcmd\" correctly configured in config file, please check your local config! Aborting..\n";}
	if (!$DBset){die "Could not find \"DBDir\" correctly configured in config file, please check your local config! Aborting..\n";}
	if (!$Tset){die "Could not find \"MFLRDir\" correctly configured in config file, please check your local config! Aborting..\n";}
	$CONFIG_HASH{"activateBase"} = "$CONDA;$CONDcmd activate $CONDAbaseEnv\n";
	$CONFIG_HASH{"CONDAbaseEnv"} = $CONDAbaseEnv;
	$CONFIG_HASH{"MFLRDir"} = $TMCpath;
	$CONFIG_HASH{"BINDir"} = $BINpath;
	$CONFIG_HASH{"DBDir"} = $DBpath;
	$CONFIG_HASH{"Rscript"} = $Rscriptcmd;
	print "  Done. ";
}


sub getProgPaths{
	my @var = @_;
	my $srchVar = $var[0] ;
	my $required=1;
	if (@var > 1){$required = $var[1];}
	#die "$required\n";

	my @multVars = ();
	if (ref $srchVar eq 'ARRAY') {
		#print "ARRAY\n";
		@multVars = @{$srchVar};
	}
	
	
	if (scalar(keys %CONFIG_HASH) == 0){
		#read in config hash _once_
		loadConfigs();
	}
	if (scalar(keys %CONFIG_HASH) == 0){
		die "Something went wrong loading MATAFILER configs.. aborting\n";
	}
	
	
	if (@multVars > 0){
		my @retA;
		for (my$i=0;$i<scalar(@multVars);$i++){if (exists($CONFIG_HASH{$multVars[$i]})) { $retA[$i] = $CONFIG_HASH{$multVars[$i]};}} 
		return \@retA;
	}
	if (exists($CONFIG_HASH{$srchVar})){
		return $CONFIG_HASH{$srchVar};
	} else {
		die "Can't find configuration for $srchVar in MATAFILER config ($CONFIG_FILE)\n" if ($required != 0);
	}
	
	return "";
}

sub activateBase{
	die "activateBase:: not used";
	my $ret = "";
	#my $condaA = getProgPaths("CONDA");
	#$ret = "$condaA\n$CONDcmd activate $1\n$reV";
	return $ret;
}


sub mapperDBbuilt( $ $){
	my ($DBbtRef, $MapperProg2) = @_;
	my $bwt2IdxFileSuffix = ".bw2";my $mini2IdxFileSuffix = ".mmi";
	my $kmaIdxFileSuffix = ".kma";
	if ($MapperProg2 == 5){return 1;} #strobealign doesn't need index..
	#print "($MapperProg2 == 1 || $MapperProg2 == -1) && !-s $DBbtRef$bwt2IdxFileSuffix.rev.2.bt2\n";
	if ( 
		($MapperProg2 ==0 && !-e "$DBbtRef$bwt2IdxFileSuffix.0.sa") 
		|| ( ($MapperProg2 == 1 || $MapperProg2 == -1) && (!-s "$DBbtRef$bwt2IdxFileSuffix.rev.2.bt2" || !-s "$DBbtRef$bwt2IdxFileSuffix.1.bt2" ) ) #bowtie2
		||( $MapperProg2 == 2 && !-s "$DBbtRef.pac" ) #bwa
		||( ($MapperProg2 == 3 || $MapperProg2 == -1 ) && !-s "$DBbtRef$mini2IdxFileSuffix" ) #minimap2
		||( ($MapperProg2 == 4 ) && !-s "$DBbtRef$kmaIdxFileSuffix.seq.b" )#kma
	) {
		return 0;
	}
	return 1;
}

sub buildMapperIdx($ $ $ $){
	my ($REF,$ncore,$lrgDB,$MapperProg) = @_;
	#1=bowtie2, 2=bwa, 3=minimap2
	if ($MapperProg == 5){return ("",$REF,$REF);} #strobealign doesn't need index..
	my $bwt2IdxFileSuffix = ".bw2";my $mini2IdxFileSuffix = ".mmi";
	my $kmaIdxFileSuffix = ".kma";
	my $bwtIdx = $REF.$bwt2IdxFileSuffix;
	my $chkFi = $bwtIdx;
	$MapperProg = decideMapper($MapperProg,"");
	if ($MapperProg==1){$chkFi .=".1.bt2";
	}elsif ($MapperProg==2){$chkFi .= $REF.".pac";
	} elsif ($MapperProg == 3){$chkFi = $REF.$mini2IdxFileSuffix;
	} elsif ($MapperProg == 4){$chkFi = $REF.$kmaIdxFileSuffix.".seq.b";
	}
	my $dbCmd ="";
	$dbCmd .= "if [ ! -s $chkFi ];then \n";
	$dbCmd .= "echo \"Building index for mapper $MapperProg\"\n";
	if ($MapperProg==1){
		my $bwt2Bin = getProgPaths("bwt2");
		$dbCmd .= $bwt2Bin."-build ";
		$dbCmd .= " --large-index "if ($lrgDB);
		$dbCmd .= "-q $REF --threads $ncore $bwtIdx\n";
		if (-s $REF."$bwt2IdxFileSuffix.1.bt2" || -s $REF."$bwt2IdxFileSuffix.1.bt2l"){$dbCmd = "";} #deactivate if already built
	} elsif($MapperProg==2) { 
		my $bwaBin  = getProgPaths("bwa");
		$dbCmd .= $bwaBin." index $REF\n";
		if (-s $REF.".pac"){$dbCmd = "";} 
		#die $dbCmd."\n";
	} elsif ($MapperProg == 3){
		$bwtIdx = $REF.$mini2IdxFileSuffix;
		my $mini2  = getProgPaths("minimap2");
		$dbCmd .= "$mini2 -t $ncore -H -d $bwtIdx $REF\n";
		$dbCmd = "" if (-e $bwtIdx);
	} elsif ($MapperProg==4){ 			
		$bwtIdx = $REF.$kmaIdxFileSuffix;
		my $kmaBin = getProgPaths("kma");
		$dbCmd .= "$kmaBin index -i $REF -o $bwtIdx 2>/dev/null \n"; #-t $ncore 
	}

	$dbCmd .= "fi\n" unless ($dbCmd eq "");
	#die "$dbCmd\n";
	#my $jobN = "_ASDB$JNUM"; my $tmpCmd;
	#($jobN, $tmpCmd) = qsubSystem($logDir."BAM2CRAMxtra.sh",$dbCmd,1,"10G",$jobN,"","",1,[],\%QSBopt);
	return ($dbCmd,$bwtIdx,$chkFi);
}


sub inputFmtSpades($ $ $ $ $){
	my ($p1ar,$p2ar,$singlAr,$logDir,$cReadTecAr) = @_;
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	my @singl = @{$singlAr};
	my @readTec = @{$cReadTecAr};
	my $doYAML = 0;
	if (@p1 > 9){$doYAML=1;}
	if (@p1 != @p2){print "Unequal paired read array lengths arrays for Spades\n"; exit(2);}
	my $sprds = "";
	if ($doYAML==0){
		for (my $i =0; $i<@p1;$i++){
			my $peTerm = "--pe";$peTerm = "--gemcode" if ($readTec[$i] =~ m/SLR/);
			$sprds .= " ${peTerm}".($i+1) ."-1 $p1[$i] ${peTerm}".($i+1) ."-2 $p2[$i]";
		}
		for (my $i=0;$i<@singl;$i++){
			my $peTerm = "--pe";$peTerm = "--gemcode" if ($readTec[$i] =~ m/SLR/);
			$sprds .= " ${peTerm}".($i+1) ."-s $singl[$i]";
		}
	} else {
		open O,">$logDir/spadesInput.yaml";
		print O "  [\n      {\n        orientation: \"fr\",\n        type: \"paired-end\",\n        left reads: [\n";
		for (my $i =0; $i<@p2;$i++){
			if (($i+1) == @p2){	print O "          \"$p1[$i]\"\n";	} else { print O "          \"$p1[$i]\",\n";}
		}
		print O "     ],\n        right reads: [\n";
		for (my $i =0; $i<@p1;$i++){
			if (($i+1) == @p1){	print O "          \"$p2[$i]\"\n";	} else { print O "          \"$p2[$i]\",\n";}
			#print O "          \"$p2[$i]\"\n";
		}
		print O "       ]\n      }";
		if (@singl > 0){
			print O ",\n      {\n        type: \"single\",\n        single reads: [\n";
			for (my $i=0;$i<@singl;$i++){
				if (($i+1) == @singl){	print O "          \"$singl[$i]\"\n";	} else { print O "          \"$singl[$i]\",\n";}
				#print O  "          \"$singl[$i]\"\n";
			}
			print O "       ]\n      }";
		}
		#end of yaml file
		print O "\n     ]\n";
		close O;
		$sprds = " --dataset $logDir/spadesInput.yaml";
	}
	#die "$sprds\n\ninputFmtSpades\n";
	return $sprds;
 }

 sub inputFmtMegahit($ $ $ $){
	my ($p1ar,$p2ar,$singlAr,$logDir) = @_;
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	my @singl = @{$singlAr};
	if (@p1 != @p2){print "Unequal paired read array lengths arrays for Spades\n"; exit(2);}
	my $sprds = "";
	if (@p1 > 0){ 
		$sprds .= "-1 ".join(",",@p1) . " -2 ".join(",",@p2)." ";
	}
	if (@singl > 0){
		$sprds .= "-r ".join(",",@singl);
	}
	return $sprds;
 }


sub jgi_depth_cmd{
	my $dirsAR = $_[0];
	my $out = $_[1];
	my $perID = $_[2];#) = @_;
	my $numCores = 1;
	$numCores = $_[3] if (@_ > 3);
	my $refFA = "";
	$refFA = $_[4] if (@_ > 4);
	#die $dirs."\n";
	my $smtBin = getProgPaths("samtools");
	my $jgiScr = getProgPaths("jgiDepth");

	my @dirSS = @{$dirsAR};#split(',',$dirs);
	#go through each dir and find sample name
	my $comBAM = "";
	my $isCram=0;
	foreach my $DDI (@dirSS){
		if ( $DDI =~ m/\/$/ ||  $DDI !~ m/bam$/ ){
			my $SmplNm = `cat $DDI/mapping/done.sto`;#$SmplNm =~ s/-smd.bam\n?//;
			chomp $SmplNm;
			my $tbam = "$DDI/mapping/$SmplNm";
			if (!-e $tbam){$isCram=1;$tbam =~ s/\.bam/\.cram/;}
			die "jgi_depth_cmd:::Can't find either bam nor cram at $DDI\n" unless (-e $tbam);
			$comBAM .= "$tbam ";
		} else {
			$comBAM .= "$DDI ";
		}
	}
	#my $comBAM = join("/mapping/Align_ment-smd.bam ",@dirSS);
	my $covCmd = "";
	$covCmd .= "rm -f $out.jgi.*\n";
	if ($isCram){
		die "jgi_depth_cmd:::No reference Fasta given for @dirSS\n" if ($refFA eq "");
		if (@dirSS == 1){
			$covCmd .= "$smtBin view -T $refFA -@ $numCores -b $comBAM | ";
			$comBAM = "-";
		} else {
			my @splSS = split /\s/,$comBAM;
			$comBAM="";
			for (my $i=0;$i<@splSS;$i++){
				my $tmpBam = "$out.jgi.tmp.$i.bam";
				$covCmd .= "$smtBin view -T $refFA -@ $numCores -b $splSS[$i] > $tmpBam\n";
				$comBAM .= "$tmpBam ";
			}
		}
	}
	$covCmd .= $jgiScr;#"/g/bork3/home/hildebra/bin/metabat/./jgi_summarize_bam_contig_depths";
	#--pairedContigs $out.jgi.pairs.sparse
	$covCmd .= " --outputDepth $out.jgi.depth.txt  --percentIdentity $perID  $comBAM\n";
	#$covCmd .= "gzip $out.jgi*\n";
	if (-e "$out.jgi.depth.txt"){$covCmd="";}

	#$covCmd .= "gzip $nxtBAM.jgi*\n";
	return $covCmd;
}

#2nd: arrray of files, paired sep by ","
sub createGapFillopt($ $ $){
 my ($ofile, $arFiles, $insertSizAr) = @_;
 my @Files = @{$arFiles};
 my @insertSiz = @{$insertSizAr};
 my $opt = "";
 for (my $cnt=0;$cnt<@Files; $cnt++){
	my $line = "LIB$cnt bwa ";
	my $lineMode = 2; 
	if ($lineMode==2){ #PE; only available option at the moment
		my @curFils = split(",",$Files[$cnt]);
		die ("Only paired reads accepatble to GapFiler.\n") if (@curFils != 2);
		$line .= $curFils[0]." ".$curFils[1];
	}
	$line .= " $insertSiz[$cnt] 0.3 FR\n";
	$opt .= $line;
 }
 #print $opt."\n";
 open O,">",$ofile; print O $opt; close O;
}
