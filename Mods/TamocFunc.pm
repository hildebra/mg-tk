package Mods::TamocFunc;
use warnings;
#use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 

#use Mods::GenoMetaAss qw(qsubSystem);

use Exporter qw(import);
our @EXPORT_OK = qw( sortgzblast uniq  readTabbed readTabbed2 bam2cram cram2bsam checkMF
					readTabbed3 readTable getSpecificDBpaths  getFileStr displayPOTUS 
					checkMFFInstall
					);
use Mods::GenoMetaAss qw(systemW readFasta renameFastHD gzipwrite gzipopen filsizeMB);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::Subm qw ( qsubSystem   );



#gets required FA genes via samtools faidx and saves them to file
#Usage: [text file with target genes] [out faa] [in faa(search genes in this fna)] [new names for output]



sub checkMFProgVers{
	my $progSet = 1;
	$progSet = $_[0] if (@_ >= 1);
	#print "Checking essential binaries ..";
	my @progs;
	if ($progSet == 1){
		@progs = ("sdm","LCA","readCov");
	} elsif ($progSet == 2){
		@progs = ("LCA","rare","clusterMAGs","canopy","MSAfix");
	} else {
		die "TamocFunc:::checkMFProgVers:: unknown \$progSet var: $progSet\n Aborting..\n";
	}
	foreach my $prog (@progs){
		#print "$prog\n";
		my $binary = getProgPaths($prog,1);
		die "Could not find $prog! Aborting\n" if ($binary eq "");
		my $cmd = "$binary -v  > /dev/null";
		#my $cmd = "$binary -v 2&>1"; #my $ver = `$cmd `;
		#if ($ver !~ m/.*$binary /){
		if (system($cmd)){die "Program $prog ($binary) failed!\nPlease check that the program is executable before continuing to run MATAFILER\n";}
		
	}
	
	#print "looks good!\n";
	
}


sub checkVersion($ $){
	my ($progN, $progI) = @_;
	
	#print "Checking $progN..";
	my $bin1 = getProgPaths($progI);
	
	my $versS = `$bin1 -v`;
	if ($versS =~ m/version ([\.\d]+)/){
		return $1;
	} elsif ($versS =~ m/([\.\d]+)/){
		return $1;
	}
	#print "$versS\n\n$progN found $1 \n";
	return "";
}

sub checkProg($ $ $){
	my ($progN, $progI, $doVer) = @_;
	
	print "Checking $progN..";
	my $bin1 = getProgPaths($progI);
	print "$bin1\n";
	my @bin = split /\n/,$bin1;
	my $idxA = 0; my $binA="";
	for (my $i=0;$i<@bin;$i++){if ($bin[$i] =~ m/micromamba\s*activate\s/){
			$idxA=$i;$binA=`$bin[$idxA];which $bin[$idxA+1]`;
			#die "\n\n".$binA."\n\n\n";
			$idxA++;last;
		} 
	}
	#print "@bin\n\nwhich $bin[$idxA]\n$idxA\n";
	$binA = `which $bin[$idxA]` if ($idxA==0);
	unless ($binA){
		die "Fatal: Could not excectute: $bin[$idxA] ! Please check MATAFILER install!\nExiting MATAFILER\n";
	}
	my $ver = "";
	if ($doVer){$ver = checkVersion($progN, $progI); $ver = "ver $ver"; }
	print "   $ver Present and excecutable\n";
}


sub checkMFFInstall{
	my ($instSto,$doExit) = @_;
	print "\n";
	print "Checking essential programs for presence and excecutability\n\n---------------------------------------------------\n";
	checkProg("simple demultiplexer (sdm)","sdm",1);
	#checkProg("neg test","noneeneoneono");
	checkProg("read Coverage estimator","readCov",1);
	checkProg("MSAfix","MSAfix",1);
	checkProg("clusterMAGs","clusterMAGs",1);
	checkProg("rarefaction toolkit2 (rtk2)","rare",1);
	checkProg("least common ancestor (LCA)","LCA",1);
	#env MGTKbinners
	checkProg("canopy clustering genome binner","canopy",0);
	
	checkProg("pigz","pigz",0);
	checkProg("bzip2","bzip2",0);
	checkProg("rsync","rsync",0);

	
	#mappers
	checkProg("bowtie2 mapper","bwt2",0);
	checkProg("strobealign mapper","strobealign",0);
	checkProg("foldseek mapper","foldseek",0);
	checkProg("diamond mapper","diamond",0);
	checkProg("lambda mapper","lambda",0);
	checkProg("minimap2 mapper","minimap2",0);
	checkProg("bwa mapper","bwa",0);
	#checkProg("vsearch mapper","vsearch",1);
	checkProg("hmmsearch","hmmsearch",0);
	
	checkProg("samtools","samtools",0);
	
	#assemblers
	checkProg("pprodigal gene prediction","pprodigal",0);
	checkProg("spades assembler","spades",1);
	checkProg("megahit assembler","megahit",1);
	checkProg("flye assembler","flye",1);
	checkProg("metaMDBG assembler","metaMDBG",1);
	
	checkProg("mmseqs2 seq clustering","mmseqs2",1);
	
	#clustering
	checkProg("cdhit seq clustering","cdhit",1);
	
	

	
	#binners
	checkProg("metabat2 genome binner","metabat2",0);
	checkProg("SemiBin2 genome binner","SemiBin2",1);
	#checkProg("MetaDecoder","MetaDecoder",1);
	
	checkProg("checkm2 MAG qual","checkm2",0);
	checkProg("GTDBtk profiler MAGs","GTDBtk",1);
	#checkVersion("GTDBtk profiler MAGs","GTDBtk",1);
	checkProg("MetaPhlan profiler","metPhl2",1);
	checkProg("mOTUs profiler","motus2",0);
	
	
	

	checkProg("clustalo MSA","clustalo",0);
	checkProg("MUSCLE5 MSA","MUSCLE5",0);

	
	
	
	
		#env MFFphylo
	checkProg("trimal","trimal",0);
	checkProg("mafft MSA","mafft",0);
	checkProg("raxmlng phylogeny","raxmlng",0);
	checkProg("fasttree phylogeny","fasttree",0);
	checkProg("iqtree phylogeny","iqtree",1);
	checkProg("ete3","ete3",0);
	checkProg("eggNOGmapper","emapper",1);
	
	checkProg("Kraken2 read classifier","kraken2",1);

	
	
	print "\n\n---------------------------------------------------\n";
	print "All essential tools seem to be present and working\n";
	
	if ($instSto eq ""){
		my $MF3Dir = getProgPaths("MFLRDir");
		$instSto ="$MF3Dir/helpers/install/progsChecked.sto";
	}
	system "touch $instSto";
	if ($doExit){
		print "Exiting MATAFILER\n";
		exit(0);
	}
}

sub checkMF{
	my $progSet = 1;
	$progSet = $_[0] if (@_ >= 1);
	#die "XX$progSet\n";
	my $expEnv = getProgPaths("CONDAbaseEnv");
	my $curEnv = `echo \$CONDA_DEFAULT_ENV`;
	chomp $curEnv;
	if ($curEnv ne $expEnv){
		print "Current Conda/Mamba environment ($curEnv) is not the expected one: $expEnv\nUse \"micromamba activate $expEnv\" and rerun MATAFILER\n";
		exit(23);
	}
	

	
	checkMFProgVers($progSet);

	my $MF3Dir = getProgPaths("MFLRDir");
	my $instSto ="$MF3Dir/helpers/install/progsChecked.sto";
	#print "$instSto\n";die;
	if (!-e $instSto){
		#die;
		print "Seems like MATAFILER did not yet check install paths, doing this now..\n";
		checkMFFInstall($instSto,0);
	}


	#die;
}




sub bam2cram($ $ $ $ $ $){#save further space: convert the bam to cram
	my ($iBAM,$REF,$del,$doCram, $stone, $numCore) = @_;
	my $ret = "";
	my $nxtCRAM = $iBAM;
	my $smtBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";

	if (!$doCram){return ($ret,"");}
	$nxtCRAM =~ s/\.bam$/\.cram/;
	if (!-e $iBAM && -e $nxtCRAM){print "CRAM exists, but no stone set\n"; system "rm -f $nxtCRAM";return $ret;}
	#my $stone = $iBAM;	$stone =~ s/\.bam$/\.cram\.sto/;
	$ret.="rm -f $nxtCRAM\n" if (-e $nxtCRAM);
	$ret.="$smtBin view -@ $numCore -T $REF -C -o $nxtCRAM $iBAM\n";
	$ret.="rm -f $iBAM\n" if ($del);
	$ret .= "touch $stone\n";
	#die $ret;
	return ($ret,$nxtCRAM);
}

sub cram2bsam{
	#doBAMSAM 1:bam 2:sam
	my ($iCRAM,$REF,$oBSAM,$doBAMSAM, $numCore) = @_;
	my $ret = "";
	my $smtBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
	if (!$doBAMSAM){return ($ret,"");}
	if ($doBAMSAM == 2){#SAM output
		$ret.="$smtBin view -h -@ $numCore -T $REF -o $oBSAM $iCRAM\n";
	} elsif ($doBAMSAM == 1){
		$ret.="$smtBin view -h -@ $numCore -T $REF -bo $oBSAM $iCRAM\n";
	}
	
	return ($ret);
}


sub getFileStr{
	my $inF = $_[0];
	my $req = 1;
	my $str = "";$req = $_[1] if (@_ > 1);
	my $tailN = -1;$tailN = $_[2] if (@_ > 2);
	if ($tailN >0){
		if (-e $inF){
			#print "Using tail: $inF\n";
			$str = `tail -n $tailN $inF`; chomp $str;
			return $str;
		}
	}
	my $fsMB = filsizeMB($inF);
	if ($fsMB > 5){
		print "getFileStr:: Very large file: ${fsMB}Mb $inF .. reading header only\n";
		$str = `head -n 2000 $inF`; chomp $str;
		return $str;
	}
	my ($I,$OK) = gzipopen($inF,"",$req,$req) ;
	die "getFileStr: Can't open $inF\n" if (!$OK && $req);
	return $str if (!$OK);
	while (my $l = <$I>){$str.=$l;}
	return $str;
}
sub displayPOTUS(){
	print "______ _____ _____ _   _ _____ \n| ___ \  _  |_   _| | | /  ___|\n| |_/ / | | | | | | | | \ `--. \n|  __/| | | | | | | | | |`--. \n";
	print "| |   \ \\_/ / | | | |_| /\\__/ /\n\\_|    \\___/  \\_/  \\___/\\____/ \n                               \n";
	#and ascii art
	print "II????????I??:...........~?I?IIIIIIIIIII\nIIII??I??I~,...,............:?IIIIIIIIII\nIIIIII??=,..,,.,..............=IIIIIIIII\nIIIII??:,,,,,..,.,..............+IIIIIII\n";
	print "IIIII,,,,,,,.,,,,.................IIIIII\nIII?..,,:======+=~::::::,,,,,.....,IIIII\nIIII:,,~=++++++++===~~~~~::::,,....?IIII\nIII.,,:=+++++++++++==~=~~~::::,,....IIII\n";
	print "II?,:~~=++?????+++++====~~::::,,,,,.IIII\nII=,~:==++????????++====~~::::,,,,..IIII\nII,,~:=+++????++++?+===~~~~:::::,,:.IIII\nII:,~+++++??????+++++====~~:::::,:,.IIII\n";
	print "II=,=+++++??????+?++?+==~~~~::,:::,,IIII\nIII,+++=+=~==~~==++==~:,......:,::..7III\n?~~:+++~,:~~,,::=+?=::~~~~~:,,.,::.,:,+I\n?+??=+=++~.=.,::=??=:,:,I..~..,:::::~:,I\n";
	print "++I?=++===?==~~==+?~:,,:~~~:,,::::~,~~,I\nI=?==++?+++??++??++~:::~~:~~~~::::~,:~7I\nI?==+=+?????????+++=,::~~~~==~::::~.,,77\nII+~+=+???+?????+++~::,,~~~~~~:::::::777\n";
	print "II?++=+????+?+=++?+=~::,,~~~~:::::~~7777\nIII+?=++???+=++?=+=:,,:::~=~::::,:,=7777\nIIIIII++++++=??++==~:::~~:~~::::,?777777\nIIIIII++++==+??++=~+~::::::::::,:7I77777\n";
	print "IIIIII?+++:+=++++=~~::,:,,:,~::,77777777\nIIIIIII+++~=.::::,:.,,,.,,:,~,,~77777777\nIIIIIIII=++=++++===~::::::,:.,:777777777\nIIIIIIIII===+?+++=~:~::::,,,,:~?77777777\n";
	print "IIIIIIIII+==++=+===~::::::,,:::I.7777777\nIIIIIIII?++=++++++==~~~:,,,,:,II..777777\nIIIIIII~.I++==++++==~~:,,,,:~III...?7777\nIII~...,.7=+++~~=~~~::,,,,,??I?,.....777\n";
	print "........?77?++++==~::::::??II?~.........\n......,.7777+++++~~=~:IIIIIII=..........\n........77777.++===~77777III+..........."
}

sub readTabbed($){
	my ($inF) = @_;
	my %ret;
	#open It,"<$inF" or die "Cant open tabbed infile $inF\n";
	my ($It,$OK) = gzipopen($inF,"Table file",1);

	while (<$It>){
		chomp;
		my @spl = split /\t/;
		$ret{$spl[0]} = $spl[1];
	}
	close $It;
	return \%ret;
}

sub readTable{
	my ($inF,$sep) = @_;
	my $merg = "";#can also merg into string..
	$merg = $_[2] if (@_>2);
	#skip lines?
	my $skip=0;
	$skip = $_[3] if (@_>3);#not used currently

	my %ret;
	my ($It,$OK) = gzipopen($inF,"Table file",1);
	#open It,"<$inF" or die "Cant open tabbed infile $inF\n";
	my $cnt=-1;
	while (my$l = <$It>){
		$cnt ++;
		next if ($cnt < $skip);
		#chomp;
		next if ($l =~ m/^#/);
		$l =~ s/\n$//;
		my @spl = split /$sep/,$l,-1;
		#if (@spl < $col){die"Requested column $col, but file $inF has only ".@spl." columns\n";}
		my $id = shift @spl;
		if ($merg eq ""){
			$ret{$id} = \@spl;
		} else {
			$ret{$id} = join($merg,@spl);
		}
	}
	close $It;
	return \%ret;
}



sub readTabbed3($ $){
	my ($inF,$col) = @_;
	my %ret;
	my ($It,$OK) = gzipopen($inF,"Tabbed file",1);
	#open It,"<$inF" or die "Cant open tabbed infile $inF\n";
	while (my$l = <$It>){
		#chomp;
		$l =~ s/\n$//;
		my @spl = split /\t/,$l,-1;
		if (@spl < $col){die"Requested column $col, but file $inF has only ".@spl." columns\n";}
		$ret{$spl[0]} = $spl[$col];
	}
	close $It;
	return \%ret;
}
sub readTabbed2{
	my $inF = $_[0];
	my $useLast = 0;
	$useLast = $_[1] if (@_>1);#not used currently
	
	my %ret;my $maxDepth=0;
	open It,"<$inF" or die "Cant open tabbed infile $inF\n";
	
	while (my $l=<It>){
		chomp $l;
		my @spl = split /\t/,$l,-1;
		my $ke;
		if ($useLast==1){
			$ke	= pop @spl;
		} else {
			$ke	= shift @spl;
		}
		#print "$ke\n";
		#print $ke." ";
		$ret{$ke} = \@spl;
		#die "@spl\n";
		if (scalar(@spl) > $maxDepth){$maxDepth=scalar(@spl);}
	}
	close It;
	return (\%ret,$maxDepth);
}


sub getSpecificDBpaths($ $){
	my ($curDB,$checkDBpreped) = @_;
	if ($curDB eq "mp3"){return "","",$curDB;}
	my $DBpath = "";	my $refDB = ""; my $shrtDB = "";
	#if ($curDB eq "NOG"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/eggNOG10/";	$refDB = "eggnog4.proteins.all.fa"; $shrtDB = $curDB;}
	#elsif ($curDB eq "MOH"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/MohFuncts/"; $refDB = "Extra_functions.fna";$shrtDB = $curDB;}
	#elsif ($curDB eq "CZy"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/CAZy/"; $refDB = "Cazys_2015.fasta";$shrtDB = $curDB;}
	#elsif ($curDB eq "ABR"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/ABR_FORS/"; $refDB = "ardb_and_reforghits.fa";$shrtDB = $curDB;}
	#elsif ($curDB eq "ABRc"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/ABR_Card/"; $refDB = "f11_and_card.faa";$shrtDB = $curDB; }
	#elsif ($curDB eq "KGE"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/KEGG/"; $refDB = "genus_eukaryotes.pep";$shrtDB = $curDB; }
	#elsif ($curDB eq "KGB"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/KEGG/"; $refDB = "species_prokaryotes.pep";$shrtDB = $curDB; }
	#elsif ($curDB eq "ACL"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/Aclame/"; $refDB = "aclame_proteins_all_0.4.fasta";$shrtDB = $curDB; }
	#elsif ($curDB eq "KGM"){$DBpath = "/g/bork3/home/hildebra/DB/FUNCT/KEGG/"; $refDB = "euk_pro.pep";$shrtDB = $curDB; }
	if ($curDB eq "NOG"){$DBpath = getProgPaths("eggNOG40_path_DB",0);	$refDB = "eggnog4.proteins.all.fa"; $shrtDB = $curDB;}
	elsif ($curDB eq "MOH"){$DBpath = getProgPaths("Moh_path_DB"); $refDB = "Extra_functions.faa";$shrtDB = $curDB;}
	elsif ($curDB eq "MOH2"){$DBpath = getProgPaths("Moh_path_DB"); $refDB = "Nitrogen_cycl_genes.faa";$shrtDB = $curDB;}
	elsif ($curDB eq "CZy"){$DBpath = getProgPaths("CAZy_path_DB"); $refDB = "Cazys_2019.fasta";$shrtDB = $curDB;}
	elsif ($curDB eq "ABR"){$DBpath = getProgPaths("ABRfors_path_DB"); $refDB = "ardb_and_reforghits.fa";$shrtDB = $curDB;}
	elsif ($curDB eq "ABRc"){$DBpath = getProgPaths("ABRcard_path_DB"); $refDB = "card.parsed.f11.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "KGE"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "genus_eukaryotes.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "KGB"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "species_prokaryotes.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "ACL"){$DBpath = getProgPaths("ACL_path_DB");$refDB = "aclame_proteins_all_0.4.fasta";$shrtDB = $curDB; }
	elsif ($curDB eq "KGM"){$DBpath = getProgPaths("KEGG_path_DB"); $refDB = "euk_pro.pep";$shrtDB = $curDB; }
	elsif ($curDB eq "TCDB"){$DBpath = getProgPaths("TCDB_path_DB"); $refDB = "tcdb.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "PTV"){$DBpath = getProgPaths("PATRIC_VIR_path_DB"); $refDB = "PATRIC_VF.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "PAB"){$DBpath = getProgPaths("ABprod_path_DB"); $refDB = "dedup_best_prod_predictions.faa";$shrtDB = $curDB; }
	elsif ($curDB eq "VDB"){$DBpath = getProgPaths("VirDB_path_DB"); $refDB = "VFDB_setB_pro.fas";$shrtDB = $curDB; }
	else {die"Unknown DB for func assignments: $curDB\n";}

	#basic file checks
	unless (-d $DBpath){die "getSpecificDBpaths:: Specified DB ($curDB) did not have valid DBpath: $DBpath\n";
	unless (-e "$DBpath/$refDB"){die "getSpecificDBpaths:: Specified DB ($curDB) did not have valid file: $DBpath/$refDB\n";
	
	if ($checkDBpreped){
		die "getSpecificDBpaths:: Can't find prepared diamond database at:\n$DBpath$refDB.db.dmnd" unless (-e "$DBpath$refDB.db.dmnd");
		die "getSpecificDBpaths:: Can't find length file at:\n$DBpath$refDB.length" unless (-e "$DBpath$refDB.length");
	}
	return ($DBpath ,$refDB ,$shrtDB );
}


sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
sub sortgzblast{ #function that checks if the diamond output was already sorted (required for paired end stuff with reads)
	my ($input) = @_;
	#print "$input\n";
	if ( $input =~ m/\.srt\.gz$/ ) { #redo srt in case there's a $trial file
		my $trial = $input; $trial =~ s/\.srt//; my $trialuse=0;
		if (!-e $input && !-e $trial){die "$input doesn't exist!\n";}
		if (-e $trial && -e $input && (-s $trial > -s $input)){$input = $trial; $trialuse=1; }#print "trial\n";
		if (-e $input && !$trialuse){return $input; }
		if (-e $trial && !-e $input){$input = $trial;}
		die "something went wrong in gzip sort 1\n$input\n" unless (-e $input);
	}
	my $tmpd=""; 
	my $cmd = "";
	my @chars = ("A".."Z", "a".."z");my $randstring;
	$randstring .= $chars[rand @chars] for 1..8;
	my $tmpDset=0;
	if (@_ >= 2){
		$tmpd = $_[1];$tmpDset=1;
	} 
	if ($tmpd eq ""){
		$input =~ m/^(.*\/)[^\/]+$/;$tmpd = $1;
	}
	my $input2=$input;
	$input2 =~ s/\.gz$//;
	my $input3=$input;
	$input3 =~ s/\.srt\.gz$//;
	if (!-e $input){ #maybe already something done here..
		if (!-e "$input2.srt.gz" && -e "$input.srt.gz"){system "mv $input.srt.gz $input2.srt.gz";}
		if (-e "$input2.srt.gz"){$input = "$input2.srt.gz";
		} elsif ( -e $input3 ){$input = $input3;
		}
	}
	#print $input."\n";
	unless ($input =~ m/\.srt\.gz$/){ #do sort (and maybe gz)
		if ($input =~ m/\.srt$/){
			$cmd = "gzip $input"; $input .= ".gz";
		} elsif ($input =~ m/\.gz$/) { #not sorted, but gz
			system "mkdir -p $tmpd" unless (-d "$tmpd");
			my $tmpf = "$tmpd/rawBLast$randstring.bla";
			$cmd = "zcat $input > $tmpf; sort $tmpf | gzip > $input2.srt.gz; rm -f $input $tmpf; ";
			if (!-e $input){die "Wrong file as input provided: $input\n";}
		} else { #not gz, not sort
			$cmd = "sort $input > $input.srt; gzip $input.srt; rm $input;";
		}
	}
	#die $cmd."\n$input\n";
	unless ($cmd eq ""){
		if (system $cmd) { die "$cmd \nfailed\n"; }
	}
	$input = "$input2.srt.gz";
	die "Something went wrong in sortgzblast 2\n" if (!-e $input);
	return $input;
}










