package Mods::FuncTools;
use warnings;
use Cwd 'abs_path';
#use File::chdir;

use strict;
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::TamocFunc qw( getSpecificDBpaths);
use Mods::GenoMetaAss qw( clenSplitFastas splitFastas systemW gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt );



use Exporter qw(import);
our @EXPORT_OK = qw(passBlast lambdaBl assignFuncPerGene calc_modules readGene2Func readGene2COG);


sub passBlast($ $ $){
	my $spl = shift;my $reqID = shift; my $inM8 = shift;
	my $lengthGood = 0;
	if (@{$spl} < 12 || !defined($spl->[13])|| !defined($spl->[12])){ #prob serious error.. delete file
		print "not enough entries:\n@{$spl}\n$inM8\n";
		return -2;
	}
	$lengthGood=1 if ($spl->[3] >= ($spl->[13]*0.1) && $spl->[3] >= ($spl->[12]*0.7)  ); #subject query
	if ($spl->[2] < ($reqID*0.985) || !$lengthGood)  {return 0;}
	return 1;
}
sub lambdaBl{
	my ($tar,$DB, $taxblastf,$BlastCores) = @_;
	my $exec = 1;
	$exec = $_[4] if (@_ > 4);

	my $cmd="";
	my $pigzBin = getProgPaths("pigz");
	my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";
	#my $lambdaIdxBin = $lambdaBin."_indexer";#getProgPaths("");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda_indexer";
#	if ($BlastCores > 20){$BlastCores = 20;}

	my $lamVer  = 0.4;
	
	my $lverTxt = `$lambdaBin --version`;
	if ( $lverTxt =~ m/.*lambda\d? version: (\d\.\d)[\.0-9]*.*/ ) {$lamVer = $1;}
	if ($lamVer < 3.0){die "Found Lambda v $lamVer at $lambdaBin\n no longer supports lambda below version 3.0\nPlease update\n";}
	
	
	if (  !-f "$DB.lba" && !-f "$DB.lba.gz" ) {
		print "\nBuilding LAMBDA index anew for $DB (this only happens the first time you use this ref DB, it may take several hours to build)..\n";
		my $cmdIdx =""; my $postCmd = "";
		my $DB1 = $DB; 
		if ($DB !~ m/\.fa$/){
			$DB =~ s/\.[^\.]+$/\.fa/;
			#this complicated stuff has to be done since lambda 3.0 doesn't accept databases not ending in .fa
			$cmdIdx .= "mv $DB1 $DB;\n";
			$postCmd .= "mv $DB $DB1;mv $DB.lba.gz $DB1.lba.gz;\n";
			#$postCmd .= "ln -s $DB1 $DB;mv $DB.lba.gz $DB1.lba.gz;\n";
		}
		$cmdIdx .=  "$lambdaBin mkindexn -t $BlastCores -d $DB ;\n";
		$cmdIdx .= "$pigzBin -p $BlastCores $DB.lba;\n";
		$cmdIdx .= $postCmd;
		$DB = $DB1;
		#die $cmdIdx."\n";
		if ($exec){
			if ( system($cmdIdx) ) {die( "Lamdba3 ref DB build failed\n$cmdIdx\n" );}
		} else {$cmd .= $cmdIdx ;}
	}

	#old lambda version.. no longer used..
	#if ( !-d "$DB.lambda" ) {
	#	print "Building LAMBDA index anew (may take up to an hour)..\n" if ($exec);
	#	my $cmdIdx = "";#"$lambdaIdxBin -p blastn -t ".int($BlastCores)." -d $DB";
		#$cmdIdx = "$lambdaIdxBin -p blastn -t $BlastCores -d $DB\n";
	#	$cmd .= $cmdIdx ;
		#if (system($cmdIdx)){die ("Lamdba ref DB build failed\n$cmdIdx\n");}
	#}
	
	my $outcols = "'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen'";
	#$cmd .= "$lambdaBin -t $BlastCores -id 30 -nm 30 -p blastn -e 1e-12 -so 7 -sl 16 -sd 1 -b 5 -pd on -q $tar -i $DB.lambda -o $taxblastf -oc $outcols;\n"; 
	$cmd .= "$lambdaBin searchn -t $BlastCores --percent-identity 30 --num-matches 30 --e-value 1e-12 -q $tar -i $DB.lba.gz -o $taxblastf --output-columns $outcols;";


	#die $cmd."\n";
	#$cmd .= "$lambdaBin -t $BlastCores -id 93  -nm 100 -p blastn -e 1e-40 -q $tar -oc \"std qlen slen\" -i $DB.lambda -o $taxblastf\n";
	#die $cmd."\n";
	systemW $cmd if (!-e $taxblastf && $exec);
	return $cmd;
}



sub readGene2Func{
	my $in = $_[0];
	my $cat = "";
	$cat = $_[1] if (@_ > 1);
	my $inF = "";
	#first try eggnogMapper
	
	
	if (-d $in && $cat ne ""){
		if ($cat =~ m/NOG|PFAM|EC|GO|KGM|BIGG/){
			$inF = "$in/Anno/Func/emapper/eggNOGmapper_". $cat .".geneAss";
		} else {
			$inF = "$in/Anno/Func/DIAass_$cat.srt.gzgeneAss.gz";
		}
	}
	print "Using reference annotation in $inF\n";
	#die "can't find input file $inF\n" unless (-e $inF);
	my %ret;
	my ($I,$OK) = gzipopen($inF,"Functional input file from gene catalog");
	while (<$I>){
		chomp;
		my @spl = split /\t/;
		$ret{$spl[0]} = $spl[1];
	}
	close $I;
	return \%ret;
}

sub readGene2COG{
	my ($inF) = @_;
	my %g2c; my %c2cat;
	open I,"<$inF" or die "Can't open $inF\n";
	print "reading gene -> COG file.. ";
	while (my $line=<I>){
		chomp $line;
		#NOG     COG5157 253     225     K       9796.ENSECAP00000013089,
		my ($cls,$COG,$p1,$p2,$cat,$genes) = split(/\t/,$line);
		$c2cat{$COG} = $cat;
		my @all = split(/,/,$genes);
		foreach (@all){
			#die "$_ $COG $cat\n";
			$g2c{$_} = $COG;
		}
	}
	print "Done\n";
	close I;
	return (\%g2c,\%c2cat);
}
sub min ($$) { $_[$_[0] > $_[1]] }
sub max ($$) { $_[$_[0] < $_[1]] }


sub buildFSdb{
	my ($DB, $DBout ,$ncore,$QSBoptHR,$qsubDir) = @_;
	my $prst5W = getProgPaths("PtostT5_Weights");
	#my $ar = splitFastas($query,$fastaSplits,$otpsHR->{splitPath}."/");
	#@subFls = @{$ar};
	#for (my $i =0 ; $i< @subFls;$i++){
	my $FSbin = getProgPaths("foldseek");
	my $DBcmd = "$FSbin createdb $DB $DBout --threads $ncore --prostt5-model $prst5W \n";
	my ($jN, $tmpCmd) = qsubSystem($qsubDir."foldseekDBprep.sh",$DBcmd,$ncore,"2G","diaDB","","",1,[],$QSBoptHR);
	return $jN;
}


#unifying script to diamond AA gene file against required DB & assign functions via the dia filter scripts of MATAFILER
#for now mostly diamond focused: NOG,KGM,KGE,KGB,ABRc,PTV,CZy,ACL,TCDB,PAB
#now-diamond: mp3 (pathogenic genes)
sub assignFuncPerGene{
	my ($query,$outD,$tmpD,$curDB) = @_;
	my $secCogBin = getProgPaths("secCogBin_scr");#MATAFILER script to filter raw mappings by diamond/blast etc
	my $diaBin = getProgPaths("diamond");#diamond binary
	my $prsMP3_scr = getProgPaths("mp3Prs");
	my $FSbin = getProgPaths("foldseek");
	my $avx2Constr =  getProgPaths("avx2_constraint",0);

	my $otpsHR = {};
	$otpsHR = $_[4] if (@_ > 4);
	my $doQsub=0; #if 0, then local execution
	my $QSBoptHR = {};
	my $qsubDir = $outD."qsubLOG/";;
	if (@_ > 5 && $_[5] != 0){
		$doQsub = 1;
		$QSBoptHR = $_[5] ;
		$qsubDir = $QSBoptHR->{qsubDir} if (exists($QSBoptHR->{qsubDir}) && $QSBoptHR->{qsubDir} ne "" );
	}
	my $exe = 0;
	if (@_ > 6 && $_[6] != 0){
		$exe = 1;
	}
	my $rstr = map { (q(a)..q(z))[rand(26)] } 1 .. 10;
	my $locTmpD = getProgPaths("nodeTmpDir")."/diaFunc/$rstr/";#specific cache for diamond

	my $localTmp = 0;#local tmp, requires different file copying
	#$localTmp = 1 if ($doQsub);
	my $fastaSplits = 1; #defaults
	my $ncore = 10;
	#my %opts = %{$otpsHR};
	$ncore = $otpsHR->{ncore} if (exists($otpsHR->{ncore}));
	$ncore = 1 if ($curDB eq "mp3"); #global override, since prog can't have more cores ..
	$fastaSplits = $otpsHR->{fastaSplits} if (exists($otpsHR->{fastaSplits}));
	$otpsHR->{splitPath} = $tmpD if (!exists($otpsHR->{splitPath}));
	$otpsHR->{redo} = 0 if (!exists($otpsHR->{redo}));
	$otpsHR->{keepSplits} = 0 if (!exists($otpsHR->{redo}));
	my $aligner = $otpsHR->{align};
	my $redo = $otpsHR->{redo};
	my $globalDiamondDependence = "";
	system "mkdir -p $outD" unless (-d $outD);
	if (!$otpsHR->{keepSplits}){
		clenSplitFastas($query,$otpsHR->{splitPath}."/");
	}
	
	
	#build DB
	my ($DBpath ,$refDB ,$shrtDB) = getSpecificDBpaths($curDB,0);
	if ($DBpath eq ""){print "Could not find DB $curDB !\n Exiting..\n"; exit(0);}
	#die "$DBpath\n";
	if ($exe && $DBpath ne "" ){
		my $DBcmd = "";
		unless (-e "$DBpath/$refDB.length"){
			my $genelengthScript = getProgPaths("genelength_scr");#= "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/geneLengthFasta.pl";
			$DBcmd .= "$genelengthScript $DBpath$refDB $DBpath$refDB.length\n";
		}
		if($aligner eq "diamond" ){
			$DBcmd .= "$diaBin makedb --in $DBpath$refDB -d $DBpath$refDB.db -p $ncore\n" ;
			if (!-e "$DBpath$refDB.db.dmnd" && $doQsub){
				my ($jN, $tmpCmd) = qsubSystem($qsubDir."DiamondDBprep.sh",$DBcmd,$ncore,"2G","diaDB","","",1,[],$QSBoptHR);
				$globalDiamondDependence = $jN;
			}
		} elsif ($aligner eq "foldseek" ){ 
			if (!-e "$DBpath$refDB.DB3di"){
				my $jns = buildFSdb("$DBpath$refDB", "$DBpath$refDB.DB3di",$ncore,$QSBoptHR,$qsubDir);
				$globalDiamondDependence = $jns;
			}
		} else {die"FuncTools.pm::assignFuncPerGene: Unknown aligner: $aligner\n";}

		if (!$doQsub){
			systemW $DBcmd;
		}
		#die "$DBcmd\n";
	}
	
	my @subFls ;
#	print "Splitting catalog into $fastaSplits\n" if ($fastaSplits>1);
	my @jdeps; my @allFiles;
	my $allAss = "$outD/DIAass_$shrtDB.srt.gz";
	my $tarAnno = "${allAss}geneAss.gz";
	system "rm -f $allAss" if ($redo);
	system "rm -f $tarAnno" if ($redo);
	my $calcDia = 1;$calcDia = 0 if (-e $allAss);
	my $interpDia = 1;$interpDia = 0 if (-e $tarAnno);
	#my $N = 20;
	my $jdep=""; my $qCmd = "";
	
	#die "!$calcDia && !$interpDia $curDB\n$allAss\n$tarAnno\n";
	return $allAss,$jdep if (!$calcDia && !$interpDia);

	
	#die "$calcDia\n$allAss\n";
	#alignment options (and defaults)
	$otpsHR->{eval} = 1e-7 if (!exists($otpsHR->{eval}));
	$otpsHR->{percID} = 25 if (!exists($otpsHR->{percID}));
	$otpsHR->{minAlignLen} = 60 if (!exists($otpsHR->{minAlignLen}));
	$otpsHR->{minBitScore} = 60 if (!exists($otpsHR->{minBitScore}));
	$otpsHR->{minPercSbjCov} = 0.3 if (!exists($otpsHR->{minPercSbjCov}));
	$otpsHR->{bacNOG} = 0 if (!exists($otpsHR->{bacNOG}));

	print "$query assigned to $curDB ($fastaSplits splits, $ncore cores)\n" if ($calcDia || $interpDia);
	my $mem = 20;
	$mem = 160 if ($shrtDB eq "NOG" || $shrtDB eq "KGM");

	my $tmpD2 = "$tmpD/$curDB/";
	if ($calcDia){
		print "Diamond pars: eval=$otpsHR->{eval}, percID=$otpsHR->{percID}, minAlLength=$otpsHR->{minAlignLen}, minBitScore=$otpsHR->{minBitScore}, minPercSbjCov=$otpsHR->{minPercSbjCov}}\n" if ($calcDia);
		my $ar = splitFastas($query,$fastaSplits,$otpsHR->{splitPath}."/");
		@subFls = @{$ar};
		for (my $i =0 ; $i< @subFls;$i++){
			my $tmpD3 = "$tmpD2/R$i/";
			my $cmd = "mkdir -p $tmpD3\nmkdir -p $locTmpD\n";
			my $outF = "";
			if ($curDB eq "mp3"){
				my $mp3Dir = getProgPaths("mp3");
				#$CWD = $mp3Dir;
				chdir $mp3Dir;	
				$cmd .= "$mp3Dir/./mp3 $subFls[$i] 1 50 0.2\n";#1=complete genes,2=metag prots
				$subFls[$i] =~ m/\/([^\/]+$)/; my $fnm=$1;
				if ($localTmp){
					$outF = "$outD/$fnm.Hybrid.result.gz"; #this is sometimes shared, sometimes local tmp
					$cmd .= "$prsMP3_scr $subFls[$i].Hybrid.result $outF\n";
					
				} else {
					$outF = "$subFls[$i].Hybrid.result.gz";
					$cmd .= "$prsMP3_scr $subFls[$i].Hybrid.result $outF\n";
				}
				$cmd .= "rm $subFls[$i].Hybrid.result\n";
			}elsif ($aligner eq "foldseek"){
				my $prst5W = getProgPaths("PtostT5_Weights");
				$cmd .= "foldseek easy-search $subFls[$i] tinycazy_db result.m8 $tmpD3 --prostt5-model $prst5W";
			} else {
				#my $outF = "$GCd/DiaAssignment.sub.$i";
				if ($localTmp){
					$outF = "$outD/DiaAs.sub.$i.$shrtDB.gz";
				} else {
					$outF = "$tmpD3/DiaAs.sub.$i.$shrtDB.gz";
				}
				$cmd .= "$diaBin blastp -f tab --compress 1 --quiet -t $tmpD3 -d $DBpath$refDB.db -q $subFls[$i] -e $otpsHR->{eval} -o $outF -p $ncore\n";#--sensitive
				#--memory-limit ". int($mem *0.8-0.8) ."
				#$cmd = "$diaBin blastp -f tab --compress 1 --sensitive --quiet -d $eggDB.db -q $subFls[$i] -k 3 -e 0.001 -o $outF -p $ncore\n";
				#$cmd .= "$diaBin view -a $outF.tmp -o $outF -f tab\nrm $outF.tmp* $subFls[$i] \n";
				$cmd .= "rm -r $tmpD3\n" if ($localTmp);
			}
			#die "$cmd\n";
			system "rm -f $outF" if ($redo);
			if ($calcDia && !-e $outF && $exe){
						#die "$cmd\n";
				if ($doQsub){
					my @preCons = @{$QSBoptHR->{constraint}};
					push(@{$QSBoptHR->{constraint}}, $avx2Constr);#--constraint=sse4
					my $preHDDspace=$QSBoptHR->{tmpSpace};
					$QSBoptHR->{tmpSpace} = ($mem*2) . "G";
					my ($jobName,$mptCmd) = qsubSystem($qsubDir."D$shrtDB.$i.sh",$cmd,$ncore,($mem/$ncore)."G","D$shrtDB$i",$globalDiamondDependence,"",1,[],$QSBoptHR); #$jdep.";".
					@{$QSBoptHR->{constraint}} = @preCons;
					$QSBoptHR->{tmpSpace} = $preHDDspace;
					push(@jdeps,$jobName);
					#die "$jobName\n";
				} else {
					systemW $cmd;
				}
				#print $qsubDir."Diamond.sh\n";
			}
			push(@allFiles,$outF);
			#if ($i==5){die;}
		}
	}
	#last job that converges all
	my $cmd = "";
	my $catcmd = "cat";
	if (!$calcDia){
	
	} elsif (@allFiles == 1){
		if ($allFiles[0] =~ m/\.gz$/){
			$cmd .= "mv $allFiles[0] $allAss\n";
		} else {
			$cmd .= "gzip -c $allFiles[0] > $allAss\nrm $allFiles[0]\n";
		}
	} else {
		if ($allAss !~ m/\.gz$/ && $allFiles[0] =~ m/\.gz$/){
			$allAss.=".gz" ;
		} elsif($allAss =~ m/\.gz$/ && $allFiles[0] !~ m/\.gz$/){ #zip output
			$catcmd = "gzip -c";
		}
		#my $cmd= "cat ".join(" ",@allFiles). " > $allAss\n";   #
		$cmd .= "$catcmd ".join(" ",@allFiles). " > $allAss\n";
		$cmd .= "rm -f ".join(" ",@allFiles) . "\n";
	}
	if ($curDB eq "mp3"){
		$cmd .= "mv $allAss ${allAss}geneAss.gz\ntouch $allAss\n";
	} else {
		$cmd .= "$secCogBin -i $allAss -DB $shrtDB -singleSpecies 1  -bacNOG $otpsHR->{bacNOG} -KOfromNOG 0 -eggNOGmap 1 -calcGeneLengthNorm 0 -lenientCardAssignments 2 ";
		$cmd .= "-mode 2 -CPU $ncore -percID $otpsHR->{percID} -LF $DBpath/$refDB.length -DButil $DBpath -tmp $tmpD2 -eggNOGmap 0 -minPercSbjCov $otpsHR->{minPercSbjCov} ";
		$cmd .= "-minBitScore $otpsHR->{minBitScore} -minAlignLen $otpsHR->{minAlignLen} -eval $otpsHR->{eval}\n";
	}

	#die "$cmd";
	$cmd .= "rm -f ".join(" ",@subFls)."\n" unless($otpsHR->{keepSplits});
	if ($exe){
		if ($interpDia && $doQsub){
			($jdep,$qCmd) = qsubSystem($qsubDir."colDIA$shrtDB.sh",$cmd,1,"30G","ColDIA",join(";",@jdeps),"",1,[],$QSBoptHR);
		} elsif ($interpDia) {
			systemW $cmd;
		}
	}
	return $allAss,$jdep;
}



sub calc_modules{
	my ($inMat,$outD,$ModCompl,$EnzCompl,$exeLoc) = @_;
	my $rareBin = getProgPaths("rare");
	my $modD = getProgPaths("Module_path_DB");
	#my $modD = "/g/bork3/home/hildebra/DB/FUNCT/myModules/Feb16/";
	my @modDBfs = ("module_new.list","module_c.list","modg.list","module_s.list","modGBM.list");
	my @modDBodir = ("modules","metaCyc","BSB","SEED","GBM");
	my @modDescr = ("mod.descr","modc.descr","modg.descr","mods.descr","modGBM.descr");
	my @modHiera = ("mod_hiera.txt","modc_hiera.txt","modg_hiera.txt","mods_hiera.txt","modGBM_hiera.txt");
	my @modShort = ("modKEEG","modMC","modGMM","modSEED","modGBM");
	my $cmdMod = "";
	for (my $k=0;$k<@modDBfs;$k++){
		next if ($k==1); #skip metacyc due to error
		my $keggDB = $modD.$modDBfs[$k];
		my $outD2 = "$outD/$modDBodir[$k]/";
		system "mkdir -p $outD2" unless (-d $outD2);
		my $outMat = $outD2.$modShort[$k];
		next if (-e $outD2."$modShort[$k].mat");
		#die "$inMat\n";
		$cmdMod .= "$rareBin module -i $inMat -o $outMat -refMods $keggDB -description $modD/$modDescr[$k] -hiera  $modD/$modHiera[$k] -redundancy 5 -writeExtraModEstimates -moduleCompl $ModCompl -enzymeCompl $EnzCompl -collapseDblModules\n";
		if ($exeLoc){
			die "Failed module calc:\n$cmdMod\n" if (system "$cmdMod");
		}
	}
	return $cmdMod;
}
