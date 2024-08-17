#!/usr/bin/env perl
#makes an abundance matrix of NOGs across samples
#./combine_DIA.pl /g/scb/bork/hildebra/SNP/GNMass2_singl/
#./combine_DIA.pl /g/scb/bork/hildebra/Tamoc/FinSoil/
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipopen readMapS);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::TamocFunc qw (getSpecificDBpaths );

sub writeMat;

#CZyparse.CZy.cat.cnts.gz
#CZyparse.CZy.gene.cnts.gz
#MOHparse.MOH.cat.cnts.gz
#MOHparse.MOH.gene.cnts.gz
#NOGparse.NOG.cat.cnts.gz
#NOGparse.NOG.gene.cnts.gz
my @DBs = ("NOG"); # ("ABRc","KGB","KGE","CZy","NOG","ABR","MOH","KGM","ACL","PAB");
my @modDBs = ("KGM");
#@DBs = ("NOG","KGB","KGE"); # ("KGB","KGE","CZy","NOG","ABR","MOH");#("NOG","CZy","MOH");
#@DBs = ("NOG");#("NOG","CZy","MOH");
#my $rareBin = "/g/bork3/home/hildebra/dev/C++/rare/./rare_fix"; #TODO temp
my $rareBin = getProgPaths("rare");#"/g/bork3/home/hildebra/dev/C++/Rarefaction/rtk/rtk";

#my $modD = "/g/bork3/home/hildebra/DB/FUNCT/myModules/Feb16/";
my $modD = getProgPaths("Module_path_DB");
my @modDBfs = ("module_new.list","module_c.list","modg.list","module_s.list");
#my @modDBfs = ("modg.list");
my @modDBodir = ("modules","metaCyc","BSB","SEED");
my @modDescr = ("mod.descr","modc.descr","modg.descr","mods.descr");
my @modHiera = ("mod_hiera.txt","modc_hiera.txt","modg_hiera.txt","mods_hiera.txt");


my $inD = $ARGV[0];
my $DBstr = "";
$DBstr = $ARGV[1] if (@ARGV>1);
chomp $DBstr;
@DBs = split/,/,$DBstr if ($DBstr ne "");
my $dieOnMissing = 0;

my $calcModules=0;

$inD .= "/" unless ($inD =~ m/\/$/);
my $outD1 = $inD."pseudoGC/FUNCT/";
system "mkdir -p $outD1" unless (-d $outD1);
my ($hrm,$hr2) = readMapS($inD."LOGandSUB/inmap.txt",-1);
my %map = %{$hrm};
my @samples = @{$map{opt}{smpl_order}};
#if (!-d $map{$samples[0]}{wrdir}){ #check for new way of reading dir
#	($hrm,$hr2) = readMapS($inD."LOGandSUB/inmap.txt",0);
#	%map = %{$hrm};
#	@samples = @{$map{smpl_order}};
#}
#die "@samples\n";
#find out if subfolders exists
my %DBsD;my $testD;
my %Taxs; #contains the subparts that the matrix is split into
#look what dirs do exists and if they have the DB signature files in them
#die "$map{$samples[0]}{dir}\n";

print "Working through ".@samples." samples\n";

foreach my $DB (@DBs){
	my @subDirs = ();
	$testD = $map{$samples[0]}{wrdir}."/diamond/";
	print $testD."\n";
	opendir D,$testD or die "Can't open dir $testD\n";
	my $cnt=0;
	while (my $dd = readdir D){
		if (-d $testD.$dd  && $dd !~ m/^\./){
			my (@FList) = glob($testD.$dd."/$DB*");
			#print "@FList\n";
			next unless (@FList > 0);
			#print $dd . " ".$testD.$dd."/$DB*\n";
			push @subDirs, $dd ;
		}
		$cnt++;
		if ($cnt > 1000){die "Abort while1\n";}
		#print $dd if (-d $testD.$dd);
	}
	closedir D;
	#die "@subDirs\n\n";
	$DBsD{$DB} = \@subDirs;

	foreach my $smpl(@samples){ #check that diamond blastx finished properly
		my $dir2rd = $map{$smpl}{wrdir};
		if (!-e "$dir2rd/diamond/dia.$DB.blast.srt.gz"){
			print "Can't find base diamond output for $dir2rd\n";
		}
	}

}
#die;
#print "D1\n";
#look in each subdir specific to DB for valid *cat file to get the TAXs available
foreach my $DB (@DBs){
	if ($DB eq "KGB"){$Taxs{$DB} = ["BAC"]; $calcModules=1; next;}
	if ($DB eq "KGE"){$Taxs{$DB} = ["EUK"]; $calcModules = 1; next;}
	if ($DB eq "KGM"){$Taxs{$DB} = ["EUK","BAC","UNC"]; $calcModules = 1; next;}
	my @tmp2 = @{$DBsD{$DB}};# die "@tmp2\n";
	print "@{$DBsD{$DB}}\n".@tmp2."\n";
	$testD = $map{$samples[0]}{wrdir}."/diamond/".${$DBsD{$DB}}[0]."/";
	#die $testD;
	opendir D,$testD or die "Can't open dir $testD\n";
	my @files = readdir D;closedir D;
	my @tmp; foreach (@files){push @tmp,$_ if (m/^$DB.*cat.*/);}
	if (@tmp == 0){#search for another pattern
		#print "XXXX\n";
		foreach (@files){push @tmp,$_ if (m/^${DB}parse*/);}
	}
	my @ins;
	foreach my $v (@tmp){
		push @ins,$1 if ($v =~ m/$DB.?parse\.$DB\.(\w+)\.cnt/);
	}
	print "@ins\n";
	$Taxs{$DB} = \@ins;
}
#die;
#print "D2\n";

#die;
#push @subDirs,"" if (@subDirs == 0);
#go over different DBs
foreach my $DB (@DBs){
	#last; #tmp deactivate
	print "---------------  $DB --------------\n";
	my $outD = $outD1."$DB/";
	system "mkdir -p $outD" unless (-d $outD);
	#sub up the contetns of each subdir
	foreach my $deepDir (@{$DBsD{$DB}}){
		my $deepDirT = $deepDir; $deepDirT =~ s/CNT_//;
		print $deepDirT."\n";
		foreach my $NORM ("cnt","GLN"){
			my %tCOGs=(); my %tCATs=();
			my %tCGsums=(); my %tCTsums = ();
			my %failC; my %failCOG;
			my $didreadCOG=0;
			#die $NORM."\n";
			foreach my $TAX (@{$Taxs{$DB}}){
				my %COGs=(); my %CATs=();
				my %CGsums=(); my %CTsums = ();
				#$genes{1}{gg} = "falk";
				#die $genes{1}{gg};
				#die "$TAX\n";
				foreach my $smpl(@samples){
					my $dir2rd = $map{$smpl}{wrdir};
					my $SmplName = $map{$smpl}{SmplID};
					#print $SmplName."\n";
					my $DiaCOGf = "$dir2rd/diamond/$deepDir/$DB"."parse.$DB.$TAX.$NORM.cat.cnts";
					my $DiaCATf = "";# = "$dir2rd/diamond/$deepDir/$DB"."parse.$TAX.$NORM.CATcnts";
					
					#print $DiaCATf."\n";
					#MOHparse.MOH.gene.cnts.gz
					if ($DB eq "CZy"){
					}elsif ($DB ne "NOG" &&$DB ne "PTV" && $DB ne "PAB" &&$DB ne "PTV" &&  $DB ne "KGB" && $DB ne "KGE" && $DB ne "KGM"&& $DB ne "ABRc"){
						$DiaCATf = "$dir2rd/diamond/$deepDir/$DB"."parse.$DB.$TAX.$NORM.gene.cnts.gz" ;
					
					} else {
						$DiaCATf = "$dir2rd/diamond/$deepDir/$DB"."parse.$TAX.$NORM.CATcnts";
					}
					if ( $DB ne "KGB" && $DB ne "KGE" && $DB ne "KGM"){
						#print "$DiaCOGf\n";
						my ($I,$ok) = gzipopen($DiaCOGf,"Diamond $DB out",0);
						$didreadCOG=1;
						if ($ok == 0){ 
							$failC{$smpl}=1;
							#remove stone so TAMOC can redo..
							print "$smpl fail\n" if ($dieOnMissing);
							system "rm -f $dir2rd/diamond/dia.$DB.blast.gz.stone";
						} else {
							while (my $l = <$I>){
								next if (length($l) < 4);
								chomp $l;
								my @spl = split /\t/,$l;
								#print "@spl\n";#		die $l;
								$COGs{$spl[0]}{$smpl} = $spl[1]; 
								if (exists $tCOGs{$spl[0]}{$smpl}){$tCOGs{$spl[0]}{$smpl} += $spl[1];} else {$tCOGs{$spl[0]}{$smpl} = $spl[1]; }
								if (exists($CGsums{$spl[0]})){	$CGsums{$spl[0]} += $spl[1];} else {$CGsums{$spl[0]} = $spl[1];	}
								if (exists($tCGsums{$spl[0]})){	$tCGsums{$spl[0]} += $spl[1];} else {$tCGsums{$spl[0]} = $spl[1];	}
							}
							close $I;
							#die;
						}
					}
					#die "$DiaCATf\n";
					next if ($DiaCATf eq "");
					my ($I,$ok) = gzipopen($DiaCATf,"Diamond cat",0);
					if (!$ok){ $failCOG{$smpl}=1;
						print "$smpl fail\n";
						system "rm -f $dir2rd/diamond/dia.$DB.blast.gz.stone";
					} else {
						while (my $l = <$I>){
							chomp $l;
							next if (length($l) < 4);
							my @spl = split /\t/,$l;
							#die "$l\n" unless ($spl[1] =~ m/^[\d\.]+$/);
							#print "@spl\n";		die $l;
							$CATs{$spl[0]}{$smpl} = $spl[1]; 
							if (exists $tCATs{$spl[0]}{$smpl}){$tCATs{$spl[0]}{$smpl} += $spl[1]; } else {$tCATs{$spl[0]}{$smpl} = $spl[1]; }
							if (exists($CTsums{$spl[0]})){$CTsums{$spl[0]} += $spl[1];} else {$CTsums{$spl[0]} = $spl[1];}
							if (exists($tCTsums{$spl[0]})){$tCTsums{$spl[0]} += $spl[1];} else {$tCTsums{$spl[0]} = $spl[1];}
						}
						close $I;
					}
				}
				
				#write category matrix
				writeMat("$outD/","$DB".".CAT.mat.$TAX.$NORM.$deepDirT.txt",\%CTsums,\%CATs,\%failC);
				#write COG matrix
				writeMat("$outD/","$DB".".mat.$TAX.$NORM.$deepDirT.txt",\%CGsums,\%COGs,\%failCOG) if ($didreadCOG);
			}#loop over DB versions
			if (@{$Taxs{$DB}} > 1){
				#write total category matrix
				writeMat("$outD/","$DB".".CAT.mat.ALL.$NORM.$deepDirT.txt",\%tCTsums,\%tCATs,\%failC);
				#write total COG matrix
				writeMat("$outD/","$DB".".mat.ALL.$NORM.$deepDirT.txt",\%tCGsums,\%tCOGs,\%failCOG)if ($didreadCOG);
			}
		} #loop over TAX
	}#loop over norm
} #loop around sub/deepDirs 
#print "$DBsD{KGE}\n$DBsD{KGB}\n";
#now do module abundance estimates
my ($ModCompl,$EnzCompl) = (0.5,0.5);
foreach my $DB (@modDBs){
	unless ($DB eq "KGB" || $DB eq "KGM" || $DB eq "KGE"){next;}
	for (my $k=0;$k<@modDBfs;$k++){
		foreach my $deepDir (@{$DBsD{$DB}}){
	#print "X\n";
			#print "DDD\n";
			my $inD = $outD1."$DB/";
			foreach my $etag (@{$Taxs{$DB}}){#= "BAC";
				my $eval = $deepDir; $eval =~ s/^CNT_//;
				$etag = "EUK" if ($DB eq "KGE");
				my $keggDB = $modD.$modDBfs[$k];
				
				my $outD = "$inD/$modDBodir[$k]/";
				system "mkdir -p $outD" unless (-d $outD);
				my $outMat = $outD."KEGG$eval.$etag.GLN.mod";
				my $inMat = "$inD/$DB.CAT.mat.$etag.GLN.$eval.txt";
				my $cmdMod = "$rareBin module -i $inMat -o $outMat -description $modD/$modDescr[$k] -hiera  $modD/$modHiera[$k] -refMods $keggDB -redundancy 5 -moduleCompl $ModCompl -enzymeCompl $EnzCompl -collapseDblModules\n";
				$outMat = $outD."KEGG$eval.$etag.cnt.mod";
				$inMat = "$inD/$DB.CAT.mat.$etag.cnt.$eval.txt";
				$cmdMod .= "$rareBin module -i $inMat -o $outMat -refMods $keggDB -description $modD/$modDescr[$k] -hiera  $modD/$modHiera[$k] -redundancy 5 -moduleCompl $ModCompl -enzymeCompl $EnzCompl -collapseDblModules\n";
				#print "$cmdMod\n";
				system $cmdMod."\n";
				#die $cmdMod."\n";
			}
			#die;
		}
	}
}

print "All done\n";
exit(0);

sub writeMat(){
	my ($odir,$oname,$hr1,$hr2,$hr3) = @_;
	my %CTsums = %{$hr1}; my %CATs = %{$hr2};
	my %fails = %{$hr3};
	open O,">$odir/$oname" or die "Can't open matrix file $odir/$oname\n";
	print O "$oname";
	foreach my $smpl (@samples){
		print O "\t$smpl";
	}
	my @srtCTsum = (sort { $CTsums{$b} <=> $CTsums{$a} } keys %CTsums); 
	foreach my $cat (@srtCTsum){
		print O "\n$cat";
		foreach my $smpl (@samples){
			if (exists $fails{$smpl}){ #sample output file not found..
				print O "\t-";
			} elsif (exists($CATs{$cat}{$smpl})){
				print O "\t".$CATs{$cat}{$smpl};
			} else {
				print O "\t0";
			}
		}
	}
	close O;
}
