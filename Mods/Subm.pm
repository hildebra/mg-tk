package Mods::Subm;

use warnings;
use strict;
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths convert2Gb);


use Exporter qw(import);
our @EXPORT_OK = qw( findQsubSys emptyQsubOpt qsubSystem qsubSystem2 qsubSystemJobAlive
		qsubSystemWaitMaxJobs MFnext add2SampleDeps numUserJobs
		
		);#





sub randStr($){ #will be prefixed to jobname, to make jobs unique to each MF run
	my ($len) = @_;
	my @letters=('A'..'Z','a'..'z',1..9);
	my @letters2=('A'..'Z','a'..'z');
	my $total=scalar(@letters);
	my $newletter ="";
	$newletter = $letters2[rand scalar(@letters)];
	for (my $i=1;$i<$len;$i++){
		$newletter .= $letters[rand $total];
	}
	return $newletter;
}


sub numPendingJobs($){
	my ($optHR) = @_;
	my $qmode = "slurm"; $qmode = $optHR->{qmode} if (defined($optHR->{qmode}));
	my $srchCmd="" ;#= "squeue -u \$USER  -t PENDING | wc -l";
	my $num = 0;
	if ($qmode eq "slurm"){
		$srchCmd = "squeue -u \$USER -t PENDING | wc -l";
	} elsif ($qmode eq "sge"){
		$srchCmd = "qstat | grep \$USER  | wc -l";
		die "Subm.pm::numPendingJobs() not implemented for sge!\n";
	} elsif ($qmode eq "bash"){
		return 0;
	} else {$srchCmd="bsub  | wc -l";
		die "Subm.pm::numPendingJobs() not implemented for bsub!\n";
	}
	$num = `$srchCmd`; chomp $num; $num -=1;
	return $num;
}
sub numUserJobs{
	my ($optHR) = $_[0];
	my $rmSelf=0; $rmSelf = $_[1] if (@_>1);
	my $qmode = "slurm"; $qmode = $optHR->{qmode} if (defined($optHR->{qmode}));
	my $srchCmd ="";#= "squeue -u \$USER   | wc -l";
	if ($qmode eq "slurm"){
		$srchCmd = "squeue -u \$USER  | wc -l";
	} elsif ($qmode eq "sge"){
		$srchCmd = "qstat | grep \$USER  | wc -l";
	} elsif ($qmode eq "bash"){
		return 0;
	} else {$srchCmd="bsub  | wc -l";
	}
	my $num = 0;
	$num = `$srchCmd`; chomp $num;$num -=1;


	if ($rmSelf && $qmode eq "slurm"){
		my $SjobID = `echo \$SLURM_JOBID`; chomp $SjobID;
		#print "\"$SjobID\"\n";
		if ($SjobID ne ""){$num --;}
	}

	return $num;
}


sub findQsubSys($){
	my $iniVal = "";
	$iniVal = $_[0] if (@_ > 0);
	#my $iniVal = "lsf";
	if ($iniVal ne ""){
		$iniVal = lc $iniVal; 
		$iniVal = "lsf" if ($iniVal eq "bsub");
		$iniVal = "sge" if ($iniVal eq "qsub");
		$iniVal = "slurm" if ($iniVal eq "sbatch");
	} else {
		$iniVal = "lsf";
		my $bpath = `which bsub  2>/dev/null`;chomp $bpath;my $bpresent=0; 
		$bpresent=1 if ($bpath !~ m/\n/ && -e $bpath);
		my $spath = `which sbatch  2>/dev/null`;chomp $spath;my $spresent=0; 
		$spresent=1 if ($spath !~ m/\n/ && -e $spath);
		my $qpath = `which qsub  2>/dev/null`; chomp $qpath;
		my $qpresent=0; $qpresent=1 if ($qpath !~ m/\n/ && -e $qpath);
		#print "$qpath\n";
		if ($spresent ){#slurm gets preference
			$iniVal="slurm";
		}elsif (!$bpresent && $qpresent){
			$iniVal = "sge";
		}elsif (!$qpresent && !$bpresent && !$spresent){
			print "Warning: No queing system found (sbatch / qsub / bsub command)\nUsing LSF (bsub), though this will likely cause errors\n";
			
		}
	print "Using qsubsystem: $iniVal\n";
	}
	#die;
	return $iniVal;
}
sub emptyQsubOpt{
	my ($doSubm) = $_[0];
	my $locChkStr = $_[1];
	my $qmode = "";
	$qmode = $_[2] if (@_ > 2);
	
	if (@_ > 2){$qmode = $_[2];}
	$qmode = findQsubSys($qmode);
	die "qsub system mode has to be \'lsf\', \'bash\', \'slurm\' or \'sge\'!\n" if ($qmode ne "lsf" &&$qmode ne "slurm" && $qmode ne "sge"&& $qmode ne "bash");
	my $MFdir = getProgPaths("MFLRDir");
	my $longQ = getProgPaths("longQueue",0); my $shortQ =  getProgPaths("shortQueue",0); my $medQ = getProgPaths("mediumQueue",1);
	#die "$shortQ\n";
	my $gpuQ = getProgPaths("gpuQueue",0);
	my $himemQ = getProgPaths("highMemQueue",0);
	if ($longQ eq ""){$longQ =  $medQ;}
	if ($medQ eq "" ){die "FATAL: no medium queue defined!\n";};
	if ($gpuQ eq "" ){$gpuQ = $medQ;};
	if ($himemQ eq "" ){$himemQ = $medQ;};
	if ($shortQ eq "" ){$shortQ = $medQ;};
	my $xtraNodeCmds = getProgPaths("subXtraCmd",0);
	$xtraNodeCmds = "" unless (defined $xtraNodeCmds);
	my $medTime = getProgPaths("medTime",0);	my $shortTime = getProgPaths("shortTime",0);
	my $longTime = getProgPaths("longTime",0);
	my $subConfig = getProgPaths("submissionConfig",0);
	my @constr = ();
	if ($subConfig =~ s/--constraint=(\S+)//){
		#print "!!! $1\n";
		push(@constr, $1);
	}
	chomp($subConfig);
	@constr = grep(/\S/, @constr);
	#die "@constr\n$subConfig\nYW\n";
	
	#if ($qmode eq "slurm"){$shortQ = "htc"; $longQ="htc";}#$shortQ = "1day"; $longQ="1month";}
	my %ret = (
		rTag => randStr(3),
		doSubmit => $doSubm,
		LocationCheckStrg => $locChkStr,
		doSync => 0,
		longQueue => $longQ,
		gpuQueue => $gpuQ,
		highMemQueue => $himemQ,
		longTime => $longTime,#7days
		medQueue => $medQ,
		medTime => $medTime,#"24:00:00",
		shortQueue => $shortQ,
		shortTime => $shortTime, #2hrs
		useLongQueue => 0,
		useGPUQueue => 0,
		useShortQueue => 0,
		useHiMemQueue => 0,
		submissionConfig => $subConfig,
		constraint => \@constr,
		qsubPEenv => getProgPaths("qsubPEenv"),
		perl5lib => "$MFdir:\$PERL5LIB",
		cpplib => "",
		tmpSpace => 15, #default was 15G; unit is G
		tmpSpaceTag => getProgPaths("nodeTmpDirTAG",0),
		LOCKfile => "",
		tmpMinG => 10,
		afterAny => 0,
		excludeNodes => "",
		xtraNodeCmds => $xtraNodeCmds,
		qmode => $qmode,
		wcKeysForJob => "",
		#LOG => undef,
	);
	#die "$MFdir\n";
	return \%ret;
}

sub qsubSystemJobAlive{
	my ($jAr,$optHR) = @_;
	my $killFailedJobs=0;
	$killFailedJobs = $_[2] if (@_ > 2);
	my @jobs = @{$jAr};
	#clean up @jobs
	my %jobsCl;  foreach my $kj (@jobs){my @spl = split /;/,$kj; for (@spl){next if ($_ eq "");$jobsCl{$_}++;}}
	@jobs = keys %jobsCl;
	
	
	my $qmode = $optHR->{qmode};
	my $cmd1="";
	my $rTag = $optHR->{rTag};

	if ($qmode eq "slurm"){
		$cmd1 = "squeue -u \$USER ";
		for (@jobs) {s/$rTag//;}
	} elsif ($qmode eq "sge"){
		$cmd1 = "qstat | grep \$USER "
	} elsif ($qmode eq "bash"){
		return;
	} else {$cmd1="bsub";
	}
	my $cmd = "$cmd1";# | grep $_ | wc -l";
	#my $num = `$cmd`; chomp $num;
	my $num  = `$cmd`;
	my $jobsCheckd= scalar @jobs;
	#print "XX\n@jobs\n\n";
	foreach (@jobs){
		$jobsCheckd--;
		my $waitCnt = 0;
		while ( $num =~ m/$_/){
			print "Waiting for $jobsCheckd/".scalar @jobs ." jobs to finish\n" if ($waitCnt==0);
			sleep (30);
			$num = `$cmd`; #chomp $num;
			$waitCnt++;
			if ($killFailedJobs){
				my $killed = qsubDepNeverKill();
				print " Killed $killed jobs with Dependency never completed\n" if ($killed > 0);
				#die;
			}

		}
	}
	#print "returning\n";
	return;
}

sub qsubDepNeverKill{
	my $srchCmd = "squeue -u \$USER -t PENDING -o \"\%8i \%.15R \%17E\"  | grep 'ependencyNe' | cut -f1 -d' ' | xargs  -t -i scancel {} | wc -l ";
	my $num = 0;
	$num = `$srchCmd`; chomp $num;
	return $num;
	
}


sub qsubSystemWaitMaxJobs{
	my ($checkMaxNumJobs) = @_;
	my $killPend = $_[1] if (@_ > 1);
	my $optHR = {}; $optHR = $_[2] if (@_ > 2);
	
	return if ($checkMaxNumJobs <= 0);
	#my $srchCmd = "squeue |grep \$USER | grep PD |wc -l";
	my $num = numPendingJobs($optHR);
	my $waitCnt = 0;
	while ($num > $checkMaxNumJobs){
		if ($killPend){
			my $killed = qsubDepNeverKill();
			print " Killed $killed jobs with Dependency never completed\n" if ($killed > 0);
			
			#die;
		}
		print "waiting for jobs to finish (>$checkMaxNumJobs, qsubSystemWaitMaxJobs)...\n" if ($waitCnt == 0);
		sleep(40);
		$num =  numPendingJobs($optHR);;#`$srchCmd`; chomp $num;
		$waitCnt++;
		#print " $num ";
	}
	return;
}

sub qsubSystem2{
	my ($tmpsh,$optHR) = @_;
	my $hxr = $_[2] if (@_ > 2);
	my %xtras = %{$hxr}; 
	my $ncores = 0; 
	if (exists($xtras{cores})){$ncores = $xtras{cores};}
	my $nthreads= $ncores;
	if ($ncores =~ m/,/){my @spl = split /,/,$ncores;$ncores = $spl[1]; $nthreads=$spl[0];}
	if ($ncores != 0){#read in file, change it..
		open I,"<$tmpsh" or die "qsubSystem2: cant open $tmpsh\n";chomp(my @lines = <I>); close I;
		for (my $i=0;$i<@lines;$i++){
			if ($lines[$i] =~ m/--cpus-per-task/ || $lines[$i] =~ m/--mincpus/){
				$lines[$i] = "#SBATCH --cpus-per-task=$ncores";
				$lines[$i] .= "\n#SBATCH --threads-per-core=1\n#SBATCH --hint=compute_bound\n" unless ($lines[$i+1] =~ m/threads-per-core/);
			}
		}

	}
	my $xtra = "";
	my $qbin = "qsub";
	my $qmode = $optHR->{qmode};
	if ($qmode eq "slurm"){$qbin="sbatch";
	} elsif ($qmode eq "sge"){
	} else {$qbin="bsub";
	}
	my $qcm = "$qbin $xtra $tmpsh \n";
	die $qcm; #DEBUG
	system $qcm;
	return $qcm;
}
sub qsubSystem($ $ $ $ $ $ $ $ $ $){
	#args: 1[file to save bash & error & output] 2[actual bash cmd] 3[cores reserved for job]
	# 4["1G": Ram usage per core in GB] 5[0/1: synchronous execution] 6[name of job] 
	# 7[name of job dependencies, separated by ";"]
	# 8[0/1: excute in cwd?] 9[0/1: return qsub cmd or submit job to cluster]
	# Falk Hildebrand, may 2015
	my ($tmpsh,$cmd,$ncores,$memory,$jname,$waitJID,$cwd,$immSubm, $restrHostsAR, $optHR) = @_;
	#$doSync, 5th arg
	#14,12G
	#die $tmpsh."\n";
	#my $jname = $tmpsh;
	#$jname =~ s/.*\///g;$jname =~ s/\.sh$//g;
	#\n#\$ -N $tmpsh
	return("") if ($cmd eq "");
	my $LSF = 0;
	my $qbin = "qsub";
	my $xtra = "";
	my $rTag = $optHR->{rTag};
	my $qmode = $optHR->{qmode};
	
	my $tmpScratchTag = $optHR->{tmpSpaceTag};
	#my $xtraNodeCmds = $optHR->{xtraNodeCmds};
	my $submissionConfig = $optHR->{submissionConfig};
	my @constrains = @{$optHR->{constraint}};# #SBATCH --constraint=
	#die "@constrains";
	my $lockFile = $optHR->{LOCKfile};
	my $nthreads= $ncores;
	if ($ncores =~ m/,/){my @spl = split /,/,$ncores;$ncores = $spl[1]; $nthreads=$spl[0];}
	#different format for bsub and slurm
	if ($memory =~ m/^[\.\d]+$/){$memory  = int($memory+0.5);}
	if ($memory =~ m/^0G$/){$memory  = "1G";} #most likely a rounding error from too many cores..
	if ($memory =~ s/G$//){$memory = int( ($memory* 1024 * $ncores ) +0.5);};
	my $tmpSpace = convert2Gb( $optHR->{tmpSpace} );
	#die " $optHR->{tmpSpace}   $tmpSpace\n";
	#my $tmpSpace2 = $optHR->{tmpMinG};
	#my $wcKeysForJob = $optHR->{wcKeysForJob};
	my $exclNodes = $optHR->{excludeNodes};
	
	
	#die ($memory."\n");
	#my $queues = "\"".$optHR->{shortQueue}."\"";#"\"medium_priority\"";
	my $queues = "\"".$optHR->{medQueue}."\"";#"\"medium_priority\"";
	my $time = $optHR->{medTime};#"24:00:00";
	if ($optHR->{useHiMemQueue} == 1){
		$queues = "\"".$optHR->{highMemQueue}."\"";$optHR->{useHiMemQueue}=0;
	} elsif ($optHR->{useLongQueue} ==1){
		$queues = "\"".$optHR->{longQueue}."\"";#"\"medium_priority\"";
		#$time = "335:00:00";
		$optHR->{useLongQueue}=0;
	} elsif ($optHR->{useGPUQueue} ==1){
		$queues = "\"".$optHR->{gpuQueue}."\"";#"\"medium_priority\"";
		#$time = "23:00:00";
		$optHR->{useGPUQueue}=0;
	} elsif ($optHR->{useShortQueue} ==1){
		$queues = "\"".$optHR->{shortQueue}."\"";#"\"medium_priority\"";
		#$time = "00:45:00";
		$optHR->{useShortQueue}=0;
	}
	my @jspl = split(";",$waitJID); @jspl = grep /\S/, @jspl;

	if ($cwd ne "" && !-d $cwd){system "mkdir -p $cwd";}
	#if ($memory > 250001){$queues = "\"scb\"";}
	$tmpsh =~ m/^(.*\/)[^\/]+$/;
	system "mkdir -p $1" unless (-d $1);
	open O,">",$tmpsh or die "Can't open qsub bash script $tmpsh\n";
	#die "$cmd\n";
	#print "$memory   $queues\n";
	#if (`hostname` !~ m/submaster/){
	if ($qmode eq "slurm"){$LSF = 2;$qbin="sbatch";
		#if ($memory > 250001){$queues = "\"bigmem\"";}
		##SBATCH --cpus-per-task=$ncores\n
		print O "#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task=$ncores\n#SBATCH -o $tmpsh.otxt\n"; #\n#SBATCH -n  $ncores
		
		if ($nthreads != $ncores ){print O "#SBATCH --threads-per-core=1\n#SBATCH --hint=compute_bound\n";} #  specifically for iqtree/raxml
		print O "#SBATCH -e $tmpsh.etxt\n#SBATCH --mem=$memory\n#SBATCH --export=ALL\n";
		#print O "#SBATCH --kill-on-invalid-dep=yes\n";
		#print O "#SBATCH --tmp=$tmpSpace\n" if ($tmpSpace>0);#SBATCH --gres=ssd\n
		foreach my $subTerm ( split /;/, $submissionConfig){
			print O "#SBATCH $subTerm\n" if ($submissionConfig ne "");
		}
		if ($tmpSpace>0 && $tmpScratchTag ne ""){
			print O "#SBATCH $tmpScratchTag". int($tmpSpace+0.5) ."\n" ;
		}
		#"#SBATCH --gres=ssd"
		print O "#SBATCH -p $queues\n";
		#print O "#SBATCH --gres=tmp:${tmpSpace2}G\n" if ($tmpSpace2>0); #50g
		print O "#SBATCH --time=$time\n" unless ($time eq "");
		print O "#SBATCH --exclude=$exclNodes\n" unless ($exclNodes eq "");
		#print O "#SBATCH --localscratch=ssd:50\n"; #for EI cluster
		print O "#SBATCH --chdir=$cwd\n" if ($cwd ne "");
		print O "#SBATCH -J $rTag$jname\n" if ($jname ne "");
		print O "#SBATCH --wc=". $optHR->{wcKeysForJob} . "\n" if ($optHR->{wcKeysForJob} ne "");
		if (@constrains){
			print O "#SBATCH --constraint=". join(",",@constrains) ."\n" if (@constrains);
		}
		#foreach (@constrains){
	#		print O "#SBATCH --constraint=$_\n" if ($_ ne "");
		#}
		if (length($waitJID) >3 && @jspl > 0) {
			for (@jspl) {s/$rTag//;}
			#$xtra .= "--dependency=afterok:".join(":",@jspl)." " if (@jspl > 0);
			if ($optHR->{afterAny}){
				print O "#SBATCH --dependency=afterany:".join(":",@jspl)."\n" ;
			} else {
				print O "#SBATCH --dependency=afterok:".join(":",@jspl)."\n" ;
			}
			#use this one for now, as slurm currently faults without a reason..
			#$xtra .= "--dependency=afterany:".join(":",@jspl)." " if (@jspl > 0);
		}

		#print O "#\$ -S /bin/bash\n#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n";##\$ -v TMPDIR=/dev/shm\n";
		#print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n"; #causes problems..
	} elsif ($qmode eq "bash"){
		$qbin="bash";$LSF=3;
		print O "#!/bin/bash\n";
	} elsif ($qmode eq "sge"){
		print O "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n#\$ -pe ".$optHR->{qsubPEenv}." $nthreads\n#\$ -o $tmpsh.otxt\n#\$ -e $tmpsh.etxt\n#\$ -l h_rss=$memory\n";#h_vmem=$mem\n";
		print O "#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n";#\$ -v TMPDIR=/dev/shm\n";
#		print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n";
		print O "#\$ -V\n";
	} else {
		$LSF = 1;$qbin="bsub";
		print O "#!/bin/bash\n";
		print O "export LD_LIBRARY_PATH=/g/bork3/home/hildebra/env/env/miniconda/lib/:/g/bork3/home/hildebra/env/zlib-1.2.8/:/g/bork8/costea/boost_1_53_0/:/shared/ibm/platform_lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/g/bork3/x86_64/lib64:/g/bork3/x86_64/lib:\${LD_LIBRARY_PATH}\n\n";
		#print O "export LD_LIBRARY_PATH=/g/bork3/home/hildebra/env/zlib-1.2.8:/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64:\${LD_LIBRARY_PATH}\n\n";#:/g/software/linux/pack/python-2.7/lib/\nexport PATH=/g/bork3/home/zeller/py-virtualenvs/py2.7_bio1/bin/:\${PATH}\n\n";
		##BSUB -n $ncores\n#BSUB -o $tmpsh.otxt\n#BSUB -e $tmpsh.etxt\n#BSUB -M $mem\n#\$ -v LD_LIBRARY_PATH=/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64\n#\$ -v TMPDIR=/dev/shm\n#BSUB -q medium_priority\n";
		my @restrHosts = @{$restrHostsAR};
		if ( @restrHosts > 0){
			$xtra .= " -m \"".join(" ",@restrHosts)."\" ";
			$queues = "\"medium_priority scb\"";
		}
		$xtra .= "-n $nthreads -oo $tmpsh.otxt -eo $tmpsh.etxt -q $queues -M $memory -R \"select[(mem>=$memory)] ";
		$xtra .= "rusage[tmp=$tmpSpace] " if ($tmpSpace>0);
		$xtra .= "span[hosts=1]\" -R \"rusage[mem=$memory]\" "; #
	}
	#set abortion on program fails
	print O "echo \$HOSTNAME;\n";
	print O "set -eo pipefail\n";
	print O "ulimit -c 0;\n";
	#any xtra commands (like module load perl?)
	print O "$optHR->{xtraNodeCmds}\n";
	#prevent core dump files
	#file location check availability
	#print O $optHR->{LocationCheckStrg};

	print O $cmd."\n";
	close O;
	#sleep (1);
	my $depSet=0;
	if ($LSF==2){#slurm
		if ($optHR->{doSync} == 1){$qbin = "srun";}
		
	} elsif ($LSF==3){ #bash
		$xtra = "";
	} elsif ($LSF==1){ #bsub #-M memLimit; -q queueName;  -m "host_name[@cluster_name]; -n minProcessors; 
		if ($optHR->{doSync} == 1){$xtra.="-K ";}
		if ($jname ne ""){$xtra.="-J $rTag$jname ";}
		if (length($waitJID) >3) {
			my @jspl = split(";",$waitJID);
			#remove empty elements
			@jspl = grep /\S/, @jspl;
			if (@jspl > 0 ){
				$waitJID = join(") && done(",@jspl);
				$xtra.="-w \"done($waitJID)\" ";
			}
		}
		$tmpsh = " < ".$tmpsh;
	} else{ #qsub
		if ($optHR->{doSync} == 1){$xtra.="-sync y ";}
		if ($jname ne ""){$xtra.="-N $rTag$jname ";}
		if (length($waitJID) >3) {
			if (@jspl > 0 ){$xtra.="-hold_jid ".join(",",@jspl) ." ";}
		}
			#$waitJID =~ s/;/,/g;$xtra.="-hold_jid $waitJID ";}
	}
	if ($cwd ne ""){if ($LSF==1) {$xtra.="-cwd $cwd"; }  elsif ($LSF == 0) {$xtra.="-wd $cwd";} }
	my $qcm = "$qbin $xtra $tmpsh \n";
	my $LOGhandle = "";
	if (exists $optHR->{LOG}){ $LOGhandle = $optHR->{LOG};}
	#if (@restrHosts > 0){die $qcm;}
	if ($optHR->{doSubmit} != 0 && $immSubm){
		system "rm -f $tmpsh.otxt $tmpsh.etxt";
		print $LOGhandle $qcm."\n" unless ($LOGhandle eq "" || !defined($LOGhandle) );
		#print("$qcm\n\n");
		print "SUB:$jname\t";
		#actual job excecution!
		my $ret = `$qcm`; 
		#take care of lockFile now.. but only if actual job submission happened
		if ($lockFile ne "" && $ret !~ m/^sbatch: error:/ && ! -e $lockFile){
			system "touch $lockFile" ;
		}
		if ($LSF == 2){#slurm get jobid
			chomp $ret; $ret =~ m/(\d+)$/; #$ret = $1;
			$jname=$1;
		}
	}
	
	#die "$qcm\n";
	my $retJName = "$rTag$jname"; $retJName = "" if (!$immSubm); #return empty (for slurm), since no fwd job predictions..

	return ($retJName,$qcm);
}


#handles deleting of lock file, if all jobs have finished for current sample
sub MFnext($ $ $ $){
	my ($lckFile,$aR,$Jnum,$QSBoptHR) = @_;
	return if (! @{$aR});
	my $logF = $lckFile; $logF =~ s/\/[^\/]+$/rmLock.sh/;
	my $cmd = "echo \"all smpl associated jobs seem to have quit. Releasing lock..\"\nrm -f $lckFile\n";
	#my @jobs = @{$aR};
	my $jDepe = join(";",@{$aR});
	my $jobN = "RMLCK$Jnum";
	#print "$logF\n$jDepe\n\n"; 
	$QSBoptHR->{afterAny}=1;
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
	$QSBoptHR->{useShortQueue} =1;
	qsubSystem($logF,$cmd,1,"1G",$jobN,$jDepe,"",1,\{},$QSBoptHR);
	$QSBoptHR->{afterAny}=0;$QSBoptHR->{useShortQueue}=0;
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
}


sub add2SampleDeps($ $){
	my ($ar1, $ar2) = @_;
	foreach (@{$ar2}){
		push (@{$ar1}, $_) if (defined $_ && $_ ne "" && $_ =~ m/[^;]/);
	}
}