#!/usr/bin/perl
#deprecated, look at getMarkersMultiSmpl for most of this functionality
#ess100 / FMGs that have been binned are extracted from gene matrix, to get their abundance / sample
#further, the corresponding name in different samples is seperately stored
#further, extracts the contig names for TECs & contigs fastas for TEC/sample
#./findTECinCatalog.pl /g/scb/bork/hildebra/SNP/GNMass2_singl/alien-11-374-0/Binning/MaxBin/e100extr /g/scb/bork/hildebra/SNP/GNMass2_singl/
use warnings;
use strict;
sub combineClstr;
sub readClstrRev;
sub readMap;


my $tabixBin = "/g/bork5/hildebra/bin/samtools-1.2/tabix-0.2.6/./tabix";
my $bgzipBin = "/g/bork5/hildebra/bin/samtools-1.2/tabix-0.2.6/./bgzip";
my $samBin = "/g/bork5/hildebra/bin/samtools-1.2/samtools";

my $rfilesD = $ARGV[0];
my $inD = $ARGV[1];
#my $outD = $ARGV[2];
unless ($inD =~ m/[\/\\]$/){$inD .= "/";}
my $GCd = $inD."GeneCatalog/";
my $GCMversion = 0.1; if (-e "$GCd/version.txt"){$GCMversion = `cat $GCd/version.txt`; chomp $GCMversion;}

my $outD = $GCd."e100Substs/";
system "rm -rf $outD; mkdir -p $outD";

my ($hrm) = readMap($inD."LOGandSUB/inmap.txt");
my %map = %{$hrm}; #$map{MM3}
#die $map{MM3}{dir}."\n";

#combineClstr("$GCd/compl.incompl.95.fna.clstr","$GCd/Matrix.genes2rows.txt") unless (-e "$GCd/compl.incompl.95.fna.clstr.idx");

#system "$bgzipBin $GCd/Matrix.mat" unless (-e "$GCd/Matrix.mat.gz");
#system "$tabixBin -S 1 -s 1 $GCd/Matrix.mat.gz" unless (-e "$GCd/Matrix.mat.gz.tbi");

my ($hr,$hr2) = readClstrRev("$GCd/compl.incompl.95.fna.clstr.idx");
my %gene2clus = %{$hr};
my %clus2gene = %{$hr2};
#sed -n '55p;1110232p;2110232p'
my $e100inTECs = "";

my @TECNOS = ("TEC2","TEC3","TEC4","TEC5","TEC6");
my %TECcon;
foreach my $curTEC (@TECNOS){
	my $TECf = $rfilesD."/$curTEC.txt";
	my %srchs;
	open I,"<$TECf";
	while (my $line = <I>){
		chomp $line;
		next if (length($line) < 3);
		my @spl = split(/\t/,$line);
		$srchs{$spl[0]} = $spl[1];
	}
	close I;
	print "Read $curTEC info\n";


	my @srIdx ; my $geneSet="";
	my @curKeys = keys %srchs;
	for  (my $j=0;$j<@curKeys;$j++){
		 my $sr = ">".$curKeys[$j];
		unless (exists $gene2clus{$sr}){
			print "can't find gene $sr\n" ; next;
		}
		my $tblxIdx = $gene2clus{$sr};
		push (@srIdx, $tblxIdx);
		if (exists($clus2gene{$tblxIdx})){
			my @genes = split(/,/,$clus2gene{$tblxIdx});
			$geneSet .= join("\t$curTEC\n",@genes)."\t$curTEC\n";
			foreach (@genes){
				$_ =~ s/_\d+//;
				$TECcon{$curTEC}{$_} = 1;
			}
		} else {
			die "Can't find $tblxIdx in mat\n";
		}
	#	print "$sr $tblIdx\n";
	#	system "$tabixBin $GCd/Matrix.mat.gz $sr"
	}
	
	#extract fna / faa seq

	print "Selected ".@srIdx." rows.. writing to $outD/$curTEC.mat\n";
	system "$samBin faidx $GCd/compl.incompl.95.fna ".join(" ",@srIdx). " > $outD/$curTEC.e100.fna" ;
	system "$samBin faidx $GCd/compl.incompl.95.prot.faa ".join(" ",@srIdx). " > $outD/$curTEC.e100.prot.faa" ;

	push(@srIdx,1);
	open O,">$outD/$curTEC.lines";print O join("\n",@srIdx);close O;
	system "/g/bork5/hildebra/dev/C++/rare/rare lineExtr $GCd/Matrix.mat $outD/$curTEC.mat $outD/$curTEC.lines; rm $outD/$curTEC.lines";


	#add +2 as 1 is head and 0 is line 1 in mat
#	if ($GCMversion <= 0.1){
#		foreach my $ele (@srIdx){$ele+=2;} 
#	}
#	my $sedStr = join("p;",@srIdx);
#	system "sed -n '1p;$sedStr"."p' $GCd/Matrix.mat > $outD/$curTEC.mat";
	$e100inTECs .= $geneSet;
}

open O,">$outD/allTECe100.genes";
print O $e100inTECs;
close O;

#list of contigs
open O,">$outD/allTECe100.contigs";
foreach my $TE (keys %TECcon){
	print O join("\t$TE\n",keys %{$TECcon{$TE}})."$TE\n";
}
close O;

#print contigs fastas, separated by sample
foreach my $TE (keys %TECcon){
	my $tOdir = $outD.$TE."/";
	print $tOdir;
	system "mkdir $tOdir" unless (-d $tOdir);
	my @ctgs = sort(keys %{$TECcon{$TE}});
	#get ctgs from correct sample with faidx
	my $curSmpl = ""; my $curScaff=""; my $list="";
	foreach my $ctg (@ctgs){
		my @tmp = split(/__/,$ctg);
		if ($tmp[0] ne $curSmpl){
			if ($curScaff ne ""){
				print "$samBin faidx $curScaff $list > $tOdir$curSmpl.$TE.scaffs.fna";
				system "$samBin faidx $curScaff $list > $tOdir$curSmpl.$TE.scaffs.fna";
			}
			$curSmpl = $tmp[0]; $list="";
			#find in map
			die "not in map $curSmpl\n" unless (exists($map{$curSmpl}{dir}));
			#$map{$curSmpl}{dir}
			my $assD = `cat $inD/$map{$curSmpl}{dir}/assemblies/metag/assembly.txt`; chomp $assD;
			#print "ASSD ".$assD."\n";
			$curScaff = $assD."scaffolds.fasta.filt";
		}
		$list .= " '".$ctg."'";
		
		#$tmp[1];
	}
}

print "Finished extraction\n";










#-------------------------------------------------------------
sub readClstrRev(){ #version for my shortened index file
	my ($inF) = @_;
	my %retR; my %retF;
	print "Reading clstr.. \n";
	open I,"<$inF"; my $curCl="";
	while (<I>){
		chomp; 
		next if (m/^#/ || length($_) < 5);
		my $pos = index($_, "\t");
		$curCl = substr($_,0,$pos);
		my $rem = substr($_,$pos+1);
		#print $curCl."XX$rem\n";
		foreach (split /,/,$rem){
			$retR{$_} = $curCl;
		}
		$retF{$curCl} = $rem;

	}
	close I;
	print "done reading\n";
	return (\%retR,\%retF);
}
sub readClstrRevOld(){#old version for real cd hit files
#reads cd-hit outfile, but indicating the genes inside to the Cluster
	my ($inF) = @_;
	my %retR; my %retF;
	print "Reading clstr.. \n";
	open I,"<$inF"; my $curCl="";
	while (<I>){
		chomp; 
		if (m/^>(\S+)/){
			$curCl = $1;
		} else {
			#my @spl = split(/\t/);
			m/>(.*)\.\.\.\s*/;
			#die $1."\n";
			$retR{$1} = $curCl;
			if (exists($retF{$curCl})){
				$retF{$curCl}.= "##".$1;
			} else {
				$retF{$curCl} = $1;
			}
		}
	}
	close I;
	print "done reading\n";
	return (\%retR,\%retF);
}
sub combineClstr(){
	my ($clstr,$idx) = @_;
	die "deprecated combineClstr\n";
	print "Combining cluster strings.. \n";
	open I,"<$clstr"; open O,">$clstr.idx"; open C,"<$idx"; 
	my $chLine = <C>;
	while (my $line = <I>){
		chomp $line;
		if ($line =~ m/^>/){
			$chLine = <C>;
			my @spl = split(/\t/,$chLine);
			if ($spl[1] ne $line){
				die "Can;t match \n$line \n$spl[1]\n";
			}
			$line = ">$spl[0]";
			#die "FND :: $line\n";
		}
		print O $line."\n";
	}
	print "done $clstr.idx\n";
	close I; close O;close C;
}


sub readMap{
	my ($inF) = @_;
	my $cnt =-1; my %ret; my @order = ();
	my $dirCol = 1; my $smplCol = 0; my $rLenCol = -1; my $SeqTech = -1; 
	my $AssGroupCol = -1; my $EstCovCol = -1;
	my %agBP;
	#print $inF."\n";
	open I,"<$inF" or die "Can't open $inF\n";
	#AssGrps
	while (<I>){
		$cnt++;
		next if (m/^#/ && $cnt > 0 );
		chomp;	my @spl = split(/\t/);
		
		if ($cnt == 0){
			my %heads; my $ccnnt = 0; foreach(@spl){$heads{$_} = $ccnnt; $ccnnt++;}
			$dirCol = -1;if(exists($heads{Path})){$dirCol= $heads{Path};}
			$smplCol = -1;if(exists($heads{"#SmplID"})){$smplCol= $heads{"#SmplID"};}
			$AssGroupCol = -1;if(exists($heads{"AssmblGrps"})){$AssGroupCol= $heads{"AssmblGrps"};}
			
			next;} #maybe later check for col labels etc
		
		#die $spl[0]." ".$spl[1]."\n";
		$ret{$spl[0]}{dir} = $spl[$dirCol];
		$ret{$spl[0]}{SmplID} = $spl[$smplCol];
		if ($SeqTech >= 0) {$ret{$spl[0]}{SeqTech} = $spl[$SeqTech]; } else {$ret{$spl[0]}{SeqTech} = "";}
		if ($rLenCol >= 0){$ret{$spl[0]}{readLength} = $spl[$rLenCol];} else {$ret{$spl[0]}{readLength} = 0;}
		if ($EstCovCol >= 0){$ret{$spl[0]}{DoEstCoverage} = $spl[$EstCovCol];} else {$ret{$spl[0]}{DoEstCoverage} = 0;}
		if ($AssGroupCol >= 0){
			$ret{$spl[0]}{AssGroup} = $spl[$AssGroupCol];
			if (!exists($agBP{$spl[$AssGroupCol]}{CntAim})){$agBP{$spl[$AssGroupCol]}{CntAim}=0;}
			$agBP{$spl[$AssGroupCol]}{CntAim}++;
			#print $agBP{$spl[$AssGroupCol]}{CntAim}. " :$spl[$AssGroupCol]\n" ;
		} else {$ret{$spl[0]}{AssGroup} = $cnt; $agBP{$cnt}{CntAim}=0;}
		
		push(@order,$spl[0]);
	}
	#die();
	close I;
	$ret{opt}{smpl_order} = \@order;
	return (\%ret);
}

sub firstIdx(){

}