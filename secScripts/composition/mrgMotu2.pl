#!/usr/bin/perl
#helper script to merge motu2 tables..

use Mods::GenoMetaAss qw(gzipopen  systemW);

my $inD = $ARGV[0];
my $mOTU2ComplCnts = $ARGV[1];
my $outD = $inD;
$outD =~ s/[^\/]+\/?$//;

my $m2mrgSto = "$outD/m2.Smpl.cnts.stone";

my @taxLvlN = ("kingdom","phylum","class","order","family","genus","species");

opendir(DIR, $inD) or die "Can't find motu dir: $inD\n";	
	my @m2f = sort ( grep { /.*\.motu2\.tab\.gz/  && -e "$inD/$_"} readdir(DIR) );	rewinddir DIR;
	#die "@m2f\n";
	print "Reading ". @m2f ." motus files..\n";
	my %mat; my %specCnt; my %tax; my @smpls; my %matP; my %taxCnt;
	foreach my $f (@m2f){
		die "no smpl id in file $f\n" unless ($f =~ m /(^.*)\.motu2\.tab\.gz/);
		my $smpl = $1;
		push @smpls,$smpl;
		my ($FH,$OK) = gzipopen("$inD/$f","motu2 file",1);
		while (my $l = <$FH>){
			next if ($l =~ m/^#/);
			chomp $l;
			my @spl = split /\t/,$l;
			next if ($spl[2] eq "0");
			$specCnt{$spl[0]} += $spl[2];
			my $taxo = $spl[1];
			$mat{$spl[0]}{$smpl} = $spl[2];
			#only needs to be done once
			$taxo =~ s/\|/;/g;
			$taxo =~s/[kpcofgs]__//g;
			next if (exists $tax{$spl[0]});
			$tax{$spl[0]} = $taxo;
			my @lphy = split /;/,$taxo;
			@lphy = ("?","?","?","?","?","?","?") if ($taxo eq "-1" || $taxo eq "unassigned");
			my $ntax = "";
			die "mrgMotu2:: not enough tax levels for entry: @lphy , motu $spl[0]\n" if (@lphy != 7 );
			for (my $lvl =0; $lvl < 7; $lvl ++){
				my $llphy = $lphy[$lvl]; $llphy = "?" if ($llphy eq "");
				$ntax .= ";" if ($lvl >0); $ntax .= $llphy; 
				$matP{$lvl}{$ntax}{$smpl} += $spl[2];
				$taxCnt{$lvl}{$ntax} += $spl[2];
			}
				#die "@lphy";
		}
		close $FH;
	}
	print "Read source motus files, writing matrix..\n";
#motu matrix
	open O,">$outD/m2.motu.txt";
	foreach my $sm (@smpls){print O "\t$sm";}
	print O "\n";
	foreach my $ta (sort { $specCnt{$b} <=> $specCnt{$a} } keys %specCnt) {
		next if ($specCnt{$ta} <= 0);
		print O $ta;
		my $tmpStr = "";
		foreach my $sm (@smpls){if (exists($mat{$ta}{$sm})){$tmpStr .= "\t".$mat{$ta}{$sm};} else {$tmpStr .= "\t0";}}
		
		print O $tmpStr;
		print O "\n";
	}
	close O;
	print "Writing higher level matrices..\n";
#higher tax, one mat each
	for (my $lvl =0; $lvl < 7; $lvl ++){
		open O,">$outD/m2.$taxLvlN[$lvl].txt";
		print O $taxLvlN[$lvl];
		foreach my $sm (@smpls){print O "\t$sm";}
		print O "\n";
		foreach my $ta (sort { $taxCnt{$lvl}{$b} <=> $taxCnt{$lvl}{$a} } keys %{$taxCnt{$lvl}}) {
			next if ($taxCnt{$lvl}{$ta} <= 0);
			print O $ta;
			foreach my $sm (@smpls){
				if (exists($matP{$lvl}{$ta}{$sm})){
					print O "\t".$matP{$lvl}{$ta}{$sm};
				} else {
					print O "\t0";
				}
			}
			print O "\n";
		}
		close O;
	}
	system "echo \"$mOTU2ComplCnts\" > $m2mrgSto";
	print  "motu2 tax tables are in: $outD\n";