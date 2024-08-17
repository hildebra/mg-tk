#!/usr/bin/perl
use warnings; use strict;

sub prepareSILVA;
#prepareSILVA("$ddir/SSURef_NR99_115_tax_silva.fasta",$DB,$DB2,"$ddir/SLVtax.csv");
my $ddir = "/g/bork5/hildebra/DB/silva_119/";
my $DB = "/g/bork5/hildebra/DB/silva_119/lts/SLV_SSU_119.fasta";
my $DB2 = "/g/bork5/hildebra/DB/silva_119/lts/SLV_SSU_119.tax";
my $DB3 = "/g/bork5/hildebra/DB/silva_119/lts/SLV_LSU_119.fasta";
my $DB4 = "/g/bork5/hildebra/DB/silva_119/lts/SLV_LSU_119.tax";

prepareSILVA("$ddir/SILVA_119_SSURef_Nr99_tax_silva.fasta",$DB,$DB2,"$ddir/tax_slv_ssu_nr_119.txt");
prepareSILVA("$ddir/SILVA_119_LSURef_tax_silva.fasta",$DB3,$DB4,"$ddir/tax_slv_lsu_119.txt");


sub prepareSILVA($ $ $ $){
#taxf3 is for 18S/28S #taxf3 is for SSU/LSU
my ($path, $SeqF,$taxF,$taxGuide) = @_;
print("Rewriting SILVA DB..\n");
my %taxG;
open I,"<",$taxGuide or die "Can't find taxguide file $taxGuide\n";
while (my $line = <I>){
	chomp($line); my @splg = split("\t",$line);
	$taxG{$splg[0]} = $splg[2];
} 
close I;
open I,"<",$path or die ("could not find SILVA file \n$path\n");
open OT,">",$taxF;
open OS,">",$SeqF;
#open OT2,">",$taxF2;open OS2,">",$SeqF2;
my @tdesign = (" k__"," p__"," c__"," o__"," f__"," g__"," s__");
my $skip = 0;
my $eukMode = 0;
while (my $line = <I>){
	chomp($line);
	if ($line =~ m/^>/){#header
		$skip=0;$eukMode = 0;
		my @spl = split("\\.",$line);
		if ($spl[0] =~ m/>AB201750/){
			$line = ">AB201750.1.1495 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 2;Anaerovirgula;Anaerovirgula multivorans";
			@spl = split("\\.",$line);
		} elsif ($spl[0] =~ m/>DQ643978/){
			$line = ">DQ643978.1.1627 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Geosporobacter;Geosporobacter subterraneus";
			@spl = split("\\.",$line);
		}elsif ($spl[0] =~ m/>X99238/){
			$line = ">X99238.1.1404 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 1;Thermobrachium;Thermobrachium celere";
			@spl = split("\\.",$line);
		} elsif ($spl[0] =~ m/>FJ481102/){
			$line = ">FJ481102.1.1423 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 1;Fervidicella;Fervidicella metallireducens AeB";
			@spl = split("\\.",$line);
		} elsif ($spl[0] =~ m/>EU443727/){
			$line = ">EU443727.1.1627 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;Thermotalea metallivorans";
			@spl = split("\\.",$line);
		}elsif ($spl[0] =~ m/>FR690973/){
			$line = ">FR690973.1.2373 Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae;Candidatus Thiopilula;Candidatus Thiopilula aggregata";
			@spl = split("\\.",$line);
		}elsif ($spl[0] =~ m/>CP002161/){
			$line = ">CP002161.5310.6845 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Candidatus Zinderia;Candidatus Zinderia insecticola CARI";
			@spl = split("\\.",$line);
		} elsif ($spl[0] =~ m/>FR690975/){
			$line = ">FR690975.1.2297 Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae;Candidatus Thiopilula;Candidatus Thiopilula aggregata";
			@spl = split("\\.",$line);
		}elsif ($spl[0] =~ m/>FR690991/){
			$line = ">FR690991.1.2147 Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae;Candidatus Thiopilula;Candidatus Marithioploca araucae";
			@spl = split("\\.",$line);
		}elsif ($spl[0] =~ m/>FR690991/){
			$line = ">AB910318.1.1553 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;uncultured bacterium";
			@spl = split("\\.",$line);
		}elsif ($spl[0] =~ m/>AB910318/){
			$line = ">AB910318.1.1553 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;uncultured bacterium";
			@spl = split("\\.",$line);
		}elsif ($spl[0] =~ m/>AY796047/){
			$line = ">AY796047.1.1592 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae 4;Thermotalea;uncultured bacterium";
			@spl = split("\\.",$line);
		}

		
		my $ID = $spl[0];
		$ID = substr($ID,1);
		$line =~ m/[^\s]+\s(.*)$/;
		my $tax = $1;
		if ($tax =~ m/^\s*Eukaryota/){$eukMode = 1;}#$skip = 1; next;}
		
		print OS ">".$ID."\n";
		@spl = split(";",$tax);
		my $tline;
		if (!$eukMode){
			if (@spl > 7 ){
				print $line."\n";
				print("too many categories\n");
			}
			for (my $i=0;$i<7; $i++){
				if ($i < scalar(@spl)){
					if ($spl[$i] =~ m/^unidentified/){$spl[$i] = "";}
					$spl[$i] = $tdesign[$i].$spl[$i];
				} else {
					$spl[$i] = $tdesign[$i];
				}
			}
			$tline = $ID ."\t".join(";",@spl);
		} else {#parse the levels out from taxguide
			my $tmpTax = "";
			my @jnd;
			my @soughtCls = ("domain","phylum","class","order","family","genus","species");
			my $soughtLvl = 0;
			for (my $i=0;$i<@spl; $i++){
				my $scanTax = $tmpTax.$spl[$i].";";
				if (exists($taxG{$scanTax}) || $soughtLvl == 6){
					#print "$taxG{$scanTax} LL\n";
					#SILVA has no species level in tax guide file
					if ($soughtLvl == 6){
						push(@jnd,$tdesign[$soughtLvl].$spl[$i]);
						$soughtLvl++;
						last;
					} elsif ($taxG{$scanTax} eq $soughtCls[$soughtLvl]){
						push(@jnd,$tdesign[$soughtLvl].$spl[$i]);
						#print $tdesign[$soughtLvl].$spl[$i]."\n";
						$soughtLvl++;
					} elsif ($taxG{$scanTax} eq $soughtCls[$soughtLvl+1]){#fill in empty levels
						push(@jnd,$tdesign[$soughtLvl]);
						$soughtLvl++;
						push(@jnd,$tdesign[$soughtLvl].$spl[$i]);
						#print "Skipped to level ".$tdesign[$soughtLvl].$spl[$i]."\n";
						$soughtLvl++;
					} elsif ($taxG{$scanTax} eq ""){#Euk in LSU file have no annotation..
						push(@jnd,$tdesign[$soughtLvl]."?".$spl[$i]);
						$soughtLvl++;
					}
				} else { #more likely to be low level species
					if ($spl[$i] =~ m/\S+\s\S+/ && int(@spl) >= ($i-1)){
						while ($soughtLvl<6){
							push(@jnd,$tdesign[$soughtLvl]."?");
							$soughtLvl++;
						}
						$soughtLvl = 6;
						push(@jnd,$tdesign[$soughtLvl]."?".$spl[$i]); 
						$soughtLvl++;
						last;
					} else {
						print $scanTax." JJ\n";
					}
				}
				 #Eukaryota;Fungi;Ascomycota;Archaeorhizomycetes;Archaeorhizomycetales;Archaeorhizomycetales_incertae_sedis
				$tmpTax .= $spl[$i].";";
			}
			for (;$soughtLvl<7;$soughtLvl++){
				push(@jnd,$tdesign[$soughtLvl]);
			}
			$tline = $ID."\t".join(";",@jnd);
			#die $tax." CC " .$tline."\n";
		}
		print OT $tline."\n";
		#die($tline);
	} elsif ($skip == 0){ #work through sequence
		$line =~ s/\s//g;
		$line =~ s/U/T/g;
		$line =~ s/u/t/g;
		#die $line;
		print OS $line."\n";
	}
}

close I; close OT; close OS; #close OT2; close OS2;
}
