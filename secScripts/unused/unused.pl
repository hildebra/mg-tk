	if (0){ #consecutive assembly -- used later for gene catalog assembly?
		#--------------------------------------------------------------
		#Mapping to reduce target space
		if ($JNUM > 0){ 
			#map with bwt on assembly, to get reads relevant to project
			#dependency is $AssJob
			my $map2Ctgs = mapReadsToRef($mapOut,\@cfp1,\@cfp2,$tmpPath."/toMGctgs/",8,$prevAssembly,$AssJob,$tmpPath."unaligned/");#$localAssembly);
			$ifastp = $tmpPath."unaligned/un-conc-mate.1,".$tmpPath."unaligned/un-conc-mate.2";
			#die();
		}
		#--------------------------------------------------------------
		#Assembly
		if ($JNUM >= 0){ #start of the whole process with an initial assembly & subsequent assemblies of "missing bits"
			
	#		if(0){			$prevAssembly = spades Assembly(\@cfp1,\@cfp2,$metagAssDir,1,"","");
	#		} elsif(0) {			my($arp1,$arp2) = SEEECER(\@cfp1,\@cfp2,$tmpPath."seqSEEC/");			$prevAssembly = spades Assembly($arp1,$arp2,$metagAssDir,0,"","");		} elsif (1) {
			my($arp1,$arp2,$sdmjN) = sdmClean($ifastp,$metagAssDir."seqClean/",$jdep) if ($JNUM >= $continue_JNUM);
			#($arp1,$arp2) = SEEECER($arp1,$arp2,$tmpPath."seqSEEC/");
			my $finalD = "TODO";
			$AssJob = spadesAssembly($arp1,$arp2,$metagAssDir,$finalD,0,$sdmjN,$shortAssembly) if ($JNUM >= $continue_JNUM);
			if ($JNUM > 1){#somehow combine prevAssembly & metaGassembly 
				#TODO
			}
			$prevAssembly = $metaGassembly;
			$shortAssembly = $metaGassembly."2";
			#$metaGassembly = filterSizeFasta($metagAssDir."scaffolds.fasta",200);
			
			die() if ($JNUM == 1);
			#spade sAssembly(\@cfp1,\@cfp2,$metagAssDir);
		}
	} else {