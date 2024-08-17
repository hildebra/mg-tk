#Requires Bowtie2, Python2.7, MetaGeneMark, Augustus, USEARCH7.0
#Designed for a GNU/Linux parallel processing system, e.g., here, 16 processors and at least 72 GB of RAM. 
#binaries are in your own defined directory


#Call bacterial genes on files containing contigs following assembly.
for f in *.fa
do
./MetaGeneMark_linux64/gmhmmp -a -d -f G -m ./bin/MetaGeneMark_linux64/MetaGeneMark_v1.mod -o $f.mgm.gff $f;
perl ./MetaGeneMark_linux64/nt_from_gff.pl < $f.mgm.gff > $f.mgm.nt.fasta;
perl ./MetaGeneMark_linux64/aa_from_gff.pl < $f.mgm.gff > $f.mgm.aa.fasta;
done

#Call fungal genes (select your desired species. Here, Ustilago maydis is the closest relative to Malassezia)
for f in *.fa
do
./augustus --species=ustilago_maydis $f --protein=off --codingseq=on > $f.augustus.gff;
perl ./bin/getAnnoFast.pl --seqfile=$f $f.augustus.gff;
done

######################################
#to generate taxonomic assignments, follow workflow from /blasting/blasting.sh
#For each sample, get final *matched.txt files
######################################


#now, where there are genes called from both kingdoms, gene assignment is made by BLAST against the nr database. Where there was no resolution in the kingdom, genes were marked as ambiguous or assigned to whichever caller generated a prediction. 

nano combiningall.sh
for f in `cat listoffiles.txt` 
do 
grep "^NODE" $f*.mgm.gff > $f.mgm.list;
grep -w "transcript" $f*.augustus.gff > $f.aug.list;
grep ">" $f*_velvet_assembly.*.fa > $f.list;

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $f*.augustus.codingseq > $f.aug.nolines.fa;
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $f*.mgm.nt.fasta > $f.mgm.nolines.fa;

perl ./bin/process_nodes.pl $f.mgm.list $f.aug.list;

perl ./bin/match_nodes_euk.pl $f.aug.list.output.txt $f.aug.nolines.fa;
perl ./bin/match_nodes_bac.pl $f.mgm.list.output.txt $f.mgm.nolines.fa;

perl ./bin/overlap_nodes.pl $f.list $f.aug.list.output.txt $f.mgm.list.output.txt $f.blast.matched.txt;
done
##end
sh combiningall.sh


#For each file, have merged calls from bacterial and eukaryotic calls. Decision tree. This also gets some basic stats on the numbers. 
for f in `cat listoffiles.txt`
do
perl ./bin/finalgeneparsing.pl $f.list.ALLmatched.txt
cat $f.list.ALLmatched.txt.decision.txt | wc -l > $f.stats.txt
grep "ambig_2call" $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="unambig_nocall" && $5=="unambig_nocall") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="ambig_bac" && $5=="ambig_2call") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="ambig_euk" && $5=="ambig_2call") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="ambig_any" && $5=="ambig_2call") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="unambig_bac" && $5=="unambig_bac") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="ambig_euk" && $5=="unambig_bac") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="ambig_any" && $5=="unambig_bac") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="unambig_euk" && $5=="unambig_euk") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="ambig_bac" && $5=="unambig_euk") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
awk '{ if($3=="ambig_any" && $5=="unambig_euk") print $0}' $f.list.ALLmatched.txt.decision.txt | wc -l >> $f.stats.txt
cut -f 2 $f.list.ALLmatched.txt.decision.txt | sort | uniq -c > $f.kingdomcounts.txt
done

#now, let's pull the appropriated gene. for each contig, need to retrieve ALL genes for that contig.
for f in `cat listoffiles.txt`
do
cut -f 4 $f.list.ALLmatched.txt.decision.txt | sort |uniq> $f.formatching.txt
perl ./bin/getgenes.pl $f.formatching.txt $f.aug.nolines.fa.renamed.fa 
perl ./bin/getgenes.pl $f.formatching.txt $f.mgm.nolines.fa.renamed.fa
cat $f.aug.nolines.fa.renamed.fa.output.txt $f.mgm.nolines.fa.renamed.fa.output.txt > $f.genecalled.fa
done

#the terminal file here is 
cat *.genecalled.fa > allgenecalled.fa

###########NOTE############
#because the skin metagenome also includes not-bacteria and not-fungi (e.g., viruses), we opted to create our metagenomic clusters using all contings <1kb, and then calling genes on contigs >1kb so as minimize information loss. 

#get contigs < 1000bp
perl ./bin/filtering_fromassemblies.pl *_assembly.*.fa

#now, get genes from contigs >1000
perl ./bin/filtering_fromformatching.pl *.formatching.txt

#now, lets get the sequences
for f in `cat listoffiles.txt`
do
perl ./bin/getgenes.pl $f.formatching.txt.gt1000genes.txt $f.aug.nolines.fa.renamed.fa 
perl ./bin/getgenes.pl $f.formatching.txt.gt1000genes.txt $f.mgm.nolines.fa.renamed.fa
cat $f*velvet_assembly.*.fa.sub1000.fa $f.aug.nolines.fa.renamed.fa.output.txt $f.mgm.nolines.fa.renamed.fa.output.txt > $f.contig1000genes.genecalled.fa
done

#the terminal file here is 
cat *1000genes.genecalled.fa > allgenecalled.1000.fa

############################

#We then created a non-redundant catalog by clustering at a sequence identity cutoff of 0.95 and a minimum coverage cutoff of 0.9 for the shorter sequence. However, these are very stringent thresholds and it is possible that multiple fragments may in fact represent a single gene unit. This step can take ~100GB of RAM depending on your catalog size. 

#our gene-called files are now 
allgenecalled.fa
allgenecalled.1000.fa

./usearch7.0.1001_i86linux64 -cluster_fast allgenecalled.1000.fa -id 0.95 -idprefix 4 -target_cov .9 -consout allgenecalled.1000.centroids.fa -uc allgenecalled.1000.clusters.uc
./usearch7.0.1001_i86linux64 -cluster_fast allgenecalled.fa -id 0.95 -idprefix 4 -target_cov .9 -consout allgenecalled.centroids.fa -uc allgenecalled.clusters.uc

#a little cleanup
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' allgenecalled.1000.centroids.fa > allgenecalled.1000.centroids.nolines.fa
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' allgenecalled.centroids.fa > allgenecalled.centroids.nolines.fa

#and finally, make a lengths file for downstream use. 
perl ./bin/fastaNamesSizes.pl allgenecalled.1000.centroids.nolines.fa > allgenecalled.1000.lengths.txt
perl ./bin/fastaNamesSizes.pl allgenecalled.centroids.nolines.fa > allgenecalled.lengths.txt


#and so we end with our non-redundant gene database and an associated lengths file
allgenecalled.centroids.nolines.fa
allgenecalled.lengths.txt

#or 
allgenecalled.1000.centroids.fa
allgenecalled.1000.lengths.txt


