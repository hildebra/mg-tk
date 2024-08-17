#!/usr/bin/env bash
# Format eggNOG-mapper output in a way suitable for rtk input
# to allow abundance tables of functions to be generated
# The only argument to this script should be the location of
# the eggNOG-mapper output

# See https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Output_format
# for details of the eggNOG-mapper output format
OUTDIR=$(dirname $1)
# Run each sequentially rather than tee. Issues on HPC with tee-d processing.
cat $1 | cut -f1,19 | awk '$2~!/-/ {print $0}' > $OUTDIR/eggNOGmapper_CAZy.geneAss
#cat $1 | cut -f1,12,14 | awk '$2~!/-/ {gsub(/ko:/, ""); gsub(/-/, "", $3); print $0}' > $OUTDIR/eggNOGmapper_KGM.geneAss
#cat $1 | cut -f1,12,13 | awk '$2~!/-/ {gsub(/ko:/, ""); gsub(/-/, "", $3); print $0}' > $OUTDIR/eggNOGmapper_KGP.geneAss
cat $1 | cut -f1,12,14 | awk '$2~!/-/ {gsub(/ko:/, ""); gsub(/-/, "", $3); print $1"\t"$2";"$3}' > $OUTDIR/eggNOGmapper_KGM.geneAss
cat $1 | cut -f1,12,13 | awk '$2~!/-/ {gsub(/ko:/, ""); gsub(/-/, "", $3); gsub(/,?(map[0-9]*,?)/, "", $3); print $1"\t"$2";"$3}' > $OUTDIR/eggNOGmapper_KGP.geneAss
 
cat $1 | cut -f1,11 | awk '$2~!/-/ {print $0}' > $OUTDIR/eggNOGmapper_EC.geneAss
cat $1 | cut -f1,10 | awk '$2~!/-/ {print $0}' > $OUTDIR/eggNOGmapper_GO.geneAss
cat $1 | cut -f1,20 | awk '$2~!/-/ {print $0}' > $OUTDIR/eggNOGmapper_BIGG.geneAss
cat $1 | cut -f1,21 | awk '$2~!/-/ {print $0}' > $OUTDIR/eggNOGmapper_PFAM.geneAss
tail -n +2 $1 | cut -f1,5 | sed -E 's/^([0-9]*\t[A-Za-z0-9]*)@.*/\1/g' > $OUTDIR/eggNOGmapper_NOG.geneAss
