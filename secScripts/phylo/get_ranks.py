import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()
target_ranks = ["superkingdom","phylum","class","order","family","genus","species"]
taxids = map(int, sys.argv[1:])
for tid in taxids:
    
    selected = [''] * len(target_ranks)
    try:
        ncbi.get_lineage(tid)
    except:
        continue
    for lin in ncbi.get_lineage(tid):
        if lin == 1: continue    
        rank = ncbi.get_rank([lin])[lin].lower()
        if rank in target_ranks:
            i = target_ranks.index(rank)
            selected[i] = ncbi.translate_to_names([lin])[0]
    print('\t'.join([str(tid)]+selected))
