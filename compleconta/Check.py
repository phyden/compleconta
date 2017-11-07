#!/usr/bin/env python

def check_genome_cc_weighted(marker_set, profile):
        # marker_set is not a list, but a EnogList object (check EnogList.py classes)

    marker_dict=marker_set.get_dict()

    found={}
    for item in profile:
        if marker_dict.get(item):
            num=found.get(item,0)
            found[item]=num+1

    found_sum=0
    multiple_sum=0
    for gene in found.keys():
        weight=marker_set.get_weight(gene)
        found_sum=found_sum+weight
        if found[gene]>1:
            multiple_sum=multiple_sum+(found[gene]-1)*weight
            

    tot=marker_set.get_total()
        
    completeness=float(found_sum)/tot

    #absolute contamination:
    contamination=float(multiple_sum)/tot

    return completeness, contamination
