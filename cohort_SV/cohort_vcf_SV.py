
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 2024

@author: Rafiga Masmaliyeva @RafigaM
@email: rmasmaliyeva@gmail.com
"""
import pysam
from collections import defaultdict
import numpy as np
#from fuc import pyvcf
import allel
import csv
import pandas as pd
import os
from pathlib import Path
import ast
#import matplotlib.pyplot as plt
#import seaborn as sns


from collections import Counter


def convert_array_to_scalar(arr):
    gt = str(arr[0]) + "/" + str(arr[1])
    return gt


def get_samples(input):
    #print(input)
    samples = []
    for subdir, dirs, files in os.walk(input):
        for file in files:
            filepath = subdir + os.sep + file
            if filepath.endswith(".vcf"):
                #print(filepath)
                samples.append(filepath)

    return samples


def find_unique_thresholds_with_indexes(thresholds, chr):
    unique_thresholds = []
    threshold_map = {}

    def has_significant_overlap(t1, t2):
        low1, high1 = int(t1[0]), int(t1[1])
        low2, high2 = int(t2[0]), int(t2[1])

        # Calculate the overlap
        overlap = max(0, min(high1, high2) - max(low1, low2))

        # Calculate percentage overlap for each interval
        percent_overlap_1 = overlap / (high1 - low1)
        percent_overlap_2 = overlap / (high2 - low2)

        # If either percentage is >= 20%, return True
        return percent_overlap_1 >= 0.8 or percent_overlap_2 >= 0.8

    for index, threshold in enumerate(thresholds):
        is_unique = True
        for unique_threshold in unique_thresholds:
            if has_significant_overlap(threshold, unique_threshold):
                # If it overlaps significantly, add it as a "mate" to the unique_threshold
                str_index = chr + ":" + str(threshold[0]) + ":" + str(threshold[1])
                threshold_map[tuple(unique_threshold)].append((str_index, threshold))
                is_unique = False
                break
        if is_unique:
            # If it's unique, add it to the unique_thresholds list
            unique_thresholds.append(threshold)
            # Initialize the map with the unique threshold as the key and the list of "mates" including its index
            str_index = chr + ":" + str(threshold[0]) + ":" + str(threshold[1])
            threshold_map[tuple(threshold)] = [(str_index, threshold)]

    return unique_thresholds, threshold_map


def generate_cohort(input):
    samples = get_samples(input)
    my_dict = defaultdict(
        lambda: {"CHROM": [], "POS": [], "ID": [], "REF": [], "ALT": [], "QUAL": [], "FILTER": [], "INFO": [], "sample": [], "num": []})
    for vcf in samples:
        if Path(vcf).is_file():
            callset = allel.read_vcf(vcf, fields='*')
            l = len(callset['variants/POS'])
            for record in range(l):
                FT = callset['variants/FILTER_PASS'][record]
                #imprecise = callset['variants/IMPRECISE'][record]
                if FT == True:
                    #if imprecise == False:
                    pos = callset['variants/POS'][record]
                    sv_len = callset['variants/SVLEN'][record]
                    end = pos + abs(sv_len)
                    #end = callset['variants/END'][record]
                    #sv_len = abs(int(pos) - int(end))  ###
                    if len(callset) > 0:
                        #sv_len = callset['variants/SVLEN'][record]
                        SV_type = str(callset['variants/SVTYPE'][record])
                        if SV_type == "DEL" and abs(sv_len) > 50:
                            chrom = callset['variants/CHROM'][record]
                            #pos = callset['variants/POS'][record]
                            ref = callset['variants/REF'][record]
                            alt = str(callset['variants/ALT'][record][0])
                            #end = callset['variants/END'][record]
                            variant_key = f"{chrom}:{pos}:{end}"
                            id = chrom + "_" + str(pos) + "_" + str(end)
                            qual = callset['variants/QUAL'][record]
                            # Dict to store all varinats
                            my_dict[variant_key]['CHROM'] = chrom
                            my_dict[variant_key]['POS'] = pos
                            my_dict[variant_key]['ID'] = id
                            my_dict[variant_key]['REF'] = ref
                            my_dict[variant_key]['ALT'] = alt
                            my_dict[variant_key]['QUAL'] = int(qual)
                            my_dict[variant_key]['FILTER'] = 'PASS'
                            end = callset['variants/END'][record]
                            svtype = callset['variants/SVTYPE'][record]
                            my_dict[variant_key]['INFO'] = "END=" + str(end) + ";SVTYPE=" + str(svtype) + ";SVLEN=" + str(sv_len) # + ";"
                            gt = convert_array_to_scalar(callset['calldata/GT'][record][0])
                            my_dict[variant_key]['sample'].append(gt)
                            vcf_name = vcf.split("/")[-1]
                            my_dict[variant_key]['num'].append(vcf_name)

    return my_dict

def write2vcf(output_file, converted_dict):
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=custom_script",
        "##INFO=<ID=AF,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##INFO=<ID=AC_het,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##INFO=<ID=AC_hom,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=FT,Number=1,Type=Float,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=GQ,Number=1,Type=String,Description=\"Filter Status\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"
    ]

    with open(output_file, 'w') as file:
        # Write header first
        for line in vcf_header:
            file.write(line + "\n")
        for inner_dict in converted_dict:
            vcf_row = [
                inner_dict['CHROM'],
                inner_dict['POS'],
                inner_dict['ID'],
                "REF",
                "<DEL>",
                inner_dict['QUAL'],
                inner_dict['FILTER'],
                str("AF=" + str(inner_dict['af'])+ ";AC=" + str(inner_dict['ac']) + ";AC_het=" + str(inner_dict['ac_het']) + ";AC_hom=" + str(inner_dict['ac_hom']) + ";END=" + str(inner_dict['end']) + ";SVTYPE=" + str(inner_dict['SVTYPE']) + ";SVLEN=" + str(inner_dict['SVLEN'])),
                "GT:FT:GQ",
                find_gt(inner_dict['GT']) + ":" + inner_dict['FILTER'] + ":" + str(inner_dict['QUAL']),
                #find_gt(inner_dict['GT']) + ":" + inner_dict['FILTER'] + ":" + str(int(inner_dict['QUAL'][0])),
            ]
            file.write('\t'.join(map(str, vcf_row)) + "\n")


def write2bed(output_file, inner_dict):
    print(output_file)
    my_df = pd.DataFrame.from_dict(inner_dict)
    my_df.to_csv(output_file, sep='\t', index=False, header=False)


def find_pos(mates):
    mates = list(mates)
    longest_threshold = 0
    longest_pair = None
    for start, end in mates:
        threshold = int(end) - int(start)
        if threshold > longest_threshold:
            longest_threshold = threshold
            longest_pair = (start, end)

    return(longest_pair)

def find_gt(gt_list):
    counts = Counter(gt_list)
    most_popular = counts.most_common(1)[0]

    return(str(most_popular[0]))


def write2bed2(output_file, converted_dict):
    with open(output_file, 'w') as file:
        for inner_dict in converted_dict:
            vcf_row = [
                inner_dict['CHROM'],
                inner_dict['POS'],
                inner_dict['end'],
                inner_dict['SVTYPE'],
                inner_dict['af']
            ]
            file.write('\t'.join(map(str, vcf_row)) + "\n")


def unique_sv(my_dict, n):
    n_samples = n
    dict_unique = defaultdict(
        lambda: {"CHROM": [], "POS": [], "end": [], "ID": [], "SVTYPE": [], "SVLEN": [], "FILTER": [], "QUAL": [],
                 "GT": [], "sample": [], "mates": [], "ac": [], "af": [], "ac_het": [], "ac_hom": []})
    chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    my_keys = list(my_dict.keys())
    for chr in chrs:
        poss = []
        for my_key in my_keys:
            #s = str.split(my_dict['ID'][j], "_")
            s = str.split(my_key, ":")
            if s[0] == chr:
                poss.append((int(s[1]), int(s[2])))
        unique_thresholds, threshold_map = find_unique_thresholds_with_indexes(poss, chr)
        mychr = list(threshold_map.items())

        for i in range(len(mychr)):
            #index = chr1[i][1][0][0]
            start= mychr[i][0][0]
            end = mychr[i][0][1]
            l_tmp = int(end) - int(start)
            if l_tmp < 1000000:
                key = chr + ":" + str(start) + ":" + str(end)
                dict_unique[key]['CHROM'] = chr
                #dict_unique[key]['ID'] = key
                dict_unique[key]['SVTYPE'] = "DEL"
                #dict_unique[key]['SVLEN'] = int(start) - int(end)
                dict_unique[key]['FILTER'] = "PASS"
                dict_unique[key]['af'] = 0
                dict_unique[key]['ac'] = 0
                dict_unique[key]['ac_hom'] = 0
                dict_unique[key]['ac_het'] = 0
                l = len(mychr[i][1])
                for ii in range(l):
                    index = mychr[i][1][ii][0]
                    dict_unique[key]['GT'].append(my_dict[index]['sample'])
                    dict_unique[key]['sample'].append(my_dict[index]['num'])
                    dict_unique[key]['mates'].append(mychr[i][1][ii][1])
                gt = dict_unique[key]['GT']
                gt = [item for sublist in gt for item in sublist]
                dict_unique[key]['GT'] = gt
                n_hom = 0
                n_het = 0
                for item in gt:
                    item = item.split("/")
                    if all(x in {'0', '1'} for x in item):
                        if int(item[0]) + int(item[1]) == 1:
                            n_het += 1
                        if int(item[0]) + int(item[1]) == 2:
                            n_hom += 1
                    else:
                        continue
                    dict_unique[key]['ac_het'] = n_het
                dict_unique[key]['ac_hom'] = n_hom
                af = (n_het + (2 * n_hom)) / (n_samples * 2)
                dict_unique[key]['af'] = af
                dict_unique[key]['ac'] = (n_het + (2 * n_hom))
                pos = find_pos(dict_unique[key]['mates'])
                dict_unique[key]['POS'] = pos[0]
                dict_unique[key]['end'] = pos[1]
                dict_unique[key]['SVLEN'] = abs(int(pos[0]) - int(pos[1]))
                dict_unique[key]['ID'] = chr + ":" + str(pos[0]) + "-" + str(pos[1])
                dict_unique[key]['QUAL'] = 999

    return(dict_unique)


def main():
    folder = "vcf_files"
    output_vcf = "cohort_del.vcf"
    output_bed0 = "cohort_del.bed"
    output_bed = "cohort_del_all.bed"
    my_dict = generate_cohort(folder)
    n = 1000 #number of samples in your cohort
    unique_dict = unique_sv(my_dict, n)
    unique_normal_dict = {key: value for key, value in unique_dict.items()}
    unique_df = pd.DataFrame.from_dict(unique_normal_dict, orient='index')
    unique_converted_dict = unique_df.to_dict(orient='records')
    #unique_converted_dict = unique_sv(my_dict)
    write2vcf(output_vcf, unique_converted_dict)
    write2bed(output_bed0, unique_converted_dict)
    write2bed2(output_bed, unique_converted_dict)


if __name__ == "__main__":
    main()
