#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 2024

@author: Rafiga Masmaliyeva @RafigaM
@email: rmasmaliyeva@gmail.com
"""
import pysam
from collections import defaultdict
import numpy as np
import allel
import csv
import pandas as pd
import os
from pathlib import Path


def convert_array_to_scalar(arr):
    gt = str(arr[0]) + "/" + str(arr[1])
    return gt


def get_samples(input):
    samples = []
    for subdir, dirs, files in os.walk(input):
        for file in files:
            filepath = subdir + os.sep + file
            if filepath.endswith(".vcf"):
                samples.append(filepath)

    return samples


def generate_cohort(input):
    samples = get_samples(input)
    #if '.DS_Store' in samples: samples.remove('.DS_Store')
    my_dict = defaultdict(lambda: {"CHROM": [], "POS": [], "ID": [], "REF": [], "ALT": [], "QUAL": [], "FILTER": [], "INFO": [], "FORMAT": [], "sample": [], "num": 0})
    for vcf in samples:
        if Path(vcf).is_file():
            callset = allel.read_vcf(vcf, fields='*')
            l = len(callset['variants/POS'])
            for record in range(l):
                FT = callset['variants/FILTER_PASS'][record]
                if FT:
                    SV_type = str(callset['variants/SVTYPE'][record])
                    sv_len = callset['variants/SVLEN'][record]
                    #if SV_type == "DEL":
                    if SV_type == "DEL":
                        chrom = callset['variants/CHROM'][record]
                        pos = callset['variants/POS'][record]
                        ref = callset['variants/REF'][record]
                        alt = str(callset['variants/ALT'][record][0])
                        end = callset['variants/END'][record]
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
                        my_dict[variant_key]['FORMAT'] = 'GT:GQ:FT'
                        gt = convert_array_to_scalar(callset['calldata/GT'][record][0])
                        gq = callset['calldata/GQ'][record][0]
                        my_dict[variant_key]['sample'] = gt + ":" + str(gq) + ":" + "PASS"
                        my_dict[variant_key]['num'] += 1

    return my_dict


def write2vcf(output_file, converted_dict):
    # VCF file header
    vcf_header = [
        "##fileformat=VCFv4.2",
        "##source=custom_script",
        "##INFO=<ID=INFO,Number=1,Type=Integer,Description=\"Custom INFO field\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter Status\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"
    ]

    with open(output_file, 'w') as file:
        # Write header
        for line in vcf_header:
            file.write(line + "\n")
        for inner_dict in converted_dict:
            vcf_row = [
                inner_dict['CHROM'],
                inner_dict['POS'],
                inner_dict['ID'],
                inner_dict['REF'],
                inner_dict['ALT'],
                inner_dict['QUAL'],
                inner_dict['FILTER'],
                inner_dict['INFO'],
                inner_dict['FORMAT'],
                inner_dict['sample'],
                f"num={inner_dict['num']}",
            ]
            # Join the elements with a tab and write to the file
            file.write('\t'.join(map(str, vcf_row)) + "\n")


def write2bed(output_file, inner_dict):
    print(output_file)
    my_df = pd.DataFrame.from_dict(inner_dict)
    my_df.to_csv(output_file, sep='\t', index=False, header=False)


def main():
    folder = "/path/to/the/folder"
    output_vcf = "cohort.vcf"
    output_bed = "cohort.bed"
    
    my_dict = generate_cohort(folder)
    normal_dict = {key: value for key, value in my_dict.items()}
    df = pd.DataFrame.from_dict(normal_dict, orient='index')
    converted_dict = df.to_dict(orient='records')

    write2vcf(output_vcf, converted_dict)
    write2bed(output_bed, converted_dict)




if __name__ == "__main__":
    main()
