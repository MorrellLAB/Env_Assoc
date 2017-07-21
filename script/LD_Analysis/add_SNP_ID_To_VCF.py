#!/usr/bin/env python3
#   Chaochih Liu - July 20, 2017

"""A script that replaces missing data in ID column of VCF file with SNP names from another VCF file. SNP names are added if CHROM and POS columns match in both VCF files provided.

Usage:
    add_SNP_ID_To_VCF.py [VCF_To_Edit.vcf] [VCF_with_SNP_names.vcf] > [outname.vcf]

Arguments:
    1) [VCF_To_Edit.vcf]
    2) [VCF_with_SNP_names.vcf]
"""

import sys

def read_vcf_to_edit(vcf):
    """Function that reads through file and skips any lines that start with '##'. Data will be stored in a dictionary."""
    vcf_dat = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            #   #CHROM line contains sample information
            elif line.startswith('#CHROM'):
                #   Split up the line on tabs
                tmp = line.strip().split('\t')
            else:
                #   Split up the lines on tabs
                tmp = line.strip().split('\t')
                #   Assign variables for clarity
                chrom = tmp[0]
                pos = tmp[1]
                snp_id = tmp[2]
                ref = tmp[3]
                #   If there are multiple alternate alleles, list them with a comma as the seperator
                alt = tmp[4].split(',')
                qual = tmp[5]
                filt = tmp[6] # FILTER column
                info = tmp[7]
                form = tmp[8] # FORMAT column
                genotypes = tmp[9:] # genotype fields
                #   Join chrom and pos
                chrom_pos = chrom + "_" + pos
                #   Store fields in a dictionary
                vcf_dat[chrom_pos] = [pos, snp_id, ref, alt, qual, filt, info, form, genotypes]
        return vcf_dat


def read_vcf_w_names(vcf):
    """Function that reads through file and skips any lines that start with '##'. Only first 3 fields will be stored in a dictionary."""
    vcf_dat = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            #   #CHROM line contains sample information
            elif line.startswith('#CHROM'):
                #   Split up the line on tabs
                tmp = line.strip().split('\t')
            else:
                #   Split up the lines on tabs
                tmp = line.strip().split('\t')
                #   Assign variables for clarity
                #   We are mostly interested in the SNP names (ID column)
                #   so we only need the first 3 fields from this vcf file
                chrom = tmp[0]
                pos = tmp[1]
                snp_id = tmp[2]
                #   Join chrom and pos
                chrom_pos = chrom + "_" + pos
                #   Store fields in a dictionary
                vcf_dat[chrom_pos] = [snp_id]
        return vcf_dat


def rename_id(v1, v2):
    #   Create empty dictionary
    renamed = {}
    for keys in v1.keys():
        #   If keys in v1 exist in v2,
        #   replace ID with SNP name
        if keys in v2.keys():
            #   store values in keys as list
            v1d = v1[keys]
            v2d = v2[keys]
            #   replace 2nd element (ID) with SNP name
            v1d[1] = v2d[0]
            #   add updated key-value pair to dictionary
            renamed[keys] = v1d
        #   If keys in v1 do not exist in v2,
        #   replace ID with chrom_pos naming scheme
        elif keys not in v2.keys():
            #   store values in keys as list
            v1d = v1[keys]
            v1d[1] = keys
            #   add updated key-value pair to dictionary
            renamed[keys] = v1d
    return renamed


def main(file1, file2):
    """Main function that runs the program."""
    vcf1 = read_vcf_to_edit(file1)
    vcf2 = read_vcf_w_names(file2)
    prd = rename_id(vcf1, vcf2)
    for key, val in list(prd.items()):
        print(
            key + '\t' + val[0] + '\t' + val[1] + '\t' + \
            val[2] + '\t' + val[3][0] + '\t' + val[4] + '\t' + \
            val[5] + '\t' + val[6] + '\t' + val[7] + '\t' \
            + '\t'.join(map(str, val[8]))
        )


main(sys.argv[1], sys.argv[2])
