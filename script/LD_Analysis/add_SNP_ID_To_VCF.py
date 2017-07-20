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
                #   Store fields in a dictionary
                vcf_dat[chrom] = [pos, snp_id, ref, alt, qual, filt, info, form, genotypes]
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
                #   Store fields in a dictionary
                vcf_dat[chrom] = [pos, snp_id]
        return vcf_dat


def main(vcf):
    """Main function that runs the program."""
    vcf_to_edit = read_vcf_to_edit(vcf)
    vcf_with_names = read_vcf_w_names(vcf)
    return

main(sys.argv[1], sys.argv[2])
