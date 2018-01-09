#!/usr/bin/env bash

#   This script uses a Python script to calculate MAF for all SNPs
#       then finds the SNP with MAF that is closest to MAF of significant
#       SNP Fumi's GWAS analysis identified.
#   This script was created because the target SNPs (ideally the significant SNP)
#       used to create LD decay plots isn't always present due to things like
#       stringent filters used when creating the VCF. So we wanted to find the next SNP
#       that is closest in MAF to the fake target SNP.

#   Usage:
#       ./closest_maf.sh [arg1] [arg2]

set -e
set -o pipefail

#   For testing only, delete after done
VCF_MAF_SCRIPT=/Users/chaochih/GitHub/Env_Assoc/script/LD_Analysis/vcf_maf_snp_id.py
COMP_SNPS_DIR=/Users/chaochih/Dropbox/test_files/test_maf_ld/ld_results
INT_VCF_DIR=/Users/chaochih/Dropbox/test_files/test_maf_ld/extracted_window
OUT_DIR=/Users/chaochih/Dropbox/test_files/test_maf_ld/test
PREFIX=Chr1-7_

#   User provided arguments
VCF_MAF_SCRIPT=$1 # Full path to vcf_maf_snp_id.py script
COMP_SNPS_DIR=$2 # Directory called ld_results containing *_compatibleSnps.txt
INT_VCF_DIR=$3 # Directory called extracted_window containing *_intersect.vcf
PREFIX= # This is the part that comes before the SNP ID
OUT_DIR=


#   Check if out directory exists, if not make it
mkdir -p "${OUT_DIR}"
#   Created sample list of compatible SNPs and store in array
COMP_SNP_ARRAY=($(find "${COMP_SNPS_DIR}"/*_compatibleSnps.txt | sort))

#   Make directory for compatible snps that are the fake target SNPs
mkdir -p "${OUT_DIR}"/fake_target_compatible
#   Create a lists of compatible fake target SNP names only
for i in "${COMP_SNP_ARRAY[@]}"
do
    #   Suffix will always be _compatibleSnps.txt
    #   Prefix is defined above, so we can use sed to remove prefix
    #   which leaves just the SNP name
    snp=$(basename ${i} _compatibleSnps.txt | sed -e s/^${PREFIX}//)
    head -n 1 ${i} | tr '\t' '\n' > "${OUT_DIR}"/fake_target_compatible/"${snp}"_comp_snp_names_only.txt
    #   If snp exists in array, print message
    #   This method is quick and dirty and can get false positives
    if grep -q "${snp}" "${OUT_DIR}"/fake_target_compatible/"${snp}"_comp_snp_names_only.txt
    then
        echo "${snp}" >> "${OUT_DIR}"/fake_target_compatible/exists_snps.txt
        rm "${OUT_DIR}"/fake_target_compatible/"${snp}"_comp_snp_names_only.txt
    #   else, save array to file with new line delimiter
    else
        echo "${snp}" >> "${OUT_DIR}"/fake_target_compatible/not_exists_snps.txt
    fi
done

#   Calculate MAF for all SNPs in each window
function calcMAF() {
    local comp_snps_vcf_dir=$1
    local out_dir
}

export -f calcMAF


#   Run functions


