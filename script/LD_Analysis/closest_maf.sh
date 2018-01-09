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
# VCF_MAF_SCRIPT=/Users/chaochih/GitHub/Env_Assoc/script/LD_Analysis/vcf_maf_snp_id.py
# COMP_SNPS_DIR=/Users/chaochih/Dropbox/test_files/test_maf_ld/ld_results
# INT_VCF_DIR=/Users/chaochih/Dropbox/test_files/test_maf_ld/extracted_window
# PREFIX=Chr1-7_
# OUT_DIR=/Users/chaochih/Dropbox/test_files/test_maf_ld/test

#   User provided arguments
VCF_MAF_SCRIPT=$1 # Full path to vcf_maf_snp_id.py script
COMP_SNPS_DIR=$2 # Directory called ld_results containing *_compatibleSnps.txt
INT_VCF_DIR=$3 # Directory called extracted_window containing *_intersect.vcf
PREFIX=$4 # This is the part that comes before the SNP ID
OUT_DIR=$5 # Full path to our output directory


#   Check if out directory exists, if not make it
mkdir -p "${OUT_DIR}"
#   Create sample list of compatible SNPs and store in array
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
    #   If SNP is found in sample names file, remove file
    #   This means the target SNP is the same as the significant SNP
    if grep -q "${snp}" "${OUT_DIR}"/fake_target_compatible/"${snp}"_comp_snp_names_only.txt
    then
        echo "${snp}" >> "${OUT_DIR}"/fake_target_compatible/exists_snps.txt
        rm "${OUT_DIR}"/fake_target_compatible/"${snp}"_comp_snp_names_only.txt
    #   else, save SNP name to not_exists_snps.txt
    else
        echo "${snp}" >> "${OUT_DIR}"/fake_target_compatible/not_exists_snps.txt
    fi
done

#   Calculate MAF for all SNPs in each window
function calcMAF() {
    local vcf_maf_script=$1
    local comp_snps_file=$2
    local int_vcf_dir=$3
    local out_dir=$4
    #   Extract SNP from filename first
    snp=$(basename ${comp_snps_file} _comp_snp_names_only.txt)
    #   Extract compatible SNPs from intersect.vcf files
    #   since intersect.vcf files have not been filtered yet
    grep -f "${out_dir}"/fake_target_compatible/"${snp}"_comp_snp_names_only.txt "${int_vcf_dir}"/*"${snp}"_intersect.vcf > "${out_dir}"/fake_target_vcf/"${snp}"_comp_intersect.vcf

    #   Calculate MAF for SNPs
    python3 "${vcf_maf_script}" "${out_dir}"/fake_target_vcf/"${snp}"_comp_intersect.vcf > "${out_dir}"/fake_target_maf/"${snp}"_comp_intersect_vcf.maf
}

export -f calcMAF


#   Run functions
#   Make output directory first
mkdir -p "${OUT_DIR}"/fake_target_vcf
mkdir -p "${OUT_DIR}"/fake_target_maf
#   Calculate MAF for all SNPs in each window
parallel calcMAF "${VCF_MAF_SCRIPT}" {} "${INT_VCF_DIR}" "${OUT_DIR}" ::: "${OUT_DIR}"/fake_target_compatible/*comp_snp_names_only.txt
