#!/bin/bash

set -e
set -o pipefail

#   Define usage message
function usage() {
    echo -e "\
$0: \n\
\n\
This script is used to prepare genotyping data for LDheatmap.R script. \n\
\n\
Usage: ./LD_data_prep.sh [snp_bac] [genotype_data] [out_prefix] [work_dir] [extraction_SNPs.pl] \n\
\n\
NOTE: arguments must be provided in this order! \n\
\n\
where: \n\
1. [snp_bac] is a file that includes the Query_SNP name and physical position \n\
2. [genotype_data] is a file that has genotype data \n\
3. [out_prefix] is what prefix will our output filename look like? \n\
4. [work_dir] is where we want to output our files? \n\
5. [extraction] is where our extraction_SNPs.pl script is located \n\
" >&2
exit 1
}

if [[ $# -lt 3 ]]; then usage; fi

if ! $(command -v perl > /dev/null 2> /dev/null); then echo "Failed to find Perl, exiting..." >&2; exit 1; fi

#   Generate sample lists
function sampleList() {
    local snp_bac=$1
    local out_prefix=$2
    local work_dir=$3
    #   Sort SNP_BAC.txt file by Query_SNP columns
    #   Only keep unique SNP names
    #   Process GBS SNPs with naming format S1H1_71079
    echo "Generating sample lists..."
    (head -n 1 "${snp_bac}" && tail -n +2 "${snp_bac}" | sort -u -k2n,2) > "${work_dir}/SNP_BAC_${out_prefix}_sorted_uniq.txt"
    #   Generate sample lists
    for i in $(awk '{ print $1 }' "${work_dir}/SNP_BAC_${out_prefix}_sorted_uniq.txt" | tail -n +2); do echo "$i"; done > "${work_dir}/${out_prefix}_sampleList.txt"
}

export -f sampleList

#   extraction_SNPs.pl to pull out markers of interest
#   Arguments required: (1) WBDC_genotype_count.txt, (2) Sample_names_list.txt
function extractSNPs() {
    local extraction=$1
    local geno=$2
    local sample_names=$3
    local out_prefix=$4
    local work_dir=$5
    local exists=$6
    #   Use Perl script to pull down markers of interest
    #   and create new dataframe
    echo "Pulling down markers of interest..."
    #   Markers in both genotype dataframe and sample list
    "${extraction}" "${geno}" "${sample_names}" | awk '!/NOT_EXISTS/' > "${work_dir}/${out_prefix}_EXISTS.txt"
    #   Markers that are not in genotype dataframe are redirected to new file
    "${extraction}" "${geno}" "${sample_names}" | awk '/NOT_EXISTS/' > "${work_dir}/${out_prefix}_NOT_EXISTS.txt"
    #   How many SNPs were not found in dataframe?
    NOT_IN_DF=$(wc -l "${work_dir}/${out_prefix}_NOT_EXISTS.txt")
    echo "${NOT_IN_DF} SNPs did not exist in dataframe"
}

export -f extractSNPs

#   Sort data prior to feeding to LDheatmap.R
function sortData() {
    local exists=$1
    local out_prefix=$2
    local snp_bac=$3
    local work_dir=$4
    #   Sort dataframe output from extractSNPs
    #   with SNPs as rows and individual names as columns
    #   sort -V allows alphanumeric sort
    echo "Sorting output data..."
    awk 'NR<2{ print $0;next }{ print $0 | "sort -V" }' "${exists}" > "${work_dir}/${out_prefix}_sorted_EXISTS.txt"
    echo "Creating temporary sample list of existing markers..."
    #   Create tmp sample list from *_EXISTS.txt
    awk '{ print $1 }' "${work_dir}/${out_prefix}_sorted_EXISTS.txt" | tail -n +2 > "${work_dir}/tmp_${out_prefix}_filtered_sampleList.txt"
    echo "Creating filtered SNP_BAC.txt file"
    #   Create header for output file
    head -n 1 "${work_dir}/SNP_BAC_${out_prefix}_sorted_uniq.txt" > "${work_dir}/SNP_BAC_${out_prefix}_filtered.txt"
    #   Pull out markers that actually exist from SNP_BAC.txt
    grep -f "${work_dir}/tmp_${out_prefix}_filtered_sampleList.txt" "${work_dir}/SNP_BAC_${out_prefix}_sorted_uniq.txt" >> "${work_dir}/SNP_BAC_${out_prefix}_filtered.txt"
    echo "Cleaning up..."
    #   Cleanup unnecessary files
    rm "${work_dir}/${out_prefix}_sampleList.txt"
    rm "${work_dir}/tmp_${out_prefix}_filtered_sampleList.txt"
    rm "${work_dir}/SNP_BAC_${out_prefix}_sorted_uniq.txt"
    rm "${work_dir}/${out_prefix}_EXISTS.txt"
    echo "Done."
}

export -f sortData

#   User provided arguments
SNP_BAC=$1 # snpBAC.txt pulled from HarvEST
GENO=$2 # Genotype Data
OUT_PREFIX=$3 # output file prefix
WORK_DIR=$4 # where is our working directory?
#   Script filepaths
#   Please provide full filepaths to scripts
EXTRACTION=$5

#   Generate sample lists
sampleList "${SNP_BAC}" "${OUT_PREFIX}" "${WORK_DIR}"
#   Outputs dataframe with SNPs as columns and sample names as rows
extractSNPs "${EXTRACTION}" "${GENO}" "${WORK_DIR}/${OUT_PREFIX}_sampleList.txt" "${OUT_PREFIX}" "${WORK_DIR}"
#   Sort dataframes and pull out unique SNPs
sortData "${WORK_DIR}/${OUT_PREFIX}_EXISTS.txt" "${OUT_PREFIX}" "${WORK_DIR}/SNP_BAC_${OUT_PREFIX}_sorted_uniq.txt" "${WORK_DIR}"
