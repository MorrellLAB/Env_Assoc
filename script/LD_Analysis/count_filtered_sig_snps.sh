#!/usr/bin/env bash

set -e
set -u
set -o pipefail

#   This script takes in a list of files containing filtered out SNPs and
#       counts the total number of significant SNPs that got filtered out
#       of LD analysis

#   Define usage message
function usage() {
    echo -e "\
$0: \n\
\n\
This script takes in a list of files containing filtered out SNPs and counts the total number of significant SNPs that got filtered out of LD analysis \n\
\n\
Usage: ./count_filtered_sig_snps.sh [filtered_snps_list.txt] [out_dir] \n\
\n\
NOTE: arguments must be provided in this order! \n\
\n\
where: \n\
1. [filtered_snps_list.txt] is a list of *_SNP_info-missing_data_cols.csv files \n\
2. [out_dir] the full filepath to our output directory \n\
" >&2
exit 1
}

if [[ $# -eq 0 ]]; then usage; fi

#   User provided arguments
MISS_LIST=$1 # list of *_SNP_info-missing_data_cols.csv files
OUT_DIR=$2 # full filepath to output directory

#   Store filepaths in a bash array
MISS_ARR=($(cat "${MISS_LIST}"))

#   Loop through each file in array
for file in "${MISS_ARR[@]}"
do
    #   Extract SNP name from filename
    snp_name=$(basename ${file} | sed -e 's/Chr1-7_//' -e 's/_SNP_info-missing_data_cols.csv//')
    if grep -Fq ${snp_name} ${file}
    then
        echo ${snp_name} "was filtered out."
        echo ${file} >> ${OUT_DIR}/sig_snps_filtered_out.txt
    else
        echo ${snp_name} "was not filtered out."
        echo ${file} >> ${OUT_DIR}/sig_snps_in_analyses.txt
    fi
done

echo "The number of significant SNPs that got filtered out of LD analyses is:"
wc -l ${OUT_DIR}/sig_snps_filtered_out.txt
