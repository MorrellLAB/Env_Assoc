#!/usr/bin/env bash

set -e
set -o pipefail

#   This script runs SNP_PhysPosition_Matching-csv.R on a list of filepaths

#   Usage: ./add_physPos.sh [sample_list.txt] [vcf] [out_dir]
#   Where:
#   1) sample_list.txt is a list containing full file paths to gwas.csv files
#   2) vcf file that contains physical positions and updated chr
#   3) out_dir is the full path to our output directory


#   User provided arguments
SAMP_LIST=$1
VCF=$2
OUT_DIR=$3

#   Check if output directory exists, if not, make it
mkdir -p "${OUT_DIR}" "${OUT_DIR}"

SAMP_ARRAY=($(cat "${SAMP_LIST}"))

for i in "${SAMP_ARRAY[@]}"
do
    #   store just phenotype extracted from filename
    #   i.e. /Users/chaochih/Downloads/GWAS.Results2/GAPIT..altitude.GWAS.Results.csv
    #   we would extract the "altitude" part of the name in this case
    out_pheno_name=$(echo $(basename ${i}) | awk -F '[.]' '{ print $3 }')
    Rscript SNP_PhysPosition_Matching-csv.R "${VCF}" "${i}" "${OUT_DIR}"/GAPIT_${out_pheno_name}_GWAS_results_physPos.csv
done
