#!/usr/bin/env bash

set -e
set -o pipefail

#   This script runs SNP_PhysPosition_Matching-csv.R on a list of filepaths

#   Usage: ./add_physPos.sh [SNP_PhysPosition_Matching-csv.R] [sample_list.txt] [vcf] [out_dir]
#   Where:
#   1) [SNP_PhysPosition_Matching-csv.R] is the full file path to this script
#   2) [sample_list.txt] is a list containing full file paths to gwas.csv files
#   3) [vcf] file that contains physical positions and updated chr
#   4) [out_dir] is the full path to our output directory


#   User provided arguments
SCRIPT=$1
SAMP_LIST=$2
VCF=$3
OUT_DIR=$4

#   Check if output directory exists, if not, make it
mkdir -p "${OUT_DIR}" "${OUT_DIR}"

SAMP_ARRAY=($(cat "${SAMP_LIST}"))

for i in "${SAMP_ARRAY[@]}"
do
    #   store just phenotype extracted from filename
    #   i.e. /Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/Fst/Results/Hierarchical_adegenet/FstLoci_out_long_withinGH_SNPorder.txt
    #   we would extract the "FstLoci_out_long_withinGH_SNPorder" part of the name in this case
    out_prefix=$(basename ${i} .txt)
    Rscript "${SCRIPT}" "${VCF}" "${i}" "${OUT_DIR}"/${out_prefix}_physPos.txt
done
