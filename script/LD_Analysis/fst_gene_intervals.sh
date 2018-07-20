#!/bin/bash

#   This script stores all functions needed to generate gene intervals
#   for fst outliers LD decay plots.

set -e
set -o pipefail

function bedClosest() {
    local bed_file=$1
    local prefix=$2
    local suffix=$3
    local transcripts=$4
    local out_dir=$5
    #   Strip prefix and suffix from BED filename
    snp_name=$(basename "${bed_file}" | sed -e "s/^${prefix}//" -e "s/${suffix}$//")
    #   Add header line to output file
    printf "chr_int\tinterval_start\tinterval_end\tchr_gene\tgene_start\tgene_end\tdistance\n" > "${out_dir}"/"${snp_name}"_gene_intervals.txt
    #   Extract gene intervals overlapping significant SNP window
    bedtools closest -a "${bed_file}" -b <(cut -f 1,2,3 ${transcripts}) -D b | sort -k 7,7n -k 4,4n -k 5,5n | awk '$7==0' >> "${out_dir}"/"${snp_name}"_gene_intervals.txt
}

export -f bedClosest

function main() {
    local bed_list=$1
    local bed_prefix=$2
    local bed_suffix=$3
    local transcripts=$4
    local out_dir=$5
    #   Check if output directory exists, if not, make it
    mkdir -p "${out_dir}"
    #   Store list of sorted BED files in an array
    bed_array=($(cat "${bed_list}"))
    #   Call on function to extract transcript regions that overlap with significant SNP windows
    parallel bedClosest {} "${bed_prefix}" "${bed_suffix}" "${transcripts}" "${out_dir}" ::: "${bed_array[@]}"
}

export -f main
