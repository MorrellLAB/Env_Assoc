#!/bin/bash

#PBS -l mem=22gb,nodes=1:ppn=16,walltime=03:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q small

set -e
set -u
set -o pipefail

#   Dependencies
module load bedtools_ML/2.23.0
module load parallel

#   User provided arguments
#   List of sorted BED files, these are the windows around the significant SNP
#   that we are interested in
BED_LIST=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.55_threshold/bed_win_list.txt
#   BED file prefix and suffix
BED_PREFIX=ld_Barley_NAM_200Kb_0.55_
BED_SUFFIX=_9k_masked_90idt_100000win.bed
#   BED file containing all transcript representatives
TRANSCRIPTS=/panfs/roc/groups/9/morrellp/llei/Envro_ass_landrace/Repr_Transcript_gtf/Repr_only_Transcripts.bed
#   Where do our output files go?
OUT_DIR=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.55_threshold/gene_intervals

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

#   Check if output directory exists, if not, make it
mkdir -p "${OUT_DIR}"
#   Store list of sorted BED files in an array
BED_ARRAY=($(cat /home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/extracted_window/bed_win_list.txt))

#   Call on function to extract transcript regions that overlap with significant SNP windows
parallel bedClosest {} "${BED_PREFIX}" "${BED_SUFFIX}" "${TRANSCRIPTS}" "${OUT_DIR}" ::: "${BED_ARRAY[@]}"
