#!/bin/env bash

#PBS -l mem=22gb,nodes=1:ppn=16,walltime=16:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q lab

set -e
set -o pipefail

module load R/3.4.0
module load python2/2.7.8 # Tom's VCF_To_Htable-TK.py script runs in Python 2
module load vcflib_ML/1.0.0

#   Define usage message
function usage() {
    echo -e "\
    $0: \n\
    \n\
    This script performs LD analysis on 100Kb windows around significant SNPs from Environmental Associations project GWAS analysis. This script pulls together all parts of the analysis performed by multiple scripts. \n\
    \n\
    Please fill out user provided argument fields and submit script as job. The script will submit each of the list of significant SNP datasets as a task array. \n\
    " >&2
    exit 1
}

if [[ $# -lt 3 ]]; then usage; fi

#   User provided arguments
#   Where is the directory containing this script?
#   Make sure all scripts needed are in same directory as this script
SCRIPT_DIR=/home/morrellp/liux1299/GitHub/Env_Assoc/script/LD_Analysis
#   List of significant SNP names based on Fumi's GWAS analysis
GWAS_SIG_SNPS=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/data/gwas_sig_snp_names_uniq.txt
#   Need sorted_all_9k_masked_90idt.vcf because 157 out of 158 significant SNPs exist
#       in this VCF file while only 95 of the significant SNPs exist in the
#       OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf file
9K_VCF=/home/morrellp/liux1299/GitHub/9k_BOPA_SNP/BOPA_9k_vcf_Morex_refv1/sorted_all_9k_masked_90idt.vcf
#   VCF file from dataset we are interested in (i.e. OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf)
MAIN_VCF=
#   window size (bp) upstream/downstream of SNP for extract_BED.R
BP=50000
#   Minor Allele Frequency threshold to use for VCF to Htable conversion
MAF=0.015
N_INDIVIDUALS=
P_MISSING=0.015
PREFIX=land_bNAM
OUT_DIR=


function extractSNPs() {
    local snp=$1
    local 9k_vcf=$2
    local prefix=$3
    local out_dir=$4
    #   Create vcf header for significant SNPs
    grep "#" ${9k_vcf} > ${out_dir}/prefix_${snp}_9k_masked_90idt.vcf

    #   Extract significant SNP from 9k masked VCF file
    grep -f "${snp}" ${9k_vcf} >> ${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf
}

export -f extractSNPs

function extractWin() {
    local snp=$1
    local extract_bed=$2
    local bp=$3
    local ss_vcf=$4
    local main_vcf=$5
    local prefix=$6
    local out_dir=$7
    #   Create BED file of n bp upstream/downstream of the significant SNP
    ${extract_bed} ${ss_vcf} ${bp} ${out_dir}/prefix_${snp}_9k_masked_90idt_${bp}win.bed

    #   Extract SNPs that fall within intervals in BED file
    vcfintersect -b ${out_dir}/prefix_${snp}_9k_masked_90idt_${bp}win.bed ${main_vcf} > ${out_dir}/${prefix}_${snp}_intersect.vcf
}

export -f extractWin

function vcfToHtable() {
    local snp=$1
    local vcf_to_htable=$2
    local maf=$3
    local transpose_data=$4
    local prefix=$5
    local out_dir=$6
    #   Convert intersect VCF to fake Hudson table format
    #   This script filters on MAF specified in user provided argument
    #   Output Htable should have marker names (i.e. 11_20909) as columns and
    #       sample names (i.e. WBDC-025) as row names
    ${vcf_to_htable} ${out_dir}/${prefix}_${snp}_intersect.vcf ${maf} > ${out_dir}/tmp_${snp}_intersect_Htable.txt

    #   Sort individuals (i.e. WBDC-025) before transposing data
    (head -n 1 ${out_dir}/tmp_${snp}_intersect_Htable.txt && tail -n +2 ${out_dir}/tmp_${snp}_intersect_Htable.txt | sort -uV -k1,1) > ${out_dir}/${prefix}_${snp}_intersect_Htable_sorted.txt

    #   Transpose data for downstream LD analysis
    ${transpose_data} ${out_dir}/${prefix}_${snp}_intersect_Htable_sorted.txt ${out_dir}

    #   Cleanup temporary files
    rm ${out_dir}/tmp_${snp}_intersect_Htable.txt
}

export -f vcfToHtable

function makeSnpBac() {
    local snp=$1
    local prefix=$2
    local out_dir=$3
    #   Create a SNP_BAC.txt file for all chromosomes
    #   This does not include headers
    #   Output file columns are in the following order: Chr, Physical Position, Marker ID
    awk '{ print $3 "\t" $2 "\t" $1 }' ${out_dir}/${prefix}_${snp}_intersect.vcf | tail -n +2 | sort -V -k2n,2 > ${out_dir}/tmp_snp_bac_all_chr_${prefix}_${snp}.txt

    #   Create tab delimited headers for chromosomes 1-7
    for i in $(seq 1 7)
    do
        printf 'Query_SNP\tPhysPos\tChr\n' > ${out_dir}/SNP_BAC_${prefix}_${snp}-Chr$i.txt
        grep "chr$i" ${out_dir}/tmp_snp_bac_all_chr_${prefix}_${snp}.txt >> ${out_dir}/SNP_BAC_${prefix}_${snp}-Chr$i.txt
    done

    #   Cleanup temporary files
    rm ${out_dir}/tmp_snp_bac_all_chr_${prefix}_${snp}.txt
}

export -f makeSnpBac

function ldDataPrep() {
    local snp=$1
    local ld_data_prep=$2
    local extraction_snps=$3
    local trans_htable=$4
    local prefix=$5
    local out_dir=$6
    #   Run LD_data_prep.sh on whole chromosome including SNP names along heatmap plot
    #   Caveats: Query_SNP must be first column because script sorts by first column
    for i in $(seq 1 7)
    do
        ${ld_data_prep} ${out_dir}/SNP_BAC_${prefix}_${snp}-Chr$i.txt ${trans_htable} Chr$i_${snp} ${out_dir} ${extraction_snps}
    done
}

export -f ldDataPrep

function ldHeatmap() {
    local snp=$1
    local ld_heatmap=$2
    local n_individuals=$3
    local p_missing=$4
    local prefix=$5
    local out_dir=$6
    #   SNP_BAC.txt file must be sorted by SNP names
    #   Genotype data (i.e. *EXISTS.txt genotype data) must be sorted by SNP names
    for i in $(seq 1 7)
    do
        ${ld_heatmap} ${out_dir}/Chr$i_${snp}_sorted_EXISTS.txt ${out_dir}/SNP_BAC_${prefix}_${snp}-Chr$i.txt "Chr$i ${snp}" Chr$i_${snp} ${out_dir} include ${n_individuals} ${p_missing}
    done
}

export -f ldHeatmap


#   Save the filepaths to scripts that we need
cd "${SCRIPT_DIR}"
#   extract_BED.R script creates a BED file that is 50Kb upstream/downstream of
#       the significant SNP
EXTRACT_BED=$(find $(pwd) -name extract_BED.R)
#   VCF_To_Htable-TK.py script reads in a VCF file and outputs a fake Hudson table
#       This script filters sites based on Minor Allele Frequency (MAF) threshold
VCF_TO_HTABLE=$(find $(pwd) -name VCF_To_Htable-TK.py)
#   transpose_data.R script transposes Htable created from VCF_To_Htable-TK.py
TRANSPOSE_DATA=$(find $(pwd) -name transpose_data.R)
#   LD_data_prep.sh script prepares genotyping data for LDheatmap.R script
#       and prevents errors that occur due to samples/markers mismatch between
#       SNP_BAC.txt file and genotype matrix
LD_DATA_PREP=$(find $(pwd) -name LD_data_prep.sh)
EXTRACTION_SNPS=$(find $(pwd) -name extraction_SNPs.pl)
#   LDheatmap.R script generates r2 and D' heatmap plots
LD_HEATMAP=$(find $(pwd) -name LDheatmap.R)

#   Build our SNP list but skip 1st header line
SNP_LIST=($(cat "${GWAS_SIG_SNPS}" | sort -uV))
GSS_LEN=${#SNP_LIST[@]}

#   Still have to figure out how to pass array to these functions
extractSNPs "${GSS_LEN}" "${9K_VCF}" "${PREFIX}" "${OUT_DIR}"
extractWin "${GSS_LEN}" "${EXTRACT_BED}" "${BP}" "${OUT_DIR}"/"${PREFIX}"_${GSS_LEN}_9k_masked_90idt.vcf
vcfToHtable "${GSS_LEN}" "${VCF_TO_HTABLE}" "${MAF}" "${TRANSPOSE_DATA}" "${PREFIX}" "${OUT_DIR}"
makeSnpBac "${GSS_LEN}" "${PREFIX}" "${OUT_DIR}"
ldDataPrep "${GSS_LEN" "${LD_DATA_PREP}" "${EXTRACTION_SNPS}" "${OUT_DIR}"/"${PREFIX}"_${GSS_LEN}_intersect_Htable_sorted_transposed.txt "${PREFIX}" "${OUT_DIR}"
ldHeatmap "${GSS_LEN}" "${LD_HEATMAP}" "${N_INDIVIDUALS}" "${P_MISSING}" "${PREFIX}" "${OUT_DIR}"

# for (( i=0; i<${GSS_LEN}; i++))
# do
#     echo ${SNP_LIST[$i]}
#     #   Create vcf file header for significant SNP

# done
