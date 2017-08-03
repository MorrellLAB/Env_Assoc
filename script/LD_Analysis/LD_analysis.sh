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
module load parallel

#   Usage message:
#   This script performs LD analysis on 100Kb windows around significant SNPs from Environmental
#   Associations project GWAS analysis. This script pulls together all parts of the analysis
#   performed by multiple scripts.
#   Please fill out user provided argument fields and submit script as job.


#   User provided arguments
#   Where is the directory containing this script?
#   Make sure all scripts needed are in same directory as this script
SCRIPT_DIR=/home/morrellp/liux1299/GitHub/Env_Assoc/script/LD_Analysis
#   List of significant SNP names based on Fumi's GWAS analysis
GWAS_SIG_SNPS=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/data/gwas_sig_snp_names_uniq.txt
#   Need sorted_all_9k_masked_90idt.vcf because 157 out of 158 significant SNPs exist
#       in this VCF file while only 95 of the significant SNPs exist in the
#       OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf file
VCF_9K=/home/morrellp/liux1299/GitHub/9k_BOPA_SNP/BOPA_9k_vcf_Morex_refv1/sorted_all_9k_masked_90idt.vcf
#   VCF file from dataset we are interested in (i.e. OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf)
MAIN_VCF=/home/morrellp/liux1299/Shared/Projects/Barley_NAM_Parents/SNP_calling/Variants/New_Filtering/OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf
#   window size (bp) upstream/downstream of SNP for extract_BED.R
BP=50000
#   Minor Allele Frequency threshold to use for VCF to Htable conversion
MAF=0.01
#   Missing data threshold to use for filtering
P_MISSING=0.15
#   What prefix do we want to use for our output files?
PREFIX=ld_Barley_NAM
#   Where is our output directory?
OUT_DIR=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps

#   Extract GWAS significant SNPs from 9k_masked_90idt.vcf
function extractSNPs() {
    local snp=$1
    local vcf_9k=$2
    local prefix=$3
    local out_dir=$4
    #   Create vcf header for significant SNPs
    grep "#" ${vcf_9k} > ${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf

    #   Extract significant SNP from 9k masked VCF file
    grep ${snp} ${vcf_9k} >> ${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf
}

export -f extractSNPs

#   Extract all SNPs that fall within window size defined
function extractWin() {
    local snp=$1
    local extract_bed=$2
    local bp=$3
    local ss_vcf=$4
    local main_vcf=$5
    local prefix=$6
    local out_dir=$7
    #   Create BED file of n bp upstream/downstream of the significant SNP
    ${extract_bed} ${ss_vcf} ${bp} ${out_dir}/${prefix}_${snp}_9k_masked_90idt_${bp}win.bed

    #   Extract SNPs that fall within intervals in BED file
    vcfintersect -b ${out_dir}/${prefix}_${snp}_9k_masked_90idt_${bp}win.bed ${main_vcf} > ${out_dir}/${prefix}_${snp}_intersect.vcf
}

export -f extractWin

#   Use VCF_To_Htable-TK.py to create fake Hudson table
#   from intersect.vcf file
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

    #   Remove "X" in marker names
    sed 's/X//g' ${out_dir}/${prefix}_${snp}_intersect_Htable_sorted_transposed.txt > ${out_dir}/${prefix}_${snp}_intersect_Htable_sorted_transposed_noX.txt

    #   Cleanup temporary files
    rm ${out_dir}/tmp_${snp}_intersect_Htable.txt
}

export -f vcfToHtable

#   Create a SNP_BAC.txt file that will be used in
#   LD_data_prep.sh and LDheatmap.R script
function makeSnpBac() {
    local snp=$1
    local prefix=$2
    local out_dir=$3
    #   Create tab delimited header for all chr
    printf 'Query_SNP\tPhysPos\tChr\n' > ${out_dir}/SNP_BAC_${prefix}_${snp}-all_chr.txt
    #   Create a SNP_BAC.txt file for all chromosomes
    #   This does not include headers
    #   Output file columns are in the following order: Chr, Physical Position, Marker ID
    awk '{ print $3 "\t" $2 "\t" $1 }' ${out_dir}/${prefix}_${snp}_intersect.vcf | tail -n +2 | sort -V -k2n,2 >> ${out_dir}/SNP_BAC_${prefix}_${snp}-all_chr.txt
}

export -f makeSnpBac

#   Prepare genotype data for LD heatmap:
#   pull out SNPs that exist, filter out SNPs that don't exist, and sort
#   SNP_BAC.txt data and genotype data
function ldDataPrep() {
    local snp=$1
    local ld_data_prep=$2
    local extraction_snps=$3
    local trans_htable=$4
    local prefix=$5
    local out_dir=$6
    #   Run LD_data_prep.sh on whole chromosome including SNP names along heatmap plot
    #   Caveats: Query_SNP must be first column because script sorts by first column
    ${ld_data_prep} ${out_dir}/SNP_BAC_${prefix}_${snp}-all_chr.txt ${trans_htable} Chr1-7_${snp} ${out_dir} ${extraction_snps}
}

export -f ldDataPrep

#   Perform LD analysis and generate LD heatmap plots
function ldHeatMap() {
    local snp=$1
    local ld_heatmap=$2
    local n_individuals=$3
    local p_missing=$4
    local prefix=$5
    local out_dir=$6
    #   SNP_BAC.txt file must be sorted by SNP names
    #   Genotype data (i.e. *EXISTS.txt genotype data) must be sorted by SNP names
    ${ld_heatmap} ${out_dir}/Chr1-7_${snp}_sorted_EXISTS.txt ${out_dir}/SNP_BAC_Chr1-7_${snp}_filtered.txt "Chr1-7 ${snp}" Chr1-7_${snp} ${out_dir} exclude ${n_individuals} ${p_missing}
}

export -f ldHeatMap


#   Number of Individuals we have data for (i.e. WBDC)
N_INDIVIDUALS=$(grep "#CHROM" "${MAIN_VCF}" | tr '\t' '\n' | tail -n +10 | wc -l)
echo "Number of individuals in data:"
echo ${N_INDIVIDUALS}
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
#   Number of GWAS Significant SNPs (GSS) in array
GSS_LEN=${#SNP_LIST[@]}
echo "Number of GWAS Significant SNPs in array:"
echo ${GSS_LEN}

#   Run program for each significant SNP in parallel
echo "Extracting significant SNPs from 9k_masked_90idt.vcf file..."
parallel -v extractSNPs {} "${VCF_9K}" "${PREFIX}" "${OUT_DIR}" ::: "${SNP_LIST[@]}"
echo "Done extracting significant SNPs."

echo "Extracting all SNPs that fall within window defined..."
parallel -v extractWin {} "${EXTRACT_BED}" "${BP}" "${OUT_DIR}"/"${PREFIX}"_{}_9k_masked_90idt.vcf "${MAIN_VCF}" "${PREFIX}" "${OUT_DIR}" ::: "${SNP_LIST[@]}"
echo "Done extracting SNPs within window."

echo "Converting VCF to fake Hudson table..."
parallel -v vcfToHtable {} "${VCF_TO_HTABLE}" "${MAF}" "${TRANSPOSE_DATA}" "${PREFIX}" "${OUT_DIR}" ::: "${SNP_LIST[@]}"
echo "Done converting VCF to fake Hudson table."

echo "Creating SNP_BAC.txt file..."
parallel -v makeSnpBac {} "${PREFIX}" "${OUT_DIR}" ::: "${SNP_LIST[@]}"
echo "Done creating SNP_BAC.txt."

echo "Preparing data for LD analysis..."
parallel -v ldDataPrep {} "${LD_DATA_PREP}" "${EXTRACTION_SNPS}" "${OUT_DIR}"/"${PREFIX}"_{}_intersect_Htable_sorted_transposed_noX.txt "${PREFIX}" "${OUT_DIR}" ::: "${SNP_LIST[@]}"
echo "Done preparing data."

echo "Running LD analysis..."
parallel -v ldHeatMap {} "${LD_HEATMAP}" "${N_INDIVIDUALS}" "${P_MISSING}" "${PREFIX}" "${OUT_DIR}" ::: "${SNP_LIST[@]}"
echo "Done."
