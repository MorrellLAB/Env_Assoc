#!/bin/bash

set -e
set -o pipefail

#   This script contains only the functions for the LD fst outliers analysis

#   Extract fstOutlier SNPs from 9k_masked_90idt.vcf
function extractSNPs() {
    local snp=$1
    local vcf_9k=$2
    local prefix=$3
    local out_dir=$4
    #   Create vcf header for fst outlier SNPs
    grep "#" "${vcf_9k}" > "${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf"

    #   Extract fst outlier SNP from 9k masked VCF file
    if grep -q "${snp}" "${vcf_9k}"; then
        #   If SNP exists, extract SNP from 9k masked VCF file
        grep "${snp}" "${vcf_9k}" >> "${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf"
    else
        #   If SNP doesn't exist, save SNP in another file
        echo "${snp} does not exist in 9k_masked.vcf file." >&2
        echo "${snp}" >> "${out_dir}/sig_snp_not_in_9k.txt"
        rm "${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf"
    fi
}

export -f extractSNPs

function filterSNPs() {
    local delete_array=$1
    local snp_list_array=$2
    local out_dir=$3
    echo ${snp_list_array[@]} | tr ' ' '\n' > "${out_dir}/temp/tmp_snp_list.txt"
    snp_list_filt=($(grep -vf "${out_dir}/extracted_sig_snps_vcf/sig_snp_not_in_9k.txt" "${out_dir}/temp/tmp_snp_list.txt"))
    #rm "${out_dir}/temp/tmp_snp_list.txt"
    echo "Done removing non-existent SNP from bash array."
    echo "Number of fst outlier SNPs that exist in 9k_masked_90idt.vcf file:"
    echo ${#snp_list_filt[@]}
}

export -f filterSNPs

#   Extract all SNPs that fall within window size defined
function extractWin() {
    local snp=$1
    local extract_bed=$2
    local bp=$3
    local ss_vcf=$4
    local main_vcf=$5
    local prefix=$6
    local out_dir=$7
    #   Create BED file of n bp upstream/downstream of the fst outlier SNP
    "${extract_bed}" "${ss_vcf}" "${bp}" "${out_dir}/${prefix}_${snp}_9k_masked_90idt_${bp}win.bed"

    #   Extract SNPs that fall within intervals in BED file
    vcfintersect -b "${out_dir}/${prefix}_${snp}_9k_masked_90idt_${bp}win.bed" "${main_vcf}" > "${out_dir}/${prefix}_${snp}_intersect.vcf"
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
    python "${vcf_to_htable}" "${out_dir}/extracted_window/${prefix}_${snp}_intersect.vcf" "${maf}" > "${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt"
    #   Sort individuals (i.e. WBDC-025) before transposing data
    (head -n 1 "${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt" && tail -n +2 "${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt" | sort -uV -k1,1) > "${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted.txt"
    #   Transpose data for downstream LD analysis
    "${transpose_data}" "${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted.txt" "${out_dir}/Htable"
    #   Remove "X" in marker names
    sed -e 's/X//g' "${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted_transposed.txt" > "${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted_transposed_noX.txt"
    #   Cleanup temporary files
    rm "${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt"
}

export -f vcfToHtable

#   Create a SNP_BAC.txt file that will be used in
#   LD_data_prep.sh and LDheatmap.R script
function makeSnpBac() {
    local snp=$1
    local prefix=$2
    local out_dir=$3
    #   Create tab delimited header for all chr
    printf 'Query_SNP\tPhysPos\tChr\n' > "${out_dir}/snp_bac/SNP_BAC_${prefix}_${snp}-all_chr.txt"
    #   Create a SNP_BAC.txt file for all chromosomes
    #   This does not include headers
    #   Output file columns are in the following order: Chr, Physical Position, Marker ID
    awk '{ print $3 "\t" $2 "\t" $1 }' "${out_dir}/extracted_window/${prefix}_${snp}_intersect.vcf" | tail -n +2 | sort -V -k2n,2 >> "${out_dir}/snp_bac/SNP_BAC_${prefix}_${snp}-all_chr.txt"
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
    "${ld_data_prep}" "${out_dir}/snp_bac/SNP_BAC_${prefix}_${snp}-all_chr.txt" "${trans_htable}" Chr1-7_"${snp}" "${out_dir}/ld_data_prep" "${extraction_snps}"
}

export -f ldDataPrep

#   Perform LD analysis and generate LD heatmap plots
function ldHeatMap() {
    local snp=$1
    local ld_heatmap=$2
    local n_individuals=$3
    local p_missing=$4
    local out_dir=$5
    #   SNP_BAC.txt file must be sorted by SNP names
    #   Genotype data (i.e. *EXISTS.txt genotype data) must be sorted by SNP names
    #   Arg 1: genotyping data
    #   Arg 2: physical positions
    #   Arg 3: heatmap plot name
    #   Arg 4: output file prefix, no space
    #   Arg 5: out directory
    #   Arg 6: include or exclude SNP names in heatmap
    #   Arg 7: number of individuals
    #   Arg 8: missing data threshold
    "${ld_heatmap}" "${out_dir}/ld_data_prep/Chr1-7_${snp}_sorted_EXISTS.txt" "${out_dir}/ld_data_prep/SNP_BAC_Chr1-7_${snp}_filtered.txt" "Chr1-7 ${snp}" "Chr1-7_${snp}" "${out_dir}/ld_results" "exclude" "${n_individuals}" "${p_missing}"
    #   Move files associated with SNP that had an error (i.e. gdat undefined column error) when running ldHeatMap function to subdirectory
    if [ -f "${out_dir}/ld_results/Chr1-7_${snp}_ldheatmap_fn_error.txt" ]
    then
        echo "Chr1-7_${snp}_ldheatmap_fn_error.txt found. Moving files to ldheatmap_error_snps directory."
        mv "${out_dir}"/ld_results/*"${snp}"* "${out_dir}/ld_results/ldheatmap_error_snps"
    else
        echo "LD heatmap function completed successfully for snp: ${snp}."
    fi
}

export -f ldHeatMap
