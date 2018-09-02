#!/bin/bash

set -e
set -o pipefail

#   This script contains only the functions for the LD GWAS significant SNPs analysis

function makeOutDirs() {
    local out_dir=$1
    #   Check if out directory exists, if not make it
    mkdir -p "${out_dir}/extracted_sig_snps_vcf" \
             "${out_dir}/temp" \
             "${out_dir}/extracted_window" \
             "${out_dir}/Htable" \
             "${out_dir}/snp_bac" \
             "${out_dir}/ld_data_prep" \
             "${out_dir}/ld_results/ldheatmap_error_snps"
}

export -f makeOutDirs

#   Extract GWAS significant SNPs from 9k_masked_90idt.vcf
function extractSNPs() {
    local snp=$1
    local vcf_9k=$2
    local prefix=$3
    local out_dir=$4
    #   Create vcf header for significant SNPs
    grep "#" "${vcf_9k}" > "${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf"

    #   Extract significant SNP from 9k masked VCF file
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

#   Driver function
function main() {
    local script_dir=$1
    local gwas_sig_snps=$2
    local vcf_9k=$3
    local main_vcf=$4
    local bp=$5
    local maf=$6
    local p_missing=$7
    local prefix=$8
    local out_dir=$9
    #   Make out directories
    makeOutDirs "${out_dir}"
    #   Number of Individuals we have data for (i.e. WBDC)
    n_individuals=$(grep "#CHROM" "${main_vcf}" | tr '\t' '\n' | tail -n +10 | wc -l)
    echo "Number of individuals in data:"
    echo "${n_individuals}"
    #   Save the filepaths to scripts that we need
    #   extract_BED.R script creates a BED file that is 50Kb upstream/downstream of
    #   the fst outlier SNP
    extract_bed="${script_dir}/extract_BED.R"
    #   VCF_To_Htable-TK.py script reads in a VCF file and outputs a fake Hudson table
    #   This script filters sites based on Minor Allele Frequency (MAF) threshold
    vcf_to_htable="${script_dir}/VCF_To_Htable-TK.py"
    #   transpose_data.R script transposes Htable created from VCF_To_Htable-TK.py
    transpose_data="${script_dir}/transpose_data.R"
    #   LD_data_prep.sh script prepares genotyping data for LDheatmap.R script
    #   and prevents errors that occur due to samples/markers mismatch between
    #   SNP_BAC.txt file and genotype matrix
    ld_data_prep="${script_dir}/LD_data_prep.sh"
    extraction_snps="${script_dir}/extraction_SNPs.pl"
    #   LDheatmap.R script generates r2 and D' heatmap plots
    ld_heatmap="${script_dir}/LDheatmap.R"

    #   Build our SNP list but skip 1st header line
    snp_list=($(cat "${gwas_sig_snps}" | sort -uV))
    #   Number of GWAS Significant SNPs (GSS) in array
    gss_len=${#snp_list[@]}
    echo "Number of GWAS Significant SNPs in array:"
    echo ${gss_len}

    echo "Extracting significant SNPs from 9k_masked_90idt.vcf file..."
    #   Running extractSNPs will output the following file:
    #       1) prefix_9k_masked_90idt.vcf file(s) contains significant SNPs from GWAS analysis
    #       2) sig_snp_not_in9k.txt file contains significant SNPs that don't exist in the sorted_all_9k_masked_90idt.vcf file
    #   We start with a list of GWAS significant SNP names,
    #   pull those SNPs from the sorted_all_9k_masked_90idt.vcf file
    #   to create VCF files containing 1 significant SNP/VCF file
    touch "${out_dir}/extracted_sig_snps_vcf/sig_snp_not_in_9k.txt"
    parallel extractSNPs {} "${vcf_9k}" "${prefix}" "${out_dir}/extracted_sig_snps_vcf" ::: "${snp_list[@]}"
    echo "Done extracting significant SNPs."

    echo "Removing non-existent SNP from bash array..."
    #   Filter out and remove SNPs that don't exist from bash array
    delete=($(cat "${out_dir}/extracted_sig_snps_vcf/sig_snp_not_in_9k.txt"))
    echo ${snp_list[@]} | tr ' ' '\n' > "${out_dir}/temp/tmp_snp_list.txt"
    snp_list_filt=($(grep -vf "${out_dir}/extracted_sig_snps_vcf/sig_snp_not_in_9k.txt" "${out_dir}/temp/tmp_snp_list.txt"))
    rm "${out_dir}/temp/tmp_snp_list.txt"
    echo "Done removing non-existent SNP from bash array."
    echo "Number of GWAS Significant SNPs that exist in 9k_masked_90idt.vcf file:"
    echo ${#snp_list_filt[@]}

    echo "Extracting all SNPs that fall within window defined..."
    #   Running extractWin will output the following files:
    #       1) BED file(s) of n Kb upstream/downstream of SNP (should have 1 line within file)
    #       2) intersect.vcf file(s) that contains all SNPs that fall within BED file interval
    #   We start with our prefix_9k_masked_90idt.vcf files
    #   and create a BED file interval n bp upstream/downstream of significant SNP.
    #   Then we use vcfintersect for BED file we created and our VCF file of
    #   interest (i.e. OnlyLandrace_biallelic_Barley_NAM_Parents_Final_renamed.vcf)
    #   to pull down all SNPs that fall within our BED file interval.
    parallel extractWin {} "${extract_bed}" "${bp}" "${out_dir}/extracted_sig_snps_vcf/${prefix}_{}_9k_masked_90idt.vcf" "${main_vcf}" "${prefix}" "${out_dir}/extracted_window" ::: "${snp_list_filt[@]}"
    echo "Done extracting SNPs within window."

    #   Filter out intersect.vcf files that are empty by removing SNP name from bash array
    intersect_vcf=($(find "${out_dir}"/extracted_window/*.vcf))
    snp_int_vcf=()
    for i in "${intersect_vcf[@]}"
    do
        #   redirect filename into wc to get integer only
        num_lines=$(wc -l < ${i})
        #   If there is only 1 line in the file (the header line),
        #   save the full filepath to file
        if [ "${num_lines}" -eq "1" ]
        then
            basename ${i} >> "${out_dir}/extracted_window/empty_intersect_vcf.txt"
            basename ${i} | sed -e s/^${prefix}_// -e s/_intersect.vcf// >> "${out_dir}/extracted_window/empty_intersect_vcf_SNPnamesOnly.txt"
        else
            #   Extract only the SNP name from filename using sed substitution
            #   to remove prefix and suffix.
            #   This works because ${PREFIX} is defined as variable at top of script
            #   and extractWin function uses ${PREFIX} in output file names.
            #   The suffix of extractWin output files is always "_intersect.vcf"
            snp=$(basename ${i} | sed -e s/^${prefix}_// -e s/_intersect.vcf//)
            #   Add SNPs we want to use in downstream functions to new array
            snp_int_vcf+=(${snp})
        fi
    done

    echo "Converting VCF to fake Hudson table..."
    #   Running vcfToHtable will filter on MAF and output the following files:
    #       1) Htable_sorted.txt file(s) which is the VCF converted to fake Hudson table format
    #       2) Htable_sorted_transposed.txt file(s) which outputs SNPs as rows and individuals as columns
    #       3) Htable_sorted_transposed_noX.txt which removes "X" in marker names
    #   We start with our intersect.vcf file(s) and convert them to a fake Hudson table format
    #   and filter based on MAF (specified above under user provided argument).
    #   The output will have marker names (i.e. 11_20909) as columns and sample naems (i.e. WBDC-025) as row names.
    parallel vcfToHtable {} "${vcf_to_htable}" "${maf}" "${transpose_data}" "${prefix}" "${out_dir}" ::: "${snp_int_vcf[@]}"
    echo "Done converting VCF to fake Hudson table."

    echo "Creating SNP_BAC.txt file..."
    #   Running makeSnpBac will output file(s) that contain 3 columns:
    #       1) Query_SNP which is the SNP name
    #       2) PhysPos which is the physical position
    #       3) Chr which is the chromosome
    #   Output files are sorted by physical position (column 2)
    parallel makeSnpBac {} "${prefix}" "${out_dir}" ::: "${snp_int_vcf[@]}"
    echo "Done creating SNP_BAC.txt."

    echo "Preparing data for LD analysis..."
    #   Running ldDataPrep will output the following files:
    #       1) sorted_EXISTS.txt which contains SNPs that exist in our genotyping data
    #       2) NOT_EXISTS.txt is a list of SNPs that do not exist in our genotyping data but exist in our SNP_BAC.txt file
    #       3) SNP_BAC_filtered.txt has all non-existent SNPs removed so it doesn't cause errors when using LDheatmap command in R
    parallel ldDataPrep {} "${ld_data_prep}" "${extraction_snps}" "${out_dir}"/Htable/"${prefix}"_{}_intersect_Htable_sorted_transposed_noX.txt "${prefix}" "${out_dir}" ::: "${snp_int_vcf[@]}"
    echo "Done preparing data."

    echo "Running LD analysis..."
    #   Running ldHeatMap will output the following files:
    #       1) SNP_info-empty_cols.csv is a list of samples with empty columns
    #       2) SNP_info-failed_snps.csv is a list of incompatible genotype columns
    #       3) SNP_info-missing_data_cols.csv is a list of SNPs that had greater than n% missing data
    #           (missing data threshold is defined under user provided arguments section)
    #       4) compatibleSnps.txt is a matrix of SNPs that will be used for LDheatmap analyses
    #       5) HM_r2.pdf is a heatmap for r2 calculation
    #       6) HM_Dprime.pdf is a heatmap for D' calculation
    #       7) HM_r2.txt is a matrix of r2 values used in heatmap
    #       8) HM_Dprime.txt is a matrix of D' values used in heatmap
    parallel ldHeatMap {} "${ld_heatmap}" "${n_individuals}" "${p_missing}" "${out_dir}" ::: "${snp_int_vcf[@]}"
    echo "Done with all analyses."
}

export -f main
