#!/usr/bin/env bash

#PBS -l mem=22gb,nodes=1:ppn=16,walltime=10:00:00
#PBS -m abe
#PBS -M liux1299@umn.edu
#PBS -q lab

set -e
set -o pipefail

module load R/3.4.3
module load python2/2.7.8 # Tom's VCF_To_Htable-TK.py script runs in Python 2
module load vcflib_ML/1.0.0
module load parallel

#   Usage message:
#       This script performs LD analysis on 100Kb windows around significant SNPs from Environmental
#       Associations project GWAS analysis. This script pulls together all parts of the analysis
#       performed by multiple scripts.
#       Please fill out user provided argument fields and submit script as job.


#   User provided arguments
#   Where is the directory containing this script?
#   Make sure all scripts needed are in same directory as this script
SCRIPT_DIR=/home/morrellp/liux1299/GitHub/Env_Assoc/script/LD_Analysis
#   List of significant SNP names based on Fumi's GWAS analysis
GWAS_SIG_SNPS=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/data/gwas_sig_snp_names-11_10380.txt
#   Need sorted_all_9k_masked_90idt.vcf because 157 out of 158 significant SNPs exist
#       in this VCF file while only 95 of the significant SNPs exist in the
#       OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf file
VCF_9K=/home/morrellp/liux1299/GitHub/9k_BOPA_SNP/BOPA_9k_vcf_Morex_refv1/sorted_all_9k_masked_90idt.vcf
#   VCF file from dataset we are interested in (i.e. OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf)
MAIN_VCF=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/data/11_10380_Only_landrace_biallelic_NAM_final_renamed_2018-03-25.vcf
#   window size (bp) upstream/downstream of SNP for extract_BED.R
BP=100000
#   Minor Allele Frequency threshold to use for VCF to Htable conversion (i.e. 0.01 for 1% MAF)
MAF=0.01
#   Missing data threshold to use for filtering (i.e. 0.15 for 15% missing data)
P_MISSING=0.15
#   What prefix do we want to use for our output files?
PREFIX=ld_Barley_NAM_200Kb
#   Where is our output directory?
OUT_DIR=/home/morrellp/liux1299/Shared/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/11_10380_only/0.15_missing

#   Extract GWAS significant SNPs from 9k_masked_90idt.vcf
function extractSNPs() {
    local snp=$1
    local vcf_9k=$2
    local prefix=$3
    local out_dir=$4
    #   Create vcf header for significant SNPs
    grep "#" ${vcf_9k} > ${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf

    #   Extract significant SNP from 9k masked VCF file
    if grep -q "${snp}" ${vcf_9k}; then
        #   If SNP exists, extract SNP from 9k masked VCF file
        grep "${snp}" ${vcf_9k} >> ${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf
    else
        #   If SNP doesn't exist, save SNP in another file
        echo "${snp} does not exist in 9k_masked.vcf file." >&2
        echo ${snp} >> ${out_dir}/sig_snp_not_in_9k.txt
        rm ${out_dir}/${prefix}_${snp}_9k_masked_90idt.vcf
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
    ${vcf_to_htable} ${out_dir}/extracted_window/${prefix}_${snp}_intersect.vcf ${maf} > ${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt

    #   Sort individuals (i.e. WBDC-025) before transposing data
    (head -n 1 ${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt && tail -n +2 ${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt | sort -uV -k1,1) > ${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted.txt

    #   Transpose data for downstream LD analysis
    ${transpose_data} ${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted.txt ${out_dir}/Htable

    #   Remove "X" in marker names
    sed 's/X//g' ${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted_transposed.txt > ${out_dir}/Htable/${prefix}_${snp}_intersect_Htable_sorted_transposed_noX.txt

    #   Cleanup temporary files
    rm ${out_dir}/Htable/tmp_${snp}_intersect_Htable.txt
}

export -f vcfToHtable

#   Create a SNP_BAC.txt file that will be used in
#   LD_data_prep.sh and LDheatmap.R script
function makeSnpBac() {
    local snp=$1
    local prefix=$2
    local out_dir=$3
    #   Create tab delimited header for all chr
    printf 'Query_SNP\tPhysPos\tChr\n' > ${out_dir}/snp_bac/SNP_BAC_${prefix}_${snp}-all_chr.txt
    #   Create a SNP_BAC.txt file for all chromosomes
    #   This does not include headers
    #   Output file columns are in the following order: Chr, Physical Position, Marker ID
    awk '{ print $3 "\t" $2 "\t" $1 }' ${out_dir}/extracted_window/${prefix}_${snp}_intersect.vcf | tail -n +2 | sort -V -k2n,2 >> ${out_dir}/snp_bac/SNP_BAC_${prefix}_${snp}-all_chr.txt
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
    ${ld_data_prep} ${out_dir}/snp_bac/SNP_BAC_${prefix}_${snp}-all_chr.txt ${trans_htable} Chr1-7_${snp} ${out_dir}/ld_data_prep ${extraction_snps}
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
    ${ld_heatmap} ${out_dir}/ld_data_prep/Chr1-7_${snp}_sorted_EXISTS.txt ${out_dir}/ld_data_prep/SNP_BAC_Chr1-7_${snp}_filtered.txt "Chr1-7 ${snp}" Chr1-7_${snp} ${out_dir}/ld_results exclude ${n_individuals} ${p_missing}
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

#   Check if out directory exists, if not make it
mkdir -p "${OUT_DIR}" "${OUT_DIR}"

echo "Extracting significant SNPs from 9k_masked_90idt.vcf file..."
#   Check if out directory exists, if not make directory
mkdir -p "${OUT_DIR}/extracted_sig_snps_vcf" "${OUT_DIR}/extracted_sig_snps_vcf"
#   Running extractSNPs will output the following file:
#       1) prefix_9k_masked_90idt.vcf file(s) contains significant SNPs from GWAS analysis
#       2) sig_snp_not_in9k.txt file contains significant SNPs that don't exist in the sorted_all_9k_masked_90idt.vcf file
#   We start with a list of GWAS significant SNP names,
#   pull those SNPs from the sorted_all_9k_masked_90idt.vcf file
#   to create VCF files containing 1 significant SNP/VCF file
touch "${OUT_DIR}"/extracted_sig_snps_vcf/sig_snp_not_in_9k.txt
parallel extractSNPs {} "${VCF_9K}" "${PREFIX}" "${OUT_DIR}"/extracted_sig_snps_vcf ::: "${SNP_LIST[@]}"
echo "Done extracting significant SNPs."


echo "Removing non-existent SNP from bash array..."
#   Check if out directory exists, if not make directory
mkdir -p "${OUT_DIR}/temp" "${OUT_DIR}/temp"
#   Filter out and remove SNPs that don't exist from bash array
DELETE=($(cat "${OUT_DIR}"/extracted_sig_snps_vcf/sig_snp_not_in_9k.txt))
echo ${SNP_LIST[@]} | tr ' ' '\n' > "${OUT_DIR}"/temp/tmp_snp_list.txt
SNP_LIST_FILT=($(grep -vf "${OUT_DIR}"/extracted_sig_snps_vcf/sig_snp_not_in_9k.txt "${OUT_DIR}"/temp/tmp_snp_list.txt))
rm "${OUT_DIR}"/temp/tmp_snp_list.txt
echo "Done removing non-existent SNP from bash array."
echo "Number of GWAS Significant SNPs that exist in 9k_masked_90idt.vcf file:"
echo ${#SNP_LIST_FILT[@]}


echo "Extracting all SNPs that fall within window defined..."
#   Check if out directory exists, if not make directory
mkdir -p "${OUT_DIR}/extracted_window" "${OUT_DIR}/extracted_window"
#   Running extractWin will output the following files:
#       1) BED file(s) of n Kb upstream/downstream of SNP (should have 1 line within file)
#       2) intersect.vcf file(s) that contains all SNPs that fall within BED file interval
#   We start with our prefix_9k_masked_90idt.vcf files
#   and create a BED file interval n bp upstream/downstream of significant SNP.
#   Then we use vcfintersect for BED file we created and our VCF file of
#   interest (i.e. OnlyLandrace_biallelic_Barley_NAM_Parents_Final_renamed.vcf)
#   to pull down all SNPs that fall within our BED file interval.
parallel extractWin {} "${EXTRACT_BED}" "${BP}" "${OUT_DIR}"/extracted_sig_snps_vcf/"${PREFIX}"_{}_9k_masked_90idt.vcf "${MAIN_VCF}" "${PREFIX}" "${OUT_DIR}"/extracted_window ::: "${SNP_LIST_FILT[@]}"
echo "Done extracting SNPs within window."


#   Filter out intersect.vcf files that are empty by removing SNP name from bash array
INTERSECT_VCF=($(find "${OUT_DIR}"/extracted_window/*.vcf))
SNP_INT_VCF=()
for i in "${INTERSECT_VCF[@]}"
do
    #   redirect filename into wc to get integer only
    num_lines=$(wc -l < ${i})
    #   If there is only 1 line in the file (the header line),
    #   save the full filepath to file
    if [ "${num_lines}" -eq "1" ]
    then
        basename ${i} >> "${OUT_DIR}"/extracted_window/empty_intersect_vcf.txt
        basename ${i} | sed -e s/^${PREFIX}_// -e s/_intersect.vcf// >> "${OUT_DIR}"/extracted_window/empty_intersect_vcf_SNPnamesOnly.txt
    else
        #   Extract only the SNP name from filename using sed substitution
        #   to remove prefix and suffix.
        #   This works because ${PREFIX} is defined as variable at top of script
        #   and extractWin function uses ${PREFIX} in output file names.
        #   The suffix of extractWin output files is always "_intersect.vcf"
        SNP=$(basename ${i} | sed -e s/^${PREFIX}_// -e s/_intersect.vcf//)
        #   Add SNPs we want to use in downstream functions to new array
        SNP_INT_VCF+=(${SNP})
    fi
done


echo "Converting VCF to fake Hudson table..."
#   Check if out directory exists, if not make directory
mkdir -p "${OUT_DIR}/Htable" "${OUT_DIR}/Htable"
#   Running vcfToHtable will filter on MAF and output the following files:
#       1) Htable_sorted.txt file(s) which is the VCF converted to fake Hudson table format
#       2) Htable_sorted_transposed.txt file(s) which outputs SNPs as rows and individuals as columns
#       3) Htable_sorted_transposed_noX.txt which removes "X" in marker names
#   We start with our intersect.vcf file(s) and convert them to a fake Hudson table format
#   and filter based on MAF (specified above under user provided argument).
#   The output will have marker names (i.e. 11_20909) as columns and sample naems (i.e. WBDC-025) as row names.
parallel vcfToHtable {} "${VCF_TO_HTABLE}" "${MAF}" "${TRANSPOSE_DATA}" "${PREFIX}" "${OUT_DIR}" ::: "${SNP_INT_VCF[@]}"
echo "Done converting VCF to fake Hudson table."


echo "Creating SNP_BAC.txt file..."
#   Check if out directory exists, if not make directory
mkdir -p "${OUT_DIR}/snp_bac" "${OUT_DIR}/snp_bac"
#   Running makeSnpBac will output file(s) that contain 3 columns:
#       1) Query_SNP which is the SNP name
#       2) PhysPos which is the physical position
#       3) Chr which is the chromosome
#   Output files are sorted by physical position (column 2)
parallel makeSnpBac {} "${PREFIX}" "${OUT_DIR}" ::: "${SNP_INT_VCF[@]}"
echo "Done creating SNP_BAC.txt."


echo "Preparing data for LD analysis..."
#   Check if out directory exists, if not make directory
mkdir -p "${OUT_DIR}/ld_data_prep" "${OUT_DIR}/ld_data_prep"
#   Running ldDataPrep will output the following files:
#       1) sorted_EXISTS.txt which contains SNPs that exist in our genotyping data
#       2) NOT_EXISTS.txt is a list of SNPs that do not exist in our genotyping data but exist in our SNP_BAC.txt file
#       3) SNP_BAC_filtered.txt has all non-existent SNPs removed so it doesn't cause errors when using LDheatmap command in R
parallel ldDataPrep {} "${LD_DATA_PREP}" "${EXTRACTION_SNPS}" "${OUT_DIR}"/Htable/"${PREFIX}"_{}_intersect_Htable_sorted_transposed_noX.txt "${PREFIX}" "${OUT_DIR}" ::: "${SNP_INT_VCF[@]}"
echo "Done preparing data."


echo "Running LD analysis..."
mkdir -p "${OUT_DIR}/ld_results" "${OUT_DIR}/ld_results"
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
parallel ldHeatMap {} "${LD_HEATMAP}" "${N_INDIVIDUALS}" "${P_MISSING}" "${PREFIX}" "${OUT_DIR}" ::: "${SNP_INT_VCF[@]}"
echo "Done."
