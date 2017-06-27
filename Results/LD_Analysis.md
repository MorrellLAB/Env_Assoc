# LD Analysis - Landrance Environmental Association

All scripts pertaining to this analysis are located in [`Insert GitHub directory`]().

Generate LD heatmaps using R package `LDheatmap`.

Note: `LDheatmap` package does not include self-comparison values in it's matrix. Self-comparison values are replaced with an "NA" automatically.

---

## Data & Data Preparation

### Data

The VCF file used to create the genotype matrix is: `/home/morrellp/liux1299/Shared/Projects/Barley_NAM_Parents/SNP_calling/Variants/Barley_NAM_Parents_Final_Fixed.vcf`

This file has already been filtered for missing data, multi-nucleotide polymorphisms, and non-biallelic SNPs.

Check total number of markers before starting analysis:

```bash
grep -v "#" Barley_NAM_Parents_Final_Fixed.vcf| cut -f 1,2,3 | wc -l

```

### Data Preparation

#### Step 1: Create BED file of 100Kb windows around significant SNPs

Create a list of significant SNPs:

```bash
cut -d "," -f 1 compiled.5e_4.0.01.v2_physPos.csv > ~/Downloads/significant_snp_names.txt
grep "#" sorted_all_9k_masked_90idt.vcf > ~/Downloads/env_assoc_sig_snps_9k.vcf
grep -f significant_snp_names.txt ~/GitHub/9k_BOPA_SNP/BOPA_9k_vcf_Morex_refv1/sorted_all_9k_masked_90idt.vcf >> ~/Downloads/env_assoc_sig_snps_9k.vcf
```

We are looking at 50Kb upstream and downstream of the significant SNP. Used `extract_BED.R` script to create bed file.

```bash
#   Extract SNPs that fall within intevals in BED file
vcfintersect -b env_assoc_sig_snps_9k.bed ~/Shared/Projects/Barley_NAM_Parents/SNP_calling/Variants/Barley_NAM_Parents_Final_Fixed.vcf > env_assoc_sig_snps_NAM.vcf
```

#### Step 2: Create a fake Hudson table from VCF file

Convert VCF to fake Hudson table format using a modified version of Tom's `VCF_To_Htable-TK.py` script located in [`Barley_Inversions/analyses/LD_Analysis` GitHub](https://github.com/MorrellLAB/Barley_Inversions/tree/master/analyses/LD_Analysis). This script filters on **MAF of 0.05**. The output Htable should have marker names (i.e. S1H1_42239) as columns and sample names (i.e. WBDC-025) as row names.

Command used to convert file format:

```bash
#   This script requires Python 2
VCF_To_Htable-TK.py WBDC_July2016_production_PTP_filt.vcf.recode.vcf > WBDC_July2016_production_PTP_Htable.txt
#   Then sort individuals (i.e. WBDC-XXX) before transposing data frame
(head -n 1 WBDC_July2016_production_PTP_Htable.txt && tail -n +2 WBDC_July2016_production_PTP_Htable.txt | sort -u -k1,1) > WBDC_July2016_production_PTP_Htable_sorted.txt
```

The total number of loci after running the `VCF_To_Htable-TK.py` script is: **39,167**

#### Step 4: Transpose data frame

Next, I transposed the `WBDC_July2016_production_Htable.txt` dataframe using `transpose_data-inv.R` script located in [`Barley_Inversions/analyses/LD_Analysis` GitHub](https://github.com/MorrellLAB/Barley_Inversions/tree/master/analyses/LD_Analysis).

```bash
transpose_data-inv.R WBDC_July2016_production_PTP_Htable_sorted.txt
```

#### Step 5: Create SNP_BAC.txt file

To use the same LD analysis scripts, I will need to create a fake SNP_BAC.txt file that has 2 columns used in the `LDheatmap.R` script:
- `Query_SNP`
- `PhysPos`

Commands used to create fake SNP_BAC.txt file:

```bash
#   Create tab delimited headers for chromosomes 1-7
for i in $(seq 1 7); do printf 'Query_SNP\tPhysPos\tChr\n' > SNP_BAC_WBDC_July2016-Chr$i.txt; done

#   Create a SNP_BAC.txt file for all chromosomes
#   This does not include headers
#   Output file columns are in following order
#       Chr, Physical Position, Marker ID
(grep -v "##" WBDC_July2016_production_PTP.vcf | awk '{ print $3 "\t" $2 "\t" $1 }' | tail -n +2 | sort -k2n,2) >> SNP_BAC_WBDC_July2016_production_all.txt
```

---

### Analyses Overview

#### Step 1: LD data prep

Run `LD_data_prep.sh` on whole chromosomes first (using *EXISTS.txt genotype data)

Caveats: `Query_SNP` must be first column because script sorts by first column.

```bash
LD_data_prep.sh SNP_BAC_WBDC_July2016-Chr2.txt WBDC_July2016_production_PTP_Htable_sorted_transposed.txt Chr2 ~/Shared/Projects/WBDC_inversions/LDheatmap/WBDC_GBS/data

#   Repeated above command for Chr 2, 4, and 5
#   Number of markers that exist in Chr 2
tail -n +2 SNP_BAC_Chr2_filtered.txt | wc -l
6342
#   Number of markers that exist in Chr 4
tail -n +2 SNP_BAC_Chr4_filtered.txt | wc -l
4858
#   Number of markers that exist in Chr 5
tail -n +2 SNP_BAC_Chr5_filtered.txt | wc -l
5645
```

#### Step 2: LD heatmap script

Run `LDheatmap.R` or `LD_translocation.R` scripts to generate heatmap.

Caveats:
- SNP_BAC.txt file must be sorted by SNP names
- Genotype data frame (i.e. *EXISTS.txt genotype data) must be sorted by SNP names

---

### Results