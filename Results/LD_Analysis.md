# LD Analysis - Landrance Environmental Association

All scripts pertaining to this analysis are located in [`Insert GitHub directory`]().

Generate LD heatmaps using R package `LDheatmap`.

Note: `LDheatmap` package does not include self-comparison values in it's matrix. Self-comparison values are replaced with an "NA" automatically.

---

## Data & Data Preparation

### Data

The VCF file used to create the genotype matrix is: `/home/morrellp/liux1299/Shared/Projects/Barley_NAM_Parents/SNP_calling/Variants/New_Filtering/OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf`

This file has already been filtered for missing data, multi-nucleotide polymorphisms, and non-biallelic SNPs.

Rename ID column of VCF file so we know which SNPs in the 100Kb window are the significant SNPs. Run `add_SNP_ID_To_VCF.py` script (job script is located in `Env_Assoc/script/LD_Analysis/job_scripts/rename_vcf_id.job`).

Check total number of markers before starting analysis:

```bash
grep -v "#" OnlyLandrace_Barley_NAM_Parents_Final_renamed.vcf | wc -l
645685
```

## LD Analysis

`LD_analysis.sh` will call on all the necessary scripts to perform analysis. Fill out "User provided arguments" variables at the top of the script (lines 26-44) and then submit script as a job.

---

### Results