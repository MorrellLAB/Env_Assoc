# LD Analysis

Scripts in this directory are used to run LD analysis on GWAS significant SNPs and Fst outlier SNPs.

---

## Running Analysis

Run the `.job` scripts to reproduce analyses. Scripts with the `.job` file extension are used to submit analyses to supercomputing system. The `.job` files store filepaths and parameters for analysis. The `*_0.40.job`, `*_0.50.job`, etc. indicate different parameters tested.

#### Example

We will use our GWAS significant SNPs LD analysis for 50% missingness threshold as an example. For other parameters and Fst outlier SNPs analysis, we will follow the same steps as below (just the job scripts will have different names).

**Step 1**

For the GWAS significant SNPs analyses, run the following on MSI:

```bash
qsub LD_analysis_gwas_0.50.job
```

Once this `.job` script begins running, it stores all filepaths and parameters in variables. Then it sources the `LD_analysis_gwas.sh` script to pull all the necessary functions for this analysis. Functions within the `LD_analysis_gwas.sh` script require many of the scripts seen in this directory. Here is a list of custom scripts (included in this directory) used in `LD_analysis_gwas.sh` (Note: there are some tools, like R, vcflib, parallel, and Python 2, that need to be installed prior to running this analysis):
- `extract_BED.R`
- `VCF_To_Htable-TK.py`
- `transpose_data.R`
- `LD_data_prep.sh`
- `extraction_SNPs.pl`
- `LDheatmap.R`

**Step 2**

After this analysis is done running, we will use additional scripts to generate the LD decay plots. First, we need to generate gene interval files to feed to our plotting script:

```bash
#   In ~/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50 directory
find $(pwd) -name "*9k_masked_90idt_100000win.bed" | sort > bed_win_list.txt

#   Submit job to MSI to generate gene intervals
qsub gwas_gene_intervals_0.50.job
```

The `gwas_gene_intervals_0.50.job` script sources the `gene_intervals.sh` script, which stores all the necessary functions for this step.

**Step 3**

Calculate minor allele frequencies for our VCF file. This only needs to be done once per VCF file used for the analysis. Basically, if we are testing different parameters for the same VCF file, we still only need to run this step once.

```bash
#   Should take <5 minutes and only need to run one time for one VCF file
#   In directory: ~/Projects/Land_Env_Assoc/Analysis/LD_Analysis/data
python3 /home/morrellp/liux1299/GitHub/Env_Assoc/script/LD_Analysis/vcf_maf_snp_id.py \
    NAM_landraces_100_DP2_renamed_2018-05-21.vcf \
    > NAM_landraces_100_DP2_renamed_2018-05-21.maf
```

**Step 4**

Generate LD decay plots.

```bash
#   LD heatmap matrix with r2 values in ld_results directory
cd ~/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold
#   Check number of significant SNPs used in analyses
find $(pwd) -name "*HM_r2.txt" | sort -V | wc -l
153
#   Check number of significant SNPs that failed during LDheatmap function
find $(pwd) -name "*ldheatmap_fn_error.txt" | sort -V | wc -l
3
#   Generate r2 sample list (and copy files to Dropbox if running locally)
find $(pwd) -name "*HM_r2.txt" | sort -V > r2_file_list.txt

#   Only extract SNPs used in analyses
#   this is important b/c some SNPs were filtered due to missingness, etc.
cp r2_file_list.txt ld_analysis_snp_names.txt
#   Strip prefix and suffix in Vim or using sed -e
vim ld_analysis_snp_names.txt
:%s,/home/morrellp/liux1299/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold/ld_results/Chr1-7_,,g
:%s,_HM_r2.txt,,g

#   From list of SNPs used in analyses, extract associated physPos file
#   In dir: ~/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold/ld_data_prep
find $(pwd) -name "*filtered.txt" | sort -V | grep -f ld_analysis_snp_names.txt | sort -V > phys_pos_file_list.txt

#   Only keep files that exist in ld_results and ld_data_prep
find $(pwd) -name "*gene_intervals.txt" | sort -V | grep -f ld_analysis_snp_names.txt | sort -V > interval_list_filtered.txt

#   Check and make sure SNPs are all sorted the same way across the 3 lists of filepaths:
sed -e 's,/home/morrellp/liux1299/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold/ld_data_prep/SNP_BAC_Chr1-7_,,' -e 's,_filtered.txt,,' phys_pos_file_list.txt > tmp_check_order_physPos.txt
sed -e 's,/home/morrellp/liux1299/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold/ld_results/Chr1-7_,,' -e 's,_HM_r2.txt,,' r2_file_list.txt > tmp_check_order_r2.txt
sed -e 's,/home/morrellp/liux1299/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold/gene_intervals/ld_Barley_NAM_200Kb_,,' -e 's,_gene_intervals.txt,,' interval_list_filtered.txt > tmp_check_order_gene_interval.txt

module load R/3.4.3
#   LD decay plots
~/GitHub/Env_Assoc/script/LD_Analysis/LD_decay_plot_with_genes.R \
    Chr1-7_ \
    r2_file_list.txt \
    phys_pos_file_list.txt \
    interval_list_filtered.txt \
    /home/morrellp/liux1299/Projects/Land_Env_Assoc/Analysis/LD_Analysis/data/NAM_landraces_100_DP2_renamed_2018-05-21.maf \
    /home/morrellp/liux1299/Projects/Land_Env_Assoc/Analysis/LD_Analysis/data/trans_myGD.v2.sorted.maf \
    200000 \
    ~/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold/r2_decay_plots_200Kb_0.50_2018-08-29 \
    &> ~/Projects/Land_Env_Assoc/Analysis/LD_Analysis/results/gwas_sig_snps_200Kb/0.50_threshold/r2_decay_plots_200Kb_0.50.log
```

We are done after this step and can now look at our plots.
