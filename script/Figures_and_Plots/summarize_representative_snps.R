#!/usr/bin/env Rscript

# This script generates a supplemental figure summarizing the distribution of MAF and physical position
# differences between target SNP and representative SNPs. These are GWAS significant SNPs and Fst outlier
# SNPs used in the LD analysis performed for the environmental associations study.
# Chaochih Liu - Feb 6, 2019

library(data.table)
library(ggplot2)
library(gridExtra)

readVcf <- function(filename) {
    df <- fread(
        input = filename,
        header = TRUE,
        skip = "#CHROM",
        sep = "\t"
    )
    return(df)
}

# Read in tab delimited file, including MAF file
readFile <- function(filename) {
    df <- fread(
        input = filename,
        header = TRUE,
        sep = "\t"
    )
    return(df)
}

associateSnpInfo <- function(rep_snps, main_maf, trans_myGD_maf, main_vcf, vcf_9k) {
    # Extract all cases where we are using representative SNP
    # This is when the original SNP and representative SNP differ
    rep_only <- rep_snps[!(rep_snps$originalSNP %in% rep_snps$representativeSNP), ]
    
    # Create placeholder columns for MAF values
    rep_only$originalMaf <- "NA"
    rep_only$representativeMaf <- "NA"
    rep_only$originalPos <- "NA"
    rep_only$representativePos <- "NA"
    
    # Fill in Maf and Pos
    for (i in 1:nrow(rep_only)) {
        # Fill in original SNP MAF column
        if (rep_only$originalSNP[i] %in% main_maf$ID) {
            # Extract MAF value and replace NA placeholder
            rep_only[i]$originalMaf <- as.character(
                main_maf[main_maf$ID == rep_only$originalSNP[i], ]$MAF
            )
        } else {
            # Extract MAF value from trans_myGD_maf and replace NA placeholder
            rep_only[i]$originalMaf <- as.character(
                trans_myGD_maf[trans_myGD_maf$SNP == rep_only$originalSNP[i], ]$Maf
            )
        }
        
        # Fill in representative SNP MAF column
        if (rep_only$representativeSNP[i] %in% main_maf$ID) {
            # Extract MAF value and replace NA placeholder
            rep_only[i]$representativeMaf <- as.character(
                main_maf[main_maf$ID == rep_only$representativeSNP[i], ]$MAF
            )
        } else {
            # Extract MAF value from trans_myGD_maf and replace NA placeholder
            rep_only[i]$representativeMaf <- as.character(
                trans_myGD_maf[trans_myGD_maf$SNP == rep_only$representativeSNP[i], ]$Maf
            )
        }
        
        # Fill in physical positions in original SNP column
        if (rep_only$originalSNP[i] %in% main_vcf$ID) {
            # Add physical position to data frame
            rep_only[i]$originalPos <- as.character(
                main_vcf[main_vcf$ID == rep_only$originalSNP[i], ]$POS
            )
        } else {
            # Extract position from vcf_9k
            rep_only[i]$originalPos <- as.character(
                vcf_9k[vcf_9k$ID == rep_only$originalSNP[i], ]$POS
            )
        }
        
        # Fill in physical positions in representative SNP column
        if (rep_only$representativeSNP[i] %in% main_vcf$ID) {
            # Add physical position to data frame
            rep_only[i]$representativePos <- as.character(
                main_vcf[main_vcf$ID == rep_only$representativeSNP[i], ]$POS
            )
        } else {
            # Extract position from vcf_9k
            rep_only[i]$representativePos <- as.character(
                vcf_9k[vcf_9k$ID == rep_only$representativeSNP[i], ]$POS
            )
        }
    }
    return(rep_only)
}

calcDiff <- function(rep_snps_df) {
    # Calculate differences in MAF values b/w original SNP and representative SNP
    rep_snps_df$mafDiff <- abs(
        as.numeric(rep_snps_df$originalMaf) - as.numeric(rep_snps_df$representativeMaf)
    )
    # Calculate differences in physical position b/w oritinal SNP and representative SNP
    rep_snps_df$posDiff <- abs(
        as.numeric(rep_snps_df$originalPos) - as.numeric(rep_snps_df$representativePos)
    )
    return(rep_snps_df)
}

main <- function() {
    # User input arguments
    rep_snps_fp <- "/Users/chaochih/Dropbox/temp/0.50_ld_analysis/data/all_gwas_representativeSNP_names.txt"
    rep_snps_fst_fp <- "/Users/chaochih/Dropbox/temp/0.50_ld_analysis/data/all_fst_outlier_representativeSNP_names.txt"
    main_vcf_fp <- "/Users/chaochih/Dropbox/temp/0.50_ld_analysis/data/NAM_landraces_100_DP2_renamed_2018-05-21.vcf"
    vcf_9k_fp <- "/Users/chaochih/GitHub/9k_BOPA_SNP/BOPA_9k_vcf_Morex_refv1/sorted_all_9k_masked_90idt.vcf"
    main_maf_fp <- "/Users/chaochih/Dropbox/temp/0.50_ld_analysis/data/NAM_landraces_100_DP2_renamed_2018-05-21.maf"
    trans_myGD_maf_fp <- "/Users/chaochih/Dropbox/temp/0.50_ld_analysis/data/trans_myGD.v2.sorted.maf"
    
    # Read in data
    rep_snps_df <- readFile(filename = rep_snps_fp)
    rep_snps_fst_df <- readFile(filename = rep_snps_fst_fp)
    main_vcf_df <- readVcf(filename = main_vcf_fp)
    vcf_9k_df <- readVcf(filename = vcf_9k_fp)
    main_maf_df <- readFile(filename = main_maf_fp)
    trans_myGD_maf_df <- readFile(filename = trans_myGD_maf_fp)
    
    # Associate MAF and position values to GWAS sig SNPs
    r <- associateSnpInfo(
        rep_snps = rep_snps_df,
        main_maf = main_maf_df,
        trans_myGD_maf = trans_myGD_maf_df,
        main_vcf = main_vcf_df,
        vcf_9k = vcf_9k_df
    )
    rd <- calcDiff(rep_snps_df = r)
    
    # Associate MAF and position values to Fst outlier SNPs
    rfst <- associateSnpInfo(
        rep_snps = rep_snps_fst_df,
        main_maf = main_maf_df,
        trans_myGD_maf = trans_myGD_maf_df,
        main_vcf = main_vcf_df,
        vcf_9k = vcf_9k_df
    )
    rdfst <- calcDiff(rep_snps_df = rfst)
    
    # Summaries of GWAS sig SNPs and Fst outlier SNPs
    # GWAS sig SNPs and Fst outlier SNPs combined MAF
    summary(c(rd$mafDiff, rdfst$mafDiff))
    sd(c(rd$mafDiff, rdfst$mafDiff))
    # GWAS sig SNPs and Fst outlier SNPs combined physical position
    summary(c(rd$posDiff, rdfst$posDiff))
    sd(c(rd$posDiff, rdfst$posDiff))
    
    # Plotting
    # Convert MAF diff values into long format for grouped violin plot
    mplot_df <- data.frame(mafDiff = c(rd$mafDiff, rdfst$mafDiff),
                           groupLab = rep("combinedMafDiff", length(c(rd$mafDiff, rdfst$mafDiff))))
    tmp_gwas_maf_df <- data.frame(mafDiff = rd$mafDiff, groupLab = rep("gwasMafDiff", length(rd$mafDiff)))
    tmp_fst_maf_df <- data.frame(mafDiff = rdfst$mafDiff, groupLab = rep("fstMafDiff", length(rdfst$mafDiff)))
    mplot_df <- rbind(mplot_df, tmp_gwas_maf_df, tmp_fst_maf_df)
    # Generate plots for GWAS and Fst outlier SNPs MAF differences
    m <- ggplot(mplot_df, aes(x = groupLab, y = mafDiff)) + geom_boxplot(width = 0.1) +
        theme_classic(base_size = 14) + 
        scale_x_discrete(labels = c(expression(paste("GWAS + ", italic('F' ['ST']))), "GWAS", expression(italic('F' ['ST'])))) +
        xlab("Representative SNPs") + ylab("MAF Differences") + ggtitle("Differences in MAF between SNPs") +
        theme(plot.title = element_text(hjust = 0.5)) + ylim(0, 0.4)
    
    # Convert position diff values into long format for grouped boxplot
    bplot_df <- data.frame(PosDiff = c(rd$posDiff, rdfst$posDiff),
                           groupLab = rep("combinedPosDiff", length(c(rd$posDiff, rdfst$posDiff))))
    tmp_gwas_df <- data.frame(PosDiff = rd$posDiff, groupLab = rep("gwasPosDiff", length(rd$posDiff)))
    tmp_fst_df <- data.frame(PosDiff = rdfst$posDiff, groupLab = rep("fstPosDiff", length(rdfst$posDiff)))
    bplot_df <- rbind(bplot_df, tmp_gwas_df, tmp_fst_df)
    # Generate plots for GWAS and Fst outlier SNPs physical position differences
    b <- ggplot(bplot_df, aes(x = groupLab, y = PosDiff)) + geom_violin() + geom_boxplot(width = 0.1) +
        theme_classic(base_size = 14) + 
        scale_x_discrete(labels = c(expression(paste("GWAS + ", italic('F' ['ST']))), "GWAS", expression(italic('F' ['ST'])))) +
        xlab("Representative SNPs") + ylab("Physical Position (bp)") + ggtitle("Distance between SNPs") +
        theme(plot.title = element_text(hjust = 0.5))
    # Arrange into 4 panels
    grid.arrange(m, b, nrow = 1)
}

main() # Run the program
