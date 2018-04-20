#!/usr/bin/env Rscript

#   This script takes an LD matrix and reformats it into long format.
#   Chaochih Liu - April 19, 2018

#   Usage:
#       ./r2_matrix_reformat.R [prefix] [r2_matrix_filepath] [phys_pos_filepath] [output_filename]
#   Where:
#   1. [prefix] is the prefix of our existing HM_r2.txt file (including underscores or hyphens)
#       i.e. filename Chr1-7_SCRI_RS_152696_HM_r2.txt would have prefix Chr1-7_ and SNP name
#   2. [r2_matrix_filepath] is the full filepath to the r2 matrix output file from LD_Analysis.sh
#   3. [phys_pos_filepath] is the full filepath to the snp_bac.txt output file from LD_Analysis.sh
#   4. [output_filename] is the full filepath (including filename and .txt file extension) to output file
#   Note: output file is tab-delimited .txt file

readMatrix <- function(filename) {
    df <- read.delim(
        file = filename,
        header = TRUE,
        sep = "\t",
        row.names = 1 # uses first column as row names
    )
    return(df)
}

#   Read in file containing SNP names and physical position info
#   Input file should already be sorted by Query_SNP
readPhysPos <- function(filename) {
    tmp.df <- read.delim(
        file = filename,
        header = TRUE,
        fill = TRUE,
        na.strings = "NA"
    )
    df <- data.frame(
        SNP = tmp.df$Query_SNP,
        PhysPos = tmp.df$PhysPos,
        Chr = tmp.df$Chr
    )
    return(df)
}

#   Extract target SNP name from full filepath
extractTargetSNP <- function(filename, p) {
    #   Given full filepath to matrix, strip path and extract only filename
    f <- basename(filename)
    #   Remove prefix from filename
    no.prefix <- sub(
        pattern = p,
        x = f,
        replacement = "",
        ignore.case = TRUE
    )
    #   Remove suffix and file extension from filename
    targetSNP <- sub(
        pattern = "_HM_r2.txt",
        x = no.prefix,
        replacement = "",
        ignore.case = TRUE
    )
    return(targetSNP)
}

r2.reformat <- function(r2.df, index) {
    #   Extract row of index we are currently working with
    current.row <- r2.df[index, ]
    
    #   Create reformmated data frame
    current.df <- data.frame(
        m_row_snps = rownames(current.row),
        #   BOPA SNPs (i.e. 11_10085), these will have "X" in front of name
        #   We want to remove the "X" because it messes up our file merge later on
        #   If an "X" is found in the column names, substitute all "X" with nothing
        m_col_snps = gsub(pattern = "X", x = colnames(current.row), replacement = ""),
        r2 = as.numeric(current.row[1, ], fill = TRUE)
    )
    return(current.df)
}

mergeFile <- function(ldData, physPosData) {
    #   Merge LD analysis data and physical positions based on matches found between SNPname and SNP columns
    m <- merge(
        x = ldData,
        y = physPosData,
        by.x = "comparisonSnp", # data with LD r2 or D' values
        by.y = "SNP", # physical position data,
        all.x = TRUE, # Rows that do not have a match will remain in dataframe
        all.y = FALSE # Rows that do not have a match with x will not be added to x
    )
    colnames(m)[4] <- "cSnpPhysPos"
    
    #   Extract physical position for target SNP from LD data
    #   This is the SNP that existed in VCF file and was used for LD calculation
    #   NOT the actual Query SNP (as included in the filename of the LD matrix)
    tsnpName <- as.character(unique(ldData$targetSnp))
    tsnpPos <- physPosData[physPosData$SNP == tsnpName, 2]
    #   Add column of target SNP physical positions
    m["tSnpPhysPos"] <- tsnpPos
    return(m)
}

writeOutFile <- function(df, outname) {
    write.table(
        x = df,
        file = outname,
        quote = FALSE,
        sep = "\t",
        eol = "\n",
        col.names = TRUE,
        row.names = FALSE
    )
}

main <- function() {
    #   Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided command line arguments
    prefix <- args[1] # File prefix
    r2.fp <- args[2] # Full filepath to LD r2 matrix
    physPos.fp <- args[3] # Full filepath to phys positions file
    out.fp <- args[4] # Full filepath to output file, include output filename
    
    #   Read in data
    r2 <- readMatrix(filename = r2.fp)
    physPos.df <- readPhysPos(filename = physPos.fp)
    
    #   Extract target SNP from filepath of SNP
    #   This works with output files from the LD_analysis.sh script written specifically
    #       for the environmental associations project.
    targetSNP <- extractTargetSNP(filename = r2.fp, p = prefix)
    
    #   Identify index of row that contains targetSNP
    target.index <- which(
        x = rownames(r2) == targetSNP, # extract row containing target SNP
        arr.ind = TRUE # return array index
    )
    
    #   Create empty list to store data frames that get reformatted
    d.list <- list()
    #   Loop over every row up until row containing pairwise comparison between
    #       target SNP and all other SNPs
    for (i in 1:target.index) {
        tmp.df <- r2.reformat(r2.df = r2, index = i)
        d.list[[i]] <- tmp.df
    }
    #   Combine data frames stored as lists into one data frame
    comb.data <- do.call(what = rbind, args = d.list) # Combined data up until target index number
    
    #   Extract all rows where row.snp column contains target SNP
    rs.tsnp <- comb.data[comb.data$m_row_snps == targetSNP, ]
    #   Rename columns of rs.tsnp
    colnames(rs.tsnp) <- c("targetSnp", "comparisonSnp", "r2")
    cs.tsnp <- comb.data[comb.data$m_col_snps == targetSNP, ]
    #   Swap m_row_snps and m_col_snps in cs.tsnp data frame to make it possible
    #   to merge by physical position in later steps
    cs.swap.tsnp <- data.frame(targetSnp = cs.tsnp$m_col_snps, comparisonSnp = cs.tsnp$m_row_snps, r2 = cs.tsnp$r2)
    #   Add check/print statement to make sure targetSNP matches targetSNP in cs.swap.tsnp
    
    #   Combine these two subsets
    int.df <- rbind(cs.swap.tsnp, rs.tsnp)
    #   Remove rows with NA in r2 column
    clean.data <- int.df[!is.na(int.df$r2), ]
    
    #   Merge LD matrix and physical positions based on matching SNP names
    merged.df <- mergeFile(ldData = clean.data, physPosData = physPos.df)
    
    #   Save to output file
    writeOutFile(df = merged.df, outname = out.fp)
}

main()
