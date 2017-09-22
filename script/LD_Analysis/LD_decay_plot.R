#!/usr/bin/env Rscript

#   [INSERT brief description of script purpose]
#   Chaochih Liu - September 22, 2017

#   Required arguments:

#   Usage:

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

r2.reformat <- function(df, snp.name) {
    #   Extract row with query SNP compared to all other SNPs
    target.row <- df[rownames(df) == snp.name, ]
    
    df <- data.frame(
        targetSNP = rownames(target.row),
        SNPname = gsub(
            pattern = "X", # replace X in SNP name, R headers have "X" added to beginning if they start with a number
            x = colnames(target.row),
            replacement = ""
        ),
        r2 = as.numeric(target.row[1, ]) # r2 values
    )
    return(df)
}

mergeFile <- function(ldData, physPosData) {
    #   Merge LD analysis data and physical positions based on matches found between SNPname and SNP columns
    m <- merge(
        x = ldData,
        y = physPosData,
        by.x = "SNPname", # data with LD r2 or D' values
        by.y = "SNP", # physical position data,
        all.x = TRUE, # Rows that do not have a match will remain in dataframe
        all.y = FALSE # Rows that do not have a match with x will not be added to x
    )
    return(m)
}

calcInterDist <- function(ldData, t.snp) {
    #   Extract physical position for target SNP from LD data
    tsnpPos <- ldData[ldData$SNPname == t.snp, 4]
    #   Calculate distances between target SNP and all other SNPs
    #   Add new column of distances
    merged.df["InterDist"] <- -(tsnpPos - merged.df$PhysPos)
}

plotLDdecay <- function(ldData, t.snp) {
    # superscript 
    plot(
        x = merged.df$InterDist,
        y = merged.df$r2,
        cex = 0.75,
        xlab = "Distance (bp)",
        ylab = expression(paste("LD estimate (", "r"^"2", ")"))
    )
    abline(v = 0)
}

main <- function() {
    #   Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    prefix <- "Chr1-7_" # what is the prefix of our existing HM_r2.txt and HM_Dprime.txt files?
    r2matrix <- "/Users/chaochih/Downloads/r2_decay_test/ld_results/Chr1-7_11_10380_HM_r2.txt"
    snpbac <- "/Users/chaochih/Downloads/r2_decay_test/ld_data_prep/SNP_BAC_Chr1-7_11_10380_filtered.txt"
    
    tmp <- readMatrix(filename = r2matrix)
    physPos <- readPhysPos(filename = snpbac)
    targetSNP <- extractTargetSNP(filename = r2matrix, p = prefix)
    r2.df <- r2.reformat(df = tmp, snp.name = targetSNP)
    merged.df <- mergeFile(ldData = r2.df, physPosData = physPos)
    calcInterDist(ldData = r2.df, t.snp = targetSNP)
    
    #   Following is all test code, will remove once script is working
    # f <- basename(r2matrix)
    # no.prefix <- sub(pattern = prefix, x = f, replacement = "", ignore.case = TRUE)
    # snp.name <- sub(pattern = "_HM_r2.txt", x = no.prefix, replacement = "", ignore.case = TRUE)
    #   Extract row that is the query SNP compared to all other SNPs
    #query.row <- tmp[rownames(tmp) == snp.name, ]
    # test <- data.frame(targetSNP = rownames(query.row), SNPname = colnames(query.row), r2 = as.numeric(query.row[1, ]))
    #   Extract physical position for target SNP from LD data
    # tsnpPos <- merged.df[merged.df$SNPname == targetSNP, 4]
    merged.df["InterDist"] <- -(tsnpPos - merged.df$PhysPos)
}
