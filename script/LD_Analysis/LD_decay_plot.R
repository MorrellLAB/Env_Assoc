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
        r2 = as.numeric(target.row[1, ], fill = TRUE) # r2 values
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
    #   Add new column of distances in bp
    ldData["InterDist"] <- -(tsnpPos - ldData$PhysPos)
    return(ldData)
}

plotLDdecay <- function(ldData, t.snp, windowSize, outputDir) {
    #   Set our x axis limits
    winStart <- -(windowSize/2)
    winEnd <- windowSize/2
    #   Make our LD decay plot
    pdf(file = paste0(outputDir, "/", ".pdf"))
    plot(
        #   Skip first row because it is self comparison and is filled with NA value
        x = tail(ldData$InterDist, -1),
        y = tail(ldData$r2, -1),
        xlim = c(winStart, winEnd),
        ylim = c(0, 1.0),
        xaxt = "n",
        cex = 0.8,
        xlab = "Physical Distance (bp)",
        ylab = expression(paste("LD estimate (", "r"^"2", ")")),
        main = paste("LD decay - ", targetSNP)
    )
    axis(side = 1, at = seq(from = winStart, to = winEnd, by = 10000), tick = TRUE)
    abline(v = 0, lty = 3, lwd = 1.5)
    #   Turn off graphics
    dev.off()
}

main <- function() {
    #   Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided command line arguments
    prefix <- args[1] # what is the prefix of our existing HM_r2.txt and HM_Dprime.txt files?
    r2matrix <- args[2]
    snpbac <- args[3]
    winSize <- as.numeric(args[4]) # what is the total size of our window? (i.e input 100000 for 50Kb upstream and 50Kb downstream)
    outDir <- args[5]
    #   arguments used for testing
    prefix <- "Chr1-7_" # what is the prefix of our existing HM_r2.txt and HM_Dprime.txt files?
    r2matrix <- "/Users/chaochih/Downloads/r2_decay_test/ld_results/Chr1-7_11_10085_HM_r2.txt"
    snpbac <- "/Users/chaochih/Downloads/r2_decay_test/ld_data_prep/SNP_BAC_Chr1-7_11_10085_filtered.txt"
    winSize <- 100000
    outDir <- "/Users/chaochih/Downloads/r2_decay_test/test_plots"
    
    tmp <- readMatrix(filename = r2matrix)
    physPos <- readPhysPos(filename = snpbac)
    targetSNP <- extractTargetSNP(filename = r2matrix, p = prefix)
    r2.df <- r2.reformat(df = tmp, snp.name = targetSNP)
    merged.df <- mergeFile(ldData = r2.df, physPosData = physPos)
    intDist.df <- calcInterDist(ldData = merged.df, t.snp = targetSNP)
    plotLDdecay(ldData = intDist.df, t.snp = targetSNP, windowSize = winSize)
    
    #   Following is all test code, will remove once script is working
    # f <- basename(r2matrix)
    # no.prefix <- sub(pattern = prefix, x = f, replacement = "", ignore.case = TRUE)
    # snp.name <- sub(pattern = "_HM_r2.txt", x = no.prefix, replacement = "", ignore.case = TRUE)
    #   Extract row that is the query SNP compared to all other SNPs
    #query.row <- tmp[rownames(tmp) == snp.name, ]
    # test <- data.frame(targetSNP = rownames(query.row), SNPname = colnames(query.row), r2 = as.numeric(query.row[1, ]))
    #   Extract physical position for target SNP from LD data
    # tsnpPos <- merged.df[merged.df$SNPname == targetSNP, 4]
    #merged.df["InterDist"] <- -(tsnpPos - merged.df$PhysPos)
}
