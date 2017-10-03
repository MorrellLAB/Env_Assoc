#!/usr/bin/env Rscript

#   [INSERT brief description of script purpose]
#   Chaochih Liu - September 22, 2017

#   Required arguments:

#   Usage:
#       ./LD_decay_plot.R [prefix] [r2_matrix_dir] [phys_pos_dir] [window_size] [output_directory]
#   Where:
#       1) [prefix] is the prefix of our existing HM_r2.txt file (including underscores or hyphens)
#           i.e. filename Chr1-7_SCRI_RS_152696_HM_r2.txt would have prefix Chr1-7_ and SNP name SCRI_RS_152696
#       2) [r2_matrix_dir] is the full filepath to the directory containing r2 matrix output files from LD_Analysis.sh
#       3) [phys_pos_dir] is the full filepath to the directory containing snp_bac.txt files from LD_Analysis.sh
#       4) [window_size] is the total size of our window (i.e. input 100000 for 50Kb upstream and 50Kb downstream)
#       5) [output_directory] is the full filepath to where we want our plots to go
#           NOTE: filepath to directory should not have last "/" as this will mess up building the filepath


#   Function to read matrix outputted from LDheatmap R package
#   The first row and first column contains SNP names used in LD calculation
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
    #   If any of the SNP names match our target SNP
    if(rownames(df) == snp.name) {
        #   Extract row with target SNP
        target.row <- df[rownames(df) == snp.name, ]
    } else {
        #   Extract 1st row with SNP used for LD calculation
        target.row <- df[1, ]
    }
    
    #   BOPA SNPs (i.e. 11_10085), these will have "X" in front of name
    #   We want to remove the "X" because it messes up our file merge later on
    #   If an "X" is found in the column names
    if(grepl(pattern = "X", x = colnames(target.row))) {
        #   Substitute all "X" with nothing
        s <- gsub(pattern = "X", x = colnames(target.row), replacement = "")
    } else {
        #   Otherwise, store column names of the target row
        s <- colnames(target.row)    
    }
    
    #   Create reformatted data frame
    df <- data.frame(
        targetSNP = rownames(target.row),
        SNPname = s,
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
    #   This is the SNP that existed in VCF file and was used for LD calculation
    #   NOT the actual Query SNP (as included in the filename of the LD matrix)
    tsnpPos <- ldData[ldData$SNPname == as.character(unique(ldData$targetSNP)), 4]
    #   Calculate distances between target SNP and all other SNPs
    #   Add new column of distances in bp
    ldData["InterDist"] <- -(tsnpPos - ldData$PhysPos)
    #   Include column of original Query SNP (included in filename of LD matrix)
    ldData["OriginalSNP"] <- t.snp
    return(ldData)
}

#   ld.type is either r2 or Dprime depending on which measure you are plotting
plotLDdecay <- function(ldData, t.snp, ld.type, ylabel, windowSize, outputDir) {
    #   Set our x axis limits
    winStart <- -(windowSize/2)
    winEnd <- windowSize/2
    #   Make our LD decay plot
    pdf(file = paste0(outputDir, "/", t.snp, "_", ld.type, "_LD_decay.pdf"))
    plot(
        #   Skip first row because it is self comparison and is filled with NA value
        x = tail(ldData$InterDist, -1),
        y = tail(ldData$r2, -1),
        xlim = c(winStart, winEnd),
        ylim = c(0, 1.0),
        xaxt = "n",
        cex = 0.8,
        xlab = "Physical Distance (bp)",
        ylab = ylabel,
        main = paste(
            "LD decay for SNPs around SNP: ",
            t.snp,
            "\nExisting SNP in VCF used for calculation: ",
            as.character(unique(ldData$targetSNP)),
            "\nWindow Center (",
            ldData[ldData$SNPname == as.character(unique(ldData$targetSNP)), 4],
            "bp)")
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
    r2.dir <- args[2]
    physPos.dir <- args[3]
    winSize <- as.numeric(args[4]) # what is the total size of our window? (i.e input 100000 for 50Kb upstream and 50Kb downstream)
    outDir <- args[5]
    
    #   Store a list of filepaths
    r2.fp <- list.files(path = r2.dir, pattern = "HM_r2.txt", full.names = TRUE)
    phys.fp <- list.files(path = physPos.dir, pattern = "filtered.txt", full.names = TRUE)
    
    #   Function that runs all functions for every sample in list
    runAll <- function(ldmatrix.fp, physPos.fp, file.prefix, window, out.directory) {
        #   Read in LD matrix and physical positions
        tmp.r2.df <- readMatrix(filename = ldmatrix.fp)
        physPos.df <- readPhysPos(filename = physPos.fp)
        
        #   Extract target SNP from filepath of SNP
        #   This works with output files from the LD_analysis.sh script written specifically
        #       for the environmental associations project.
        targetSNP <- extractTargetSNP(filename = ldmatrix.fp, p = file.prefix)
        
        #   Reformat LD matrix for compatibility with downstream functions
        r2.df <- r2.reformat(df = tmp.r2.df, snp.name = targetSNP)
        
        #   Merge LD matrix and physical positions based on matching SNP names
        merged.df <- mergeFile(ldData = r2.df, physPosData = physPos.df)
        
        #   Calculate distances between SNPs
        interDist.df <- calcInterDist(ldData = merged.df, t.snp = targetSNP)
        
        #   Plot LD decay and save to out directory
        plotLDdecay(
            ldData = interDist.df,
            t.snp = targetSNP,
            ld.type = "r2",
            ylabel = expression(paste("LD estimate (", "r"^"2", ")")),
            windowSize = window,
            outputDir = out.directory
        )
    }
    
    #   Run all functions on list of files
    #   mapply allows me to iterate over two lists of filepaths: r2.fp and phys.fp
    mapply(
        FUN = runAll,
        r2.fp,
        phys.fp,
        MoreArgs = list(file.prefix = prefix, window = winSize, out.directory = outDir)
    )
}

main() # run program
