#!/usr/bin/env Rscript

#   This script takes LD matrix output and filtered physical position data (snp_bac.txt) output from LD_Analysis.sh script
#   and creates r2 decay plots with marks for genes that fall within window.
#   Chaochih Liu - Jan 3, 2018

#   Required arguments:

#   Usage:
#       ./LD_decay_plot.R [prefix] [r2_matrix_dir] [phys_pos_dir] [gene_intervals_dir] [window_size] [output_directory]
#   Where:
#       1) [prefix] is the prefix of our existing HM_r2.txt file (including underscores or hyphens)
#           i.e. filename Chr1-7_SCRI_RS_152696_HM_r2.txt would have prefix Chr1-7_ and SNP name SCRI_RS_152696
#       2) [r2_matrix_dir] is the full filepath to the directory containing r2 matrix output files from LD_Analysis.sh
#           These files should end in _HM_r2.txt
#       3) [phys_pos_dir] is the full filepath to the directory containing snp_bac.txt files from LD_Analysis.sh
#           The filtered snp_bac.txt files are located in ld_data_prep directory outputted from LD_Analysis.sh
#           and should have a naming scheme similar to: SNP_BAC_Chr1-7_11_10143_filtered.txt
#       4) [gene_intervals_dir] is the full filepath to the directory containing gene intervals
#           These should be tab-delimited files.
#       5) [window_size] is the total size of our window (i.e. input 100000 for 50Kb upstream and 50Kb downstream)
#       6) [output_directory] is the full filepath to where we want our plots to go
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

#   Read in gene intervals
readGeneInt <- function(filename) {
    df <- read.table(
        file = filename,
        header = TRUE,
        sep = "\t",
        row.names = NULL
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

calcGeneInterDist <- function(ldData, geneInt.df, t.snp) {
    #   Extract physical position of target SNP
    t.snp.phys <- ldData[ldData$SNPname == t.snp, ]$PhysPos

    #   Add new column of interdist start and end positions
    geneInt.df["InterDist.start"] <- -(t.snp.phys - geneInt.df$gene_start)
    geneInt.df["InterDist.end"] <- -(t.snp.phys - geneInt.df$gene_end)
    return(geneInt.df)
}

#   ld.type is either r2 or Dprime depending on which measure you are plotting
plotLDdecay <- function(ldData, geneInt.df, t.snp, ld.type, ylabel, windowSize, outputDir) {
    #   Set our x axis limits
    winStart <- -(windowSize/2)/1000
    winEnd <- (windowSize/2)/1000
    #   Make our LD decay plot
    pdf(file = paste0(outputDir, "/", t.snp, "_", ld.type, "_LD_decay_with_genes.pdf"))
    plot(
        #   Skip first row because it is self comparison and is filled with NA value
        x = tail(ldData$InterDist/1000, -1),
        y = tail(ldData$r2, -1),
        xlim = c(winStart, winEnd),
        ylim = c(-0.1, 1.0),
        xaxt = "n",
        yaxt = "n",
        bty = "n",
        cex = 0.8,
        xlab = "Physical Distance (Kb)",
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
    axis(side = 1, at = seq(from = winStart, to = winEnd, by = 10), tick = TRUE, pos = 0)
    axis(side = 2, at = seq(from = 0, to = 1.0, by = 0.2), tick = TRUE)
    segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, lty = 3, lwd = 1.5)
    #   Add rectangle for every gene in data frame
    for (i in 1:length(geneInt.df$InterDist.start)) {
        rect(
            xleft = geneInt.df$InterDist.start[i]/1000,
            xright = geneInt.df$InterDist.end[i]/1000,
            ybottom = -0.18,
            ytop = -0.1,
            density = NA, # NA is to suppress shading so we can fill rectangle with color
            col = adjustcolor("dodgerblue", alpha.f = 0.5),
            border = NA)
    }
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
    geneInt.dir <- args[4]
    winSize <- as.numeric(args[5]) # what is the total size of our window? (i.e input 100000 for 50Kb upstream and 50Kb downstream)
    outDir <- args[6]

    #   Store a list of filepaths
    r2.fp <- list.files(path = r2.dir, pattern = "HM_r2.txt", full.names = TRUE)
    phys.fp <- list.files(path = physPos.dir, pattern = "filtered.txt", full.names = TRUE)
    geneInt.fp <- list.files(path = geneInt.dir, pattern = "gene_intervals.txt", full.names = TRUE)

    #   Function that runs all functions for every sample in list
    runAll <- function(ldmatrix.fp, geneIntervals.fp, physPos.fp, file.prefix, window, out.directory) {
        #   Read in LD matrix and physical positions
        tmp.r2.df <- readMatrix(filename = ldmatrix.fp)
        physPos.df <- readPhysPos(filename = physPos.fp)
        geneInterval.df <- readGeneInt(filename = geneIntervals.fp)

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
        #   Calculate distances between gene interval and target SNP
        #   This ensures we have the correct positions when plotting genes
        geneInterDist.df <- calcGeneInterDist(ldData = merged.df, geneInt.df = geneInterval.df, t.snp = targetSNP)

        #   Plot LD decay and save to out directory
        plotLDdecay(
            ldData = interDist.df,
            geneInt.df = geneInterDist.df,
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
        geneInt.fp,
        phys.fp,
        MoreArgs = list(file.prefix = prefix, window = winSize, out.directory = outDir)
    )
}

main() # run program
