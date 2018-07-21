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
#       2) [r2_matrix_list] is a SORTED list of full filepaths to the directory containing r2 matrix output files from LD_Analysis.sh
#           These files should end in _HM_r2.txt
#       3) [phys_pos_list] is a SORTED list of full filepaths to the directory containing snp_bac.txt files from LD_Analysis.sh
#           The filtered snp_bac.txt files are located in ld_data_prep directory outputted from LD_Analysis.sh
#           and should have a naming scheme similar to: SNP_BAC_Chr1-7_11_10143_filtered.txt
#       4) [gene_intervals_list] is a SORTED list of full filepaths to the directory containing gene intervals
#           These should be tab-delimited files.
#       5) [window_size] is the total size of our window (i.e. input 100000 for 50Kb upstream and 50Kb downstream)
#       6) [output_directory] is the full filepath to where we want our plots to go
#           NOTE: filepath to directory should not have last "/" as this will mess up building the filepath

#   Required packages
library(data.table)

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

#   Read in MAF file containing all MAF values for VCF file used in LD analyses
readMAF <- function(filename) {
    df <- fread(
        input = filename,
        header = TRUE,
        sep = "\t"
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

calcInterDist <- function(ldData, t.snp, physPosData) {
    #   Extract physical position for target SNP from LD data
    #   This is the SNP that existed in VCF file and was used for LD calculation
    #   NOT the actual Query SNP (as included in the filename of the LD matrix)
    # tsnpName <- as.character(unique(ldData$targetSnp))
    # tsnpPos <- physPosData[physPosData$SNP == tsnpName, 2]
    tsnpPos <- as.numeric(unique(ldData$tSnpPhysPos))
    #   Calculate distances between target SNP and all other SNPs
    #   Add new column of distances in bp
    ldData["InterDist"] <- -(tsnpPos - ldData$cSnpPhysPos)
    #   Include column of original Query SNP (included in filename of LD matrix)
    ldData["OriginalSNP"] <- t.snp
    return(ldData)
}

calcGeneInterDist <- function(physPosData, geneInt.df, t.snp) {
    #   Extract physical position of target SNP from physical positions data
    t.snp.phys <- physPosData[physPosData$SNP == t.snp, 2]

    #   Add new column of interdist start and end positions
    geneInt.df["InterDist.start"] <- -(t.snp.phys - geneInt.df$gene_start)
    geneInt.df["InterDist.end"] <- -(t.snp.phys - geneInt.df$gene_end)
    return(geneInt.df)
}

#   ld.type is either r2 or Dprime depending on which measure you are plotting
plotLDdecay <- function(ldData, geneInt.df, physPosData, t.snp, original.t.snp, ld.type, ylabel, windowSize, outputDir) {
    #   Set our x axis limits
    winStart <- -(windowSize/2)/1000
    winEnd <- (windowSize/2)/1000
    #   Make our LD decay plot
    pdf(
        file = paste0(outputDir, "/", original.t.snp, "_", ld.type, "_LD_decay_with_genes.pdf"),
        width = 8
    )
    #   Margin of form: c(bottom, left, top, right)
    par(mar = c(5, 5, 5, 2)) # default is c(5, 4, 4, 2)
    plot(
        #   Skip first row because it is self comparison and is filled with NA value
        x = tail(ldData$InterDist/1000, -1),
        y = tail(ldData$r2, -1),
        xlim = c(winStart, winEnd),
        ylim = c(-0.1, 1.0),
        xaxt = "n",
        yaxt = "n",
        bty = "n",
        cex = 0.9,
        lwd = 1.5,
        col = "gray25",
        cex.main = 1.7,
        cex.lab = 1.6,
        xlab = "Physical Distance (Kb)",
        ylab = ylabel,
        main = paste(
            "Original significant SNP:",
            original.t.snp,
            "\nLD decay for SNPs around SNP: ",
            t.snp,
            "\n(", unique(physPosData$Chr), ":",
            as.character(unique(ldData$tSnpPhysPos)),
            "bp )")
    )
    axis(side = 1, at = seq(from = winStart, to = winEnd, by = 20), tick = TRUE, pos = 0, cex.axis = 1.2)
    axis(side = 2, at = seq(from = 0, to = 1.0, by = 0.2), tick = TRUE, cex.axis = 1.2)
    segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, lty = 3, lwd = 1.5)
    if (length(geneInt.df$InterDist.start) >= 1) {
        #   Add rectangle for every gene in data frame
        for (i in 1:length(geneInt.df$InterDist.start)) {
            rect(
                xleft = as.numeric(geneInt.df$InterDist.start[i])/1000,
                xright = as.numeric(geneInt.df$InterDist.end[i])/1000,
                ybottom = -0.18,
                ytop = -0.1,
                density = NA, # NA is to suppress shading so we can fill rectangle with color
                col = adjustcolor("dodgerblue", alpha.f = 0.5),
                border = NA)
        }
    } else {
        cat("No genes found for SNP ")
        cat(t.snp)
        cat("\nNo genes added to plot.")
    }
    #   Turn off graphics
    dev.off()
}

main <- function() {
    #   Take command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided command line arguments
    prefix <- args[1] # what is the prefix of our existing HM_r2.txt and HM_Dprime.txt files?
    r2.list <- args[2]
    physPos.list <- args[3]
    geneInt.list <- args[4]
    maf.fp <- args[5]
    trans.maf.fp <- args[6]
    winSize <- as.numeric(args[7]) # what is the total size of our window? (i.e input 100000 for 50Kb upstream and 50Kb downstream)
    outDir <- args[8]

    #   Store a list of filepaths
    r2.fp <- scan(file = r2.list, what = "", sep = "\n")
    phys.fp <- scan(file = physPos.list, what = "", sep = "\n")
    geneInt.fp <-  scan(file = geneInt.list, what = "", sep = "\n")

    #   Read in data, this is place outside of runAll function to reduce number of times
    #   files need to be read in
    maf <- readMAF(filename = maf.fp)
    trans.maf <- readMAF(filename = trans.maf.fp)

    #   Function that runs all functions for every sample in list
    runAll <- function(ldmatrix.fp, geneIntervals.fp, physPos.fp, maf.df, trans.maf.df, file.prefix, window, out.directory) {
        #   Read in LD matrix and physical positions
        tmp.r2.df <- readMatrix(filename = ldmatrix.fp)
        physPos.df <- readPhysPos(filename = physPos.fp)
        geneInterval.df <- readGeneInt(filename = geneIntervals.fp)

        #   Extract target SNP from filepath of SNP
        #   This works with output files from the LD_analysis.sh script written specifically
        #       for the environmental associations project.
        targetSNP <- extractTargetSNP(filename = ldmatrix.fp, p = file.prefix)
        original.targetSNP <- targetSNP # Save for use in plot title later

        #   Create subset maf.df containing only SNPs in matrix we are currently working with
        sub.maf.df <- maf.df[maf.df$ID %in% rownames(x = tmp.r2.df), ]

        if (targetSNP %in% rownames(tmp.r2.df)) {
            #   Identify index of row that contains targetSNP
            target.index <- which(
                x = rownames(tmp.r2.df) == targetSNP, # extract row containing target SNP
                arr.ind = TRUE # return array index
            )
        } else if (is.na(as.numeric(maf.df[maf.df$ID == targetSNP, 7])) && is.na(as.numeric(trans.maf.df[trans.maf.df$SNP == targetSNP, 2]))) {
            #   This edge case is for when MAF doesn't exist in either MAF file provided.
            cat("Target SNP:", targetSNP, "doesn't exist in either VCF file or 9k genotyping data for 803 landrace accessions. Therefore, no MAF or physical position exists for this SNP. We will save this SNP in it's own file and exit script.\n")
            #   Set target.index to be NA if this case is met
            target.index <- NA
        } else if (is.na(as.numeric(maf.df[maf.df$ID == targetSNP, 7]))) {
            cat("Target SNP: ")
            cat(targetSNP)
            cat(" doesn't exist in VCF file, can't pick closest MAF SNP from 62 NAM Parents VCF.\n
                Will pick closest MAF SNP based on significant SNP MAF in 9k genotyping data for
                803 landrace accessions.\n")
            #   Extract MAF for significant SNP from 803 landrace significant SNPs data frame
            #   containing all MAF values for all SNPs. This is because, some edge cases the target
            #   SNP doesn't exist in the r2 data frame that's the whole reason this section is needed
            tsnp.maf <- as.numeric(trans.maf.df[trans.maf.df$SNP == targetSNP, 2])
            #   Find index of closest MAF SNP
            closest.MAF.index <- which.min(abs(as.numeric(sub.maf.df$MAF) - tsnp.maf))
            #   Extract SNP name of closest MAF SNP
            targetSNP <- sub.maf.df[closest.MAF.index]$ID
            #   Find our "fake" target SNP
            target.index <- which(
                x = rownames(tmp.r2.df) == targetSNP, # extract row containing fake target SNP
                arr.ind = TRUE # return array index
            )
        } else {
            cat("Target SNP:", targetSNP, "isn't in LD analysis matrix, will pick closest MAF SNP from 62 NAM Parents VCF to use as fake target SNP.\n")
            #   Extract MAF for significant SNP from master df containing all MAF values for all SNPs
            #   This is because, some edge cases the target SNP doesn't exist in the r2 data frame
            #   that's the whole reason this section is needed
            tsnp.maf <- as.numeric(maf.df[maf.df$ID == targetSNP, 7])
            #   Find index of closest MAF SNP
            closest.MAF.index <- which.min(abs(as.numeric(sub.maf.df$MAF) - tsnp.maf))
            #   Extract SNP name of closest MAF SNP
            targetSNP <- sub.maf.df[closest.MAF.index]$ID
            #   Find our "fake" target SNP
            target.index <- which(
                x = rownames(tmp.r2.df) == targetSNP, # extract row containing fake target SNP
                arr.ind = TRUE # return array index
            )
        }
        
        #   If target SNP does not exist in either VCF file or 9k genotyping data...
        if (is.na(target.index)) {
            #   save snp to file
            write.table(
                x = targetSNP,
                file = paste0(out.directory, "/", "not_in_VCFor9k_", targetSNP, ".txt"),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE
            )
        #   else, make LD decay plot
        } else {
            #   Reformat LD matrix for compatibility with downstream functions
            r2.df <- r2.reformat(r2.df = tmp.r2.df, index = target.index)
            
            #   Create empty list to store data frames that get reformatted
            d.list <- list()
            #   Loop over every row up until row containing pairwise comparison between
            #       target SNP and all other SNPs
            for (i in 1:target.index) {
                tmp.df <- r2.reformat(r2.df = tmp.r2.df, index = i)
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
            
            #   Calculate distances between SNPs
            interDist.df <- calcInterDist(ldData = merged.df, t.snp = original.targetSNP, physPosData = physPos.df)
            #   Calculate distances between gene interval and target SNP
            #   This ensures we have the correct positions when plotting genes
            geneInterDist.df <- calcGeneInterDist(physPosData = physPos.df, geneInt.df = geneInterval.df, t.snp = targetSNP)
            
            #   Plot LD decay and save to out directory
            plotLDdecay(
                ldData = interDist.df,
                geneInt.df = geneInterDist.df,
                physPosData = physPos.df,
                t.snp = targetSNP,
                original.t.snp = original.targetSNP,
                ld.type = "r2",
                ylabel = expression(paste("LD estimate (", italic("r"^"2"), ")")),
                windowSize = window,
                outputDir = out.directory
            )
        }
    }

    #   Run all functions on list of files
    #   mapply allows me to iterate over two lists of filepaths: r2.fp and phys.fp
    mapply(
        FUN = runAll,
        r2.fp,
        geneInt.fp,
        phys.fp,
        MoreArgs = list(maf.df = maf, trans.maf.df = trans.maf, file.prefix = prefix, window = winSize, out.directory = outDir)
    )
}

main() # run program
