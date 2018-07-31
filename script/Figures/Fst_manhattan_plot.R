#!/usr/bin/env Rscript
#   Chaochih Liu - April 27, 2018

#   This script takes in Fst output files with a column of physical positions and generates
#   a manhattan plot. This script highlights gene hits of interest.

#   Usage:
#   ./Fst_manhattan_plot_noGeneHit.R [fst.txt] [gene_hits.txt] [outlier_threshold] [plot_title] [out_dir] [out_name]

#   Where:
#   1) [fst.txt] contains our Fst results and the following columns:
#       SNP, Chromosome, Cumulative_cM, FST, Chr_2016, and PhysPos_2016
#   2) [gene_hits.txt] is a file with SNPs that are gene hits. The column name must be "SNPs_9k" for script to work
#   3) [outlier_threshold] is our fst outlier threshold (i.e. for 97.5% threshold, use 0.975)
#   4) [plot_title] must be wrapped in double quotes if there are spaces
#   5) [out_dir] is the full filepath to our output directory. NOTE: remove ending "/"
#       otherwise filename will have one too many slashes
#   6) [out_name] is our output file names including file extension ".pdf" (i.e. fst_lat.pdf)

#   Define function to read in .txt
readData <- function(filename) {
    df <- read.table(
        file = filename,
        header = TRUE,
        sep = "\t",
        na.strings = "NA"
    )
    #   Remove rows with missing physical position
    df.filter <- df[!is.na(df$PhysPos_2016), ]
    return(df.filter)
}

readGeneHits <- function(filename) {
    df <- read.table(
        file = filename,
        header = TRUE,
        sep = "\t",
        na.strings = "NA"
    )
    return(df)
}

#   Function to replace "chr" and "H" with nothing
replace <- function(df) {
    df$Chr_2016 <- sub("chr", "", df$Chr_2016)
    df$Chr_2016 <- sub("H", "", df$Chr_2016)
    return(df)
}

chrTicks <- function(chr1.whole, chr2.whole, chr3.whole, chr4.whole, chr5.whole, chr6.whole, chr7.whole) {
    #   Row 1 is the scaled position of the center of the chromosome
    #   Row 2 is the scaled position of the end of the chromosome
    chr1.t <- c(chr1.whole/2, chr1.whole)
    chr2.t <- c(chr2.whole/2 + chr1.whole, chr1.whole + chr2.whole)
    chr3.t <- c(chr3.whole/2 + chr1.whole + chr2.whole, chr1.whole + chr2.whole + chr3.whole)
    chr4.t <- c(chr4.whole/2 + chr1.whole + chr2.whole + chr3.whole, chr4.whole + chr1.whole + chr2.whole + chr3.whole)
    chr5.t <- c(chr5.whole/2 + chr1.whole + chr2.whole + chr3.whole + chr4.whole, chr5.whole + chr1.whole + chr2.whole + chr3.whole + chr4.whole)
    chr6.t <- c(chr6.whole/2 + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole, chr6.whole + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole)
    chr7.t <- c(chr7.whole/2 + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole + chr6.whole, chr7.whole + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole + chr6.whole)
    ticks <- cbind(chr1.t, chr2.t, chr3.t, chr4.t, chr5.t, chr6.t, chr7.t)
    return(ticks)
}

#   Function to scale physical positions for plotting
scaledPos <- function(df, c1.w, c2.w, c3.w, c4.w, c5.w, c6.w, c7.w) {
    #   Subset data
    df.chr1 <- df[df$Chr_2016 == 1, ]
    df.chr2 <- df[df$Chr_2016 == 2, ]
    df.chr3 <- df[df$Chr_2016 == 3, ]
    df.chr4 <- df[df$Chr_2016 == 4, ]
    df.chr5 <- df[df$Chr_2016 == 5, ]
    df.chr6 <- df[df$Chr_2016 == 6, ]
    df.chr7 <- df[df$Chr_2016 == 7, ]
    #   Add new column of scaled positions for plotting
    df.chr1["scaled.BP"] <- df.chr1$PhysPos_2016
    df.chr2["scaled.BP"] <- df.chr2$PhysPos_2016 + c1.w
    df.chr3["scaled.BP"] <- df.chr3$PhysPos_2016 + c1.w + c2.w
    df.chr4["scaled.BP"] <- df.chr4$PhysPos_2016 + c1.w + c2.w + c3.w
    df.chr5["scaled.BP"] <- df.chr5$PhysPos_2016 + c1.w + c2.w + c3.w + c4.w
    df.chr6["scaled.BP"] <- df.chr6$PhysPos_2016 + c1.w + c2.w + c3.w + c4.w + c5.w
    df.chr7["scaled.BP"] <- df.chr7$PhysPos_2016 + c1.w + c2.w + c3.w + c4.w + c5.w + c6.w
    #   Combine data frames
    df.all <- rbind(
        df.chr1,
        df.chr2,
        df.chr3,
        df.chr4,
        df.chr5,
        df.chr6,
        df.chr7
    )
    return(df.all)
}

plot.manhattan <- function(df, fst.outlier.threshold, plot.title, ticks) {
    #   Extract all rows where Chromosome is an odd number
    df.odd <- df[!as.numeric(as.character(df$Chr_2016)) %% 2 == 0, ]
    df.even <- df[as.numeric(as.character(df$Chr_2016)) %% 2 == 0, ]
    
    par(mar = c(5, 5, 5, 2))
    #   All SNPs from odd chromosomes will be black
    plot(
        x = df.odd$scaled.BP,
        y = df.odd$FST,
        col = adjustcolor(col = "gray10", alpha.f = 0.9),
        xlim = c(0, 4569031868),
        ylim = c(0, 0.8),
        pch = 1,
        cex = 0.6,
        cex.lab = 1.4,
        cex.main = 1.6,
        cex.axis = 1.2,
        xlab = "Chromosome",
        ylab = expression('F'[ST]), # subscript
        xaxt = "n",
        main = plot.title
    )
    par(new = TRUE)
    #   All SNPs from even chromosomes will be gray
    plot(
        x = df.even$scaled.BP,
        y = df.even$FST,
        col = adjustcolor(col = "gray70", alpha.f = 0.9),
        xlim = c(0, 4569031868),
        ylim = c(0, 0.8),
        pch = 1,
        cex = 0.6,
        cex.lab = 1.4,
        cex.main = 1.6,
        cex.axis = 1.2,
        xlab = "Chromosome",
        ylab = expression('F'[ST]), # subscript
        xaxt = "n",
        main = plot.title
    )
    axis(side = 1, at = c(0, ticks[2, ]), labels = FALSE, tick = TRUE)
    axis(side = 1, at = ticks[1, ], labels = c("1H", "2H", "3H", "4H", "5H", "6H", "7H"), tick = FALSE)
    #   Horizontal line for Fst outlier threshold
    abline(h = fst.outlier.threshold, lty = 3, lwd = 1.5, col = "black")
    #   Uncomment if you want legend to be added
    #   This is commented out b/c we are putting 4 plots into one figure
    #   and we only want a single legend. This is here to generate
    #   the legend, so uncomment only when necessary.
    # legend(
    #     "topright",
    #     bty = 'n',
    #     c("Gene hits", "Fst Outlier Threshold"),
    #     lty = c(0, 3), # first slot put nothing, second slot use dotted line
    #     lwd = 1.8,
    #     pch = c(19, NA), # first slot put filled circle, second slot put nothing
    #     col = c("blue", "black"),
    #     cex = 1.4,
    #     pt.cex = 1.6,
    #     pt.lwd = 1.6
    # )
}

highlight.gene.hits <- function(df, plot.title) {
    #   Plot SNPs that hit interesting genes
    par(mar = c(5, 5, 5, 2))
    plot(
        x = df$scaled.BP,
        y = df$FST,
        col = adjustcolor(col = "blue", alpha.f = 0.7),
        xlim = c(0, 4569031868),
        ylim = c(0, 0.8),
        pch = 20,
        cex = 1,
        cex.lab = 1.4,
        cex.main = 1.6,
        cex.axis = 1.2,
        xlab = "Chromosome",
        ylab = expression('F'[ST]), # subscript
        xaxt = "n",
        main = plot.title
    )
}

main <- function() {
    #   Take in user command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    fst.fp <- args[1] # full filepath to Fst file
    genehits.fp <- args[2] # full filepath to file containing gene hits
    outlier.threshold <- args[3] # threshold (i.e. 97.5% would be input as 0.975)
    title <- args[4] # Title for our plot (i.e. "Fst - Elevation below 3000m vs above 3000m")
    out.dir <- args[5] # full filepath to output directory
    out.name <- args[6] # output filename including ".pdf" file extension (i.e. fst_lat.pdf)
    
    #   Read in Fst and gene hits files
    fst.df <- readData(filename = fst.fp)
    genehits.df <- readGeneHits(filename = genehits.fp)
    
    #   Compute cutoff fst value for n% outlier threshold
    #   This is the value that will be used in the plots later on
    #   Example: for the Fst_less3000_vsmore3000_no_NA_physPos.txt file
    #   The 97.5% Fst outlier threshold has a cutoff value of 0.3014466887
    outlier.cutoff <- quantile(x = fst.df$FST, probs = seq(0, 1, as.numeric(outlier.threshold)), na.rm = TRUE)[2]
    
    #   Key for barley pseudomolecular parts positions
    chr1H_part1 <- 312837513
    chr1H_part2 <- 245697919
    chr2H_part1 <- 393532674
    chr2H_part2 <- 374542350
    chr3H_part1 <- 394310633
    chr3H_part2 <- 305400481
    chr4H_part1 <- 355061206
    chr4H_part2 <- 291998952
    chr5H_part1 <- 380865482
    chr5H_part2 <- 289164678
    chr6H_part1 <- 294822070
    chr6H_part2 <- 288558443
    chr7H_part1 <- 325797516
    chr7H_part2 <- 331426484

    #   Barley whole chromosome size
    chr1.w <- chr1H_part1 + chr1H_part2
    chr2.w <- chr2H_part1 + chr2H_part2
    chr3.w <- chr3H_part1 + chr3H_part2
    chr4.w <- chr4H_part1 + chr4H_part2
    chr5.w <- chr5H_part1 + chr5H_part2
    chr6.w <- chr6H_part1 + chr6H_part2
    chr7.w <- chr7H_part1 + chr7H_part2
    
    #   Generate ticks used for plots
    t.chr <- chrTicks(chr1.whole = chr1.w, chr2.whole = chr2.w, chr3.whole = chr3.w, chr4.whole = chr4.w, chr5.whole = chr5.w, chr6.whole = chr6.w, chr7.whole = chr7.w)
    
    #   In Chr_2016 column, replace "chr" with nothing and replace "H" with nothing
    fst.dfr <- replace(df = fst.df)
    
    #   Scale physical positions
    fst.dfrs <- scaledPos(
        df = fst.dfr, c1.w = chr1.w, c2.w = chr2.w,
        c3.w = chr3.w, c4.w = chr4.w, c5.w = chr5.w,
        c6.w = chr6.w, c7.w = chr7.w
    )

    #   Generate Plots
    #   Create subset containing SNPs that hit interesting genes
    genehits.subset <- fst.dfrs[fst.dfrs$SNP %in% genehits.df$SNPs_9k, ]
    
    #   Make plot
    pdf(file = paste0(out.dir, "/", out.name), width = 12, height = 8)
    plot.manhattan(
        df = fst.dfrs,
        fst.outlier.threshold = outlier.cutoff,
        plot.title = title,
        ticks = t.chr
    )
    par(new = TRUE)
    highlight.gene.hits(
        df = genehits.subset,
        plot.title = title
    )
    dev.off()
}

main()
