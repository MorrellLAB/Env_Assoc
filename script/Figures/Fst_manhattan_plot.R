#!/usr/bin/env Rscript

#   Define function to read in .txt
readData <- function(filename) {
    df <- read.table(
        file = filename,
        header = TRUE,
        sep = "\t",
        na.strings = "NA"
    )
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
    axis(side = 1, at = ticks[1, ], labels = c("1", "2", "3", "4", "5", "6", "7"), tick = FALSE)
    #   Horizontal line for Fst outlier threshold
    abline(h = fst.outlier.threshold, lty = 3, lwd = 1.5, col = "black")
    legend(
        "topright",
        bty = 'n',
        c("Gene hits", "Fst Outlier Threshold"),
        lty = c(0, 3), # first slot put nothing, second slot use dotted line
        pch = c(20, NA), # first slot put filled circle, second slot put nothing
        col = c("blue", "black"),
        cex = 0.75,
        pt.cex = 1.2
    )
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
    #   Set path to files (user provided arguments)
    #   Output file directory
    out.dir <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/Fst/Plots"
    #   Elevation
    e.b3000vsa3000.fp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/Fst/Results/Elevation/Fst_less3000_vsmore3000_no_NA_physPos.txt"
    #   Latitude
    wild.range.bN30vsN30_40.fp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/Fst/Results/Latitude/Fst_moreoreq30to40vsless30_no_NA_physPos.txt"
    wild.range.N30_40vsaN40.fp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/Fst/Results/Latitude/Fst_wild_range30_40_vs_higherLat40_no_NA_physPos.txt"

    #   Read in data
    #   Elevation outliers
    e.df.b3000vsa3000 <- readData(filename = e.b3000vsa3000.fp)
    wild.range.df.bN30vsN30_40 <- readData(filename = wild.range.bN30vsN30_40.fp)
    wild.range.df.N30_40vsaN40 <- readData(filename = wild.range.N30_40vsaN40.fp)

    #   Key for pseudomolecular parts positions
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

    #   Whole chromosome size
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
    e.dfr.b3000vsa3000 <- replace(df = e.df.b3000vsa3000)
    wild.range.dfr.bN30vsN30_40 <- replace(df = wild.range.df.bN30vsN30_40)
    wild.range.dfr.N30_40vsaN40 <- replace(df = wild.range.df.N30_40vsaN40)
    
    #   Scale physical positions
    e.dfrs.b3000vsa3000 <- scaledPos(
        df = e.dfr.b3000vsa3000, c1.w = chr1.w, c2.w = chr2.w,
        c3.w = chr3.w, c4.w = chr4.w, c5.w = chr5.w,
        c6.w = chr6.w, c7.w = chr7.w
    )
    wild.range.dfrs.bN30vsN30_40 <- scaledPos(
        df = wild.range.dfr.bN30vsN30_40, c1.w = chr1.w, c2.w = chr2.w,
        c3.w = chr3.w, c4.w = chr4.w, c5.w = chr5.w,
        c6.w = chr6.w, c7.w = chr7.w
    )
    wild.range.dfrs.N30_40vsaN40 <- scaledPos(
        df = wild.range.dfr.N30_40vsaN40, c1.w = chr1.w, c2.w = chr2.w,
        c3.w = chr3.w, c4.w = chr4.w, c5.w = chr5.w,
        c6.w = chr6.w, c7.w = chr7.w
    )

    #   Generate Plots
    #   Elevation - below 3000m vs above 3000m
    #   Create subset containing SNPs that hit interesting genes
    #   Cold genes
    snp12_30880 <- e.dfrs.b3000vsa3000[e.dfrs.b3000vsa3000$SNP == "12_30880", ]
    snp11_11328 <- e.dfrs.b3000vsa3000[e.dfrs.b3000vsa3000$SNP == "11_11328", ]
    snp12_20187 <- e.dfrs.b3000vsa3000[e.dfrs.b3000vsa3000$SNP == "12_20187", ]
    #   Flowering time genes
    snp12_30867 <- e.dfrs.b3000vsa3000[e.dfrs.b3000vsa3000$SNP == "12_30867", ]
    
    #   Combine genes into single df
    e.gene.hits <- rbind(snp12_30880, snp11_11328, snp12_20187, snp12_30867)
    
    #   Make plot
    pdf(file = paste0(out.dir, "/Fst_elevation_below3000_vs_above3000.pdf"), width = 12, height = 8)
    plot.manhattan(
        df = e.dfrs.b3000vsa3000,
        fst.outlier.threshold = 0.3014466887,
        plot.title = "Fst - Elevation below 3000m vs above 3000m",
        ticks = t.chr
    )
    par(new = TRUE)
    highlight.gene.hits(
        df = e.gene.hits,
        plot.title = "Fst - Elevation below 3000m vs above 3000m"
    )
    dev.off()
    
    #   Latitude
    #   Create subset containing SNPs that hit interesting genes in both wild range data frames
    snp11_20265 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "11_20265", ]
    snpSCRI_RS_142618 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "SCRI_RS_142618", ]
    #   Flowering time genes
    snpBK_12 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "BK_12", ]
    snpBK_16 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "BK_16", ]
    snpBK_15 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "BK_15", ]
    snp12_30871 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "12_30871", ]
    snp12_30872 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "12_30872", ]
    snpBK_13 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "BK_13", ]
    snpBK_14 <- wild.range.dfrs.bN30vsN30_40[wild.range.dfrs.bN30vsN30_40$SNP == "BK_14", ]
    
    #   Combine genes into single data frame
    l.gene.hits <- rbind(snp11_20265, snpSCRI_RS_142618, snpBK_12, snpBK_16, snpBK_15, snp12_30871, snp12_30872, snpBK_13, snpBK_14)
    
    #   Wild range below N30 vs N30-N40
    pdf(file = paste0(out.dir, "/Fst_latitude_belowN30_vs_N30-N40.pdf"), width = 12, height = 8)
    plot.manhattan(
        df = wild.range.dfrs.bN30vsN30_40,
        fst.outlier.threshold = 0.3124423020,
        plot.title = "Fst - Latitude below N30 vs N30-N40",
        ticks = t.chr
    )
    par(new = TRUE)
    highlight.gene.hits(
        df = l.gene.hits,
        plot.title = "Fst - Latitude below N30 vs N30-N40"
    )
    dev.off()
    
    #   Wild range N30-N40 vs above N40
    pdf(file = paste0(out.dir, "/Fst_latitude_N30-N40_vs_aboveN40.pdf"), width = 12, height = 8)
    plot.manhattan(
        df = wild.range.dfrs.N30_40vsaN40,
        fst.outlier.threshold = 0.2601220106,
        plot.title = "Fst - Latitude N30-N40 vs above N40",
        ticks = t.chr
    )
    par(new = TRUE)
    highlight.gene.hits(
        df = l.gene.hits,
        plot.title = "Fst - Latitude N30-N40 vs above N40"
    )
    dev.off()
}

main()
