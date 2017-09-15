#!/usr/bin/env Rscript

#   Function to read in Csv file used for plotting
#   Might remove in future, this function is only for compiled.5e_4.0.01.v2_physPos.csv file
readCsv <- function(filename) {
    data.file <- read.csv(
        file = filename,
        header = TRUE,
        na.strings = "NA"
    )

    #   In Chr_2016 column, replace "chr" with nothing and replace "H" with nothing
    data.file$Chr_2016 <- sub("chr", "", data.file$Chr_2016)
    data.file$Chr_2016 <- sub("H", "", data.file$Chr_2016)

    #   Pull out columns of interest
    gwas.data <- data.frame(
        SNP = data.file$SNP,
        PHENO = data.file$Phenotype,
        CHR = data.file$Chr_2016,
        BP = data.file$PhysPos_2016,
        P = data.file$P.value
    )
    return(gwas.data)
}

readGAPIT <- function(filename) {
    data.file <- read.csv(
        file = filename,
        header = TRUE,
        na.strings = "NA"
    )

    #   In Chr_2016 column, replace "chr" with nothing and replace "H" with nothing
    data.file$Chr_2016 <- sub("chr", "", data.file$Chr_2016)
    data.file$Chr_2016 <- sub("H", "", data.file$Chr_2016)
    
    #   extract phenotype from filename
    #   NOTE: this only works for the following naming scheme
    #       GAPIT_pheno_GWAS_results_physPos.csv
    #       GAPIT_altitude_GWAS_results_physPos.csv
    pheno.names <- strsplit(basename(path = filename), split = "_")[[1]][2] # after string split, pheno is 2nd element in list
    #   Add extra column of phenotype associated with dataset
    #   We will use apply to read in multiple CSV files into a list of lists
    #   This allows us to track which phenotypes are associated with which lists
    data.file["pheno"] <- pheno.names

    #   Pull out columns of interest
    full.gwas.data <- data.frame(
        SNP = data.file$SNP,
        PHENO = data.file$pheno,
        CHR = data.file$Chr_2016,
        BP = data.file$PhysPos_2016,
        P = data.file$P.value,
        MAF = data.file$maf
    )
    return(full.gwas.data)
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

scalePsuedo <- function(df, chr1.whole, chr2.whole, chr3.whole, chr4.whole, chr5.whole, chr6.whole, chr7.whole) {
    #   Subset chromosomes
    df.chr1 <- subset(df, subset = df$CHR == 1)
    df.chr2 <- subset(df, subset = df$CHR == 2)
    df.chr3 <- subset(df, subset = df$CHR == 3)
    df.chr4 <- subset(df, subset = df$CHR == 4)
    df.chr5 <- subset(df, subset = df$CHR == 5)
    df.chr6 <- subset(df, subset = df$CHR == 6)
    df.chr7 <- subset(df, subset = df$CHR == 7)

    #   Add new column of scaled positions for plotting
    df.chr1["scaled.BP"] <- df.chr1$BP
    df.chr2["scaled.BP"] <- df.chr2$BP + chr1.whole
    df.chr3["scaled.BP"] <- df.chr3$BP + chr1.whole + chr2.whole
    df.chr4["scaled.BP"] <- df.chr4$BP + chr1.whole + chr2.whole + chr3.whole
    df.chr5["scaled.BP"] <- df.chr5$BP + chr1.whole + chr2.whole + chr3.whole + chr4.whole
    df.chr6["scaled.BP"] <- df.chr6$BP + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole
    df.chr7["scaled.BP"] <- df.chr7$BP + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole + chr6.whole

    #   Combine data frames
    d.all.chr <- rbind(
        df.chr1[order(df.chr1$scaled.BP), ],
        df.chr2[order(df.chr2$scaled.BP), ],
        df.chr3[order(df.chr3$scaled.BP), ],
        df.chr4[order(df.chr4$scaled.BP), ],
        df.chr5[order(df.chr5$scaled.BP), ],
        df.chr6[order(df.chr6$scaled.BP), ],
        df.chr7[order(df.chr7$scaled.BP), ]
    )
    return(d.all.chr)
}

#   Plot function
plot.manhattan <- function(x.val, y.val, color.pheno, p.sym, sym.size, plot.title, ticks) {
    p <- plot(
        x = x.val,
        y = y.val,
        col = color.pheno,
        xlim = c(0, 4569031868),
        ylim = c(0, 12),
        pch = p.sym,
        cex = sym.size,
        xlab = "Chromosome",
        ylab = expression(-log[10](italic(p))),
        xaxt = "n",
        main = plot.title
    )
    axis(side = 1, at = c(0, ticks[2, ]), labels = FALSE, tick = TRUE)
    axis(side = 1, at = ticks[1, ], labels = c("1", "2", "3", "4", "5", "6", "7"), tick = FALSE)
    abline(h = -log10(5e-4), lty = 3, lwd = 1.5, col = "black")
    # rect(
    #     xleft = c2.inv.start,
    #     xright = c2.inv.end,
    #     ybottom = 0,
    #     ytop = 12,
    #     col = adjustcolor("gray50", alpha.f = 0.25),
    #     border = NA
    # )
    # rect(
    #     xleft = c5.inv1.start,
    #     xright = c5.inv1.end,
    #     ybottom = 0,
    #     ytop = 12,
    #     col = adjustcolor("gray50", alpha.f = 0.25),
    #     border = NA
    # )
    # rect(
    #     xleft = c5.inv2.start,
    #     xright = c5.inv2.end,
    #     ybottom = 0,
    #     ytop = 12,
    #     col = adjustcolor("gray50", alpha.f = 0.25),
    #     border = NA
    # )
}

main <- function() {
    #   User provided arguments
    #   This file contains only the significant SNPs
    #   Commented out section that only contained the significant SNPs data
    # gwas.filepath <- "/Users/chaochih/Dropbox/Landrace_Environmental_Association/Analyses/GWAS-GAPIT/compiled.5e_4.0.01.v2_physPos.csv"
    #   Directory path to full dataset csv files
    gwas.all.dir <- "/Users/chaochih/Dropbox/Landrace_Environmental_Association/Analyses/GWAS-GAPIT/GWAS_Results_Full_with_PhysPos"

    #   Read in GWAS significant SNPs data
    #   Commented out section that only contained the significant SNPs data
    # gwas.df <- readCsv(filename = gwas.filepath)

    #   Read in list of all GWAS data files
    #   These files contain the full dataset including the significant SNPs for each bioclim variable
    full.d.fp <- list.files(path = gwas.all.dir, pattern = ".csv", full.names = TRUE)
    #   Apply function to read in all files containing GAPIT GWAS data
    full.d <- lapply(
        X = full.d.fp,
        FUN = readGAPIT
    )

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

    #   Define chr2 and chr5 inverted regions
    c2.inv.start <- 267303750 + chr1.w
    c2.inv.end <- 508786535 + chr1.w
    c5.inv1.start <- 126746171 + chr1.w + chr2.w + chr3.w + chr4.w
    c5.inv1.end <- 305528375 + chr1.w + chr2.w + chr3.w + chr4.w
    c5.inv2.start <- 598975647 + chr1.w + chr2.w + chr3.w + chr4.w
    c5.inv2.end <- 609073831 + chr1.w + chr2.w + chr3.w + chr4.w

    #   Generate ticks used for plots
    t.chr <- chrTicks(chr1.whole = chr1.w, chr2.whole = chr2.w, chr3.whole = chr3.w, chr4.whole = chr4.w, chr5.whole = chr5.w, chr6.whole = chr6.w, chr7.whole = chr7.w)

    #   Scale physical positions in data for accurate representation in plots
    #   Commented out section that only contained the significant SNPs data
    # d.ss <- scalePsuedo(
    #     df = gwas.df,
    #     chr1.whole = chr1.w,
    #     chr2.whole = chr2.w,
    #     chr3.whole = chr3.w,
    #     chr4.whole = chr4.w,
    #     chr5.whole = chr5.w,
    #     chr6.whole = chr6.w,
    #     chr7.whole = chr7.w
    # )
    #   Scale physical positions for each list for each phenotype value
    #   Output variable naming scheme: d.all.pheno
    #   i.e. d.all.altitude, d.all.bio1
    for (i in 1:25) {
        assign(
            paste0("d.all.", as.character(unique(full.d[[i]]$PHENO))),
            scalePsuedo(
                df = full.d[[i]],
                chr1.whole = chr1.w,
                chr2.whole = chr2.w,
                chr3.whole = chr3.w,
                chr4.whole = chr4.w,
                chr5.whole = chr5.w,
                chr6.whole = chr6.w,
                chr7.whole = chr7.w
            )
        )
    }

    #   Generate Plots
    #   MAF 0.01 filter was used by Fumi to get list of significant SNPs from GWAS analysis
    #   See: ReadMe.GWAS.landraces in ~/Dropbox/Landrace_Environmental_Association/GWAS-GAPIT
    
    
    ########## Bio1, 6, and 11 in one panel ##########
    par(mfrow = c(3, 1))
    #   All SNPs MAF > 0.01 Bio1
    plot.manhattan(
        x.val = d.all.bio1$scaled.BP[d.all.bio1$MAF > 0.01],
        y.val = -log10(d.all.bio1$P[d.all.bio1$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO1 - Annual Mean Temperature",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio1
    plot.manhattan(
        x.val = d.all.bio1$scaled.BP[d.all.bio1$MAF < 0.01],
        y.val = -log10(d.all.bio1$P[d.all.bio1$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO1 - Annual Mean Temperature",
        ticks = t.chr
    )
    #   Legend only in first plot panel of 3
    #   Commented out because there already exists a screenshot of the legend
    #   The concern was legend would be covering some of the significant SNPs
    #   This makes it easier to control how the final publication figure looks and where to put the legend
    # legend("topright", legend = c("Sig SNPs Threshold", "MAF > 0.01", "MAF < 0.01"), lty = c(3, 0, 0), pch = c(NA, 1, 20), col = c("black", adjustcolor(col = "gray20", alpha.f = 0.5), "red"), bty = "n", y.intersp = 0.25)
    
    #   All SNPs MAF > 0.01 Bio6
    plot.manhattan(
        x.val = d.all.bio6$scaled.BP[d.all.bio6$MAF > 0.01],
        y.val = -log10(d.all.bio6$P[d.all.bio6$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO6 - Min Temperature of Coldest Month",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio6
    plot.manhattan(
        x.val = d.all.bio6$scaled.BP[d.all.bio6$MAF < 0.01],
        y.val = -log10(d.all.bio6$P[d.all.bio6$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO6 - Min Temperature of Coldest Month",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio11
    plot.manhattan(
        x.val = d.all.bio11$scaled.BP[d.all.bio11$MAF > 0.01],
        y.val = -log10(d.all.bio11$P[d.all.bio11$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO11 - Mean Temperature of Coldest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio11
    plot.manhattan(
        x.val = d.all.bio11$scaled.BP[d.all.bio11$MAF < 0.01],
        y.val = -log10(d.all.bio11$P[d.all.bio11$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO11 - Mean Temperature of Coldest Quarter",
        ticks = t.chr
    )
    
    
    ########## Bio2, 3, 4, and 5 in one panel ##########
    par(mfrow = c(4,1))
    #   All SNPs MAF > 0.01 Bio2
    plot.manhattan(
        x.val = d.all.bio2$scaled.BP[d.all.bio2$MAF > 0.01],
        y.val = -log10(d.all.bio2$P[d.all.bio2$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO2 - Mean Diurnal Range (Mean of monthly (max temp - min temp))",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio2
    plot.manhattan(
        x.val = d.all.bio2$scaled.BP[d.all.bio2$MAF < 0.01],
        y.val = -log10(d.all.bio2$P[d.all.bio2$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO2 - Mean Diurnal Range (Mean of monthly (max temp - min temp))",
        ticks = t.chr
    )
    #   Legend only in first plot panel of 4
    #   Commented out because there already exists a screenshot of the legend
    #   The concern was legend would be covering some of the significant SNPs
    #   This makes it easier to control how the final publication figure looks and where to put the legend
    # legend("topright", legend = c("Sig SNPs Threshold", "MAF > 0.01", "MAF < 0.01"), lty = c(3, 0, 0), pch = c(NA, 1, 20), col = c("black", adjustcolor(col = "gray20", alpha.f = 0.5), "red"), bty = "n", y.intersp = 0.10)
    
    #   All SNPs MAF > 0.01 Bio3
    plot.manhattan(
        x.val = d.all.bio3$scaled.BP[d.all.bio3$MAF > 0.01],
        y.val = -log10(d.all.bio3$P[d.all.bio3$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO3 - Isothermality (BIO2/BIO7) (* 100)",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio3
    plot.manhattan(
        x.val = d.all.bio3$scaled.BP[d.all.bio3$MAF < 0.01],
        y.val = -log10(d.all.bio3$P[d.all.bio3$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO3 - Isothermality (BIO2/BIO7) (* 100)",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio4
    plot.manhattan(
        x.val = d.all.bio4$scaled.BP[d.all.bio4$MAF > 0.01],
        y.val = -log10(d.all.bio4$P[d.all.bio4$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO4 - Temperature Seasonality (standard deviation *100)",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio4
    plot.manhattan(
        x.val = d.all.bio4$scaled.BP[d.all.bio4$MAF < 0.01],
        y.val = -log10(d.all.bio4$P[d.all.bio4$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO4 - Temperature Seasonality (standard deviation *100)",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio5
    plot.manhattan(
        x.val = d.all.bio5$scaled.BP[d.all.bio5$MAF > 0.01],
        y.val = -log10(d.all.bio5$P[d.all.bio5$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO5 - Max Temperature of Warmest Month",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio5
    plot.manhattan(
        x.val = d.all.bio5$scaled.BP[d.all.bio5$MAF < 0.01],
        y.val = -log10(d.all.bio5$P[d.all.bio5$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO5 - Max Temperature of Warmest Month",
        ticks = t.chr
    )

    
    ########## Bio7, 8, 9, and 10 in one panel ##########
    par(mfrow = c(4,1))
    #   All SNPs MAF > 0.01 Bio7
    plot.manhattan(
        x.val = d.all.bio7$scaled.BP[d.all.bio7$MAF > 0.01],
        y.val = -log10(d.all.bio7$P[d.all.bio7$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO7 - Temperature Annual Range (BIO5-BIO6)",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio7
    plot.manhattan(
        x.val = d.all.bio7$scaled.BP[d.all.bio7$MAF < 0.01],
        y.val = -log10(d.all.bio7$P[d.all.bio7$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO7 - Temperature Annual Range (BIO5-BIO6)",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio8
    plot.manhattan(
        x.val = d.all.bio8$scaled.BP[d.all.bio8$MAF > 0.01],
        y.val = -log10(d.all.bio8$P[d.all.bio8$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO8 - Mean Temperature of Wettest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio8
    plot.manhattan(
        x.val = d.all.bio8$scaled.BP[d.all.bio8$MAF < 0.01],
        y.val = -log10(d.all.bio8$P[d.all.bio8$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO8 - Mean Temperature of Wettest Quarter",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio9
    plot.manhattan(
        x.val = d.all.bio9$scaled.BP[d.all.bio9$MAF > 0.01],
        y.val = -log10(d.all.bio9$P[d.all.bio9$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO9 - Mean Temperature of Driest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio9
    plot.manhattan(
        x.val = d.all.bio9$scaled.BP[d.all.bio9$MAF < 0.01],
        y.val = -log10(d.all.bio9$P[d.all.bio9$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO9 - Mean Temperature of Driest Quarter",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio10
    plot.manhattan(
        x.val = d.all.bio10$scaled.BP[d.all.bio10$MAF > 0.01],
        y.val = -log10(d.all.bio10$P[d.all.bio10$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO10 - Mean Temperature of Warmest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio10
    plot.manhattan(
        x.val = d.all.bio10$scaled.BP[d.all.bio10$MAF < 0.01],
        y.val = -log10(d.all.bio10$P[d.all.bio10$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO10 - Mean Temperature of Warmest Quarter",
        ticks = t.chr
    )

    
    ########## Bio12, 13, 14, and 15 in one panel ##########
    par(mfrow = c(4, 1))
    #   All SNPs MAF > 0.01 Bio12
    plot.manhattan(
        x.val = d.all.bio12$scaled.BP[d.all.bio12$MAF > 0.01],
        y.val = -log10(d.all.bio12$P[d.all.bio12$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO12 - Annual Precipitation",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio12
    plot.manhattan(
        x.val = d.all.bio12$scaled.BP[d.all.bio12$MAF < 0.01],
        y.val = -log10(d.all.bio12$P[d.all.bio12$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO12 - Annual Precipitation",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio13
    plot.manhattan(
        x.val = d.all.bio13$scaled.BP[d.all.bio13$MAF > 0.01],
        y.val = -log10(d.all.bio13$P[d.all.bio13$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO13 - Precipitation of Wettest Month",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio13
    plot.manhattan(
        x.val = d.all.bio13$scaled.BP[d.all.bio13$MAF < 0.01],
        y.val = -log10(d.all.bio13$P[d.all.bio13$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO13 - Precipitation of Wettest Month",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio14
    plot.manhattan(
        x.val = d.all.bio14$scaled.BP[d.all.bio14$MAF > 0.01],
        y.val = -log10(d.all.bio14$P[d.all.bio14$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO14 - Precipitation of Driest Month",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio14
    plot.manhattan(
        x.val = d.all.bio14$scaled.BP[d.all.bio14$MAF < 0.01],
        y.val = -log10(d.all.bio14$P[d.all.bio14$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO14 - Precipitation of Driest Month",
        ticks = t.chr
    )
    #   Need the following code when including inverted region intervals
    # arrows(
    #     x0 = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"],
    #     y1 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"]),
    #     y0 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"] + (10 * 0.0001636497)),
    #     length = 0.1, lwd = 2
    # )
    # text(
    #     x = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"],
    #     y = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"] + (10 * 0.0001636497) + 0.05),
    #     labels = "11_20784\n(BIO14, BIO17)"
    # )
    
    #   All SNPs MAF > 0.01 Bio15
    plot.manhattan(
        x.val = d.all.bio15$scaled.BP[d.all.bio15$MAF > 0.01],
        y.val = -log10(d.all.bio15$P[d.all.bio15$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO15 - Precipitation Seasonality (Coefficient of Variation)",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio15
    plot.manhattan(
        x.val = d.all.bio15$scaled.BP[d.all.bio15$MAF < 0.01],
        y.val = -log10(d.all.bio15$P[d.all.bio15$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO15 - Precipitation Seasonality (Coefficient of Variation)",
        ticks = t.chr
    )


    ########## Bio16, 17, 18, and 19 in one panel ##########
    par(mfrow = c(4, 1))
    #   All SNPs MAF > 0.01 Bio16
    plot.manhattan(
        x.val = d.all.bio16$scaled.BP[d.all.bio16$MAF > 0.01],
        y.val = -log10(d.all.bio16$P[d.all.bio16$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO16 - Precipitation of Wettest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio16
    plot.manhattan(
        x.val = d.all.bio16$scaled.BP[d.all.bio16$MAF < 0.01],
        y.val = -log10(d.all.bio16$P[d.all.bio16$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO16 - Precipitation of Wettest Quarter",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio17
    plot.manhattan(
        x.val = d.all.bio17$scaled.BP[d.all.bio17$MAF > 0.01],
        y.val = -log10(d.all.bio17$P[d.all.bio17$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO17 - Precipitation of Driest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio17
    plot.manhattan(
        x.val = d.all.bio17$scaled.BP[d.all.bio17$MAF < 0.01],
        y.val = -log10(d.all.bio17$P[d.all.bio17$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO17 - Precipitation of Driest Quarter",
        ticks = t.chr
    )
    #   Need the following code when including inverted region intervals
    # arrows(
    #     x0 = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"],
    #     y1 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"]),
    #     y0 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"] + (10 * 0.0001636497)),
    #     length = 0.1, lwd = 2
    # )
    #     x = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"],
    #     y = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"] + (10 * 0.0001636497) + 0.05),
    #     labels = "11_20784\n(BIO14, BIO17)"
    # )
    
    #   All SNPs MAF > 0.01 Bio18
    plot.manhattan(
        x.val = d.all.bio18$scaled.BP[d.all.bio18$MAF > 0.01],
        y.val = -log10(d.all.bio18$P[d.all.bio18$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO18 - Precipitation of Warmest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio18
    plot.manhattan(
        x.val = d.all.bio18$scaled.BP[d.all.bio18$MAF < 0.01],
        y.val = -log10(d.all.bio18$P[d.all.bio18$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO18 - Precipitation of Warmest Quarter",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio19
    plot.manhattan(
        x.val = d.all.bio19$scaled.BP[d.all.bio19$MAF > 0.01],
        y.val = -log10(d.all.bio19$P[d.all.bio19$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS BIO19 - Precipitation of Coldest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio19
    plot.manhattan(
        x.val = d.all.bio19$scaled.BP[d.all.bio19$MAF < 0.01],
        y.val = -log10(d.all.bio19$P[d.all.bio19$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS BIO19 - Precipitation of Coldest Quarter",
        ticks = t.chr
    )

    
    ########## Latitude, longitude, and altitude in one panel ##########
    par(mfrow = c(3, 1))
    #   All SNPs MAF > 0.01 Latitude
    plot.manhattan(
        x.val = d.all.Latitude$scaled.BP[d.all.Latitude$MAF > 0.01],
        y.val = -log10(d.all.Latitude$P[d.all.Latitude$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS - Latitude",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Latitude
    plot.manhattan(
        x.val = d.all.Latitude$scaled.BP[d.all.Latitude$MAF < 0.01],
        y.val = -log10(d.all.Latitude$P[d.all.Latitude$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS - Latitude",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Longitude
    plot.manhattan(
        x.val = d.all.Longitude$scaled.BP[d.all.Longitude$MAF > 0.01],
        y.val = -log10(d.all.Longitude$P[d.all.Longitude$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS - Longitude",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Longitude
    plot.manhattan(
        x.val = d.all.Longitude$scaled.BP[d.all.Longitude$MAF < 0.01],
        y.val = -log10(d.all.Longitude$P[d.all.Longitude$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS - Longitude",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Altitude
    plot.manhattan(
        x.val = d.all.altitude$scaled.BP[d.all.altitude$MAF > 0.01],
        y.val = -log10(d.all.altitude$P[d.all.altitude$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS - Altitude",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Altitude
    plot.manhattan(
        x.val = d.all.altitude$scaled.BP[d.all.altitude$MAF < 0.01],
        y.val = -log10(d.all.altitude$P[d.all.altitude$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS - Altitude",
        ticks = t.chr
    )

    
    ########## IC1, 2, and 3 in one panel ##########
    par(mfrow = c(3, 1))
    #   All SNPs MAF > 0.01 IC1
    plot.manhattan(
        x.val = d.all.IC1$scaled.BP[d.all.IC1$MAF > 0.01],
        y.val = -log10(d.all.IC1$P[d.all.IC1$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS IC1 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 IC1
    plot.manhattan(
        x.val = d.all.IC1$scaled.BP[d.all.IC1$MAF < 0.01],
        y.val = -log10(d.all.IC1$P[d.all.IC1$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS IC1 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 IC2
    plot.manhattan(
        x.val = d.all.IC2$scaled.BP[d.all.IC2$MAF > 0.01],
        y.val = -log10(d.all.IC2$P[d.all.IC2$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS IC2 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 IC2
    plot.manhattan(
        x.val = d.all.IC2$scaled.BP[d.all.IC2$MAF < 0.01],
        y.val = -log10(d.all.IC2$P[d.all.IC2$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS IC2 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 IC3
    plot.manhattan(
        x.val = d.all.IC3$scaled.BP[d.all.IC3$MAF > 0.01],
        y.val = -log10(d.all.IC3$P[d.all.IC3$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 1,
        sym.size = 0.75,
        plot.title = "GWAS IC3 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 IC3
    plot.manhattan(
        x.val = d.all.IC3$scaled.BP[d.all.IC3$MAF < 0.01],
        y.val = -log10(d.all.IC3$P[d.all.IC3$MAF < 0.01]),
        color.pheno = "red",
        p.sym = 20,
        sym.size = 1,
        plot.title = "GWAS IC3 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
}

main()
