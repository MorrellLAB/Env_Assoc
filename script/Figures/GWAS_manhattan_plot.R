#!/usr/bin/env Rscript

#   Function to read in Csv file used for plotting
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
plot.manhattan <- function(x.val, y.val, color.pheno, p.sym, plot.title, ticks) {
    p <- plot(
        x = x.val,
        y = y.val,
        col = color.pheno,
        xlim = c(0, 4569031868),
        ylim = c(0, 12),
        pch = p.sym,
        cex = 1.5,
        xlab = "Chromosome",
        ylab = expression(-log[10](italic(p))),
        xaxt = "n",
        main = plot.title
    )
    axis(side = 1, at = c(0, ticks[2, ]), labels = FALSE, tick = TRUE)
    axis(side = 1, at = ticks[1, ], labels = c("1", "2", "3", "4", "5", "6", "7"), tick = FALSE)
    abline(h = -log10(5e-4), lty = 3, lwd = 1.5, col = "black")
    rect(
        xleft = c2.inv.start,
        xright = c2.inv.end,
        ybottom = 0,
        ytop = 12,
        col = adjustcolor("gray50", alpha.f = 0.25),
        border = NA
    )
    rect(
        xleft = c5.inv1.start,
        xright = c5.inv1.end,
        ybottom = 0,
        ytop = 12,
        col = adjustcolor("gray50", alpha.f = 0.25),
        border = NA
    )
    rect(
        xleft = c5.inv2.start,
        xright = c5.inv2.end,
        ybottom = 0,
        ytop = 12,
        col = adjustcolor("gray50", alpha.f = 0.25),
        border = NA
    )
}

main <- function() {
    #   User provided arguments
    #   This file contains only the significant SNPs
    gwas.filepath <- "/Users/chaochih/Dropbox/Landrace_Environmental_Association/Analyses/GWAS-GAPIT/compiled.5e_4.0.01.v2_physPos.csv"
    #   Directory path to full dataset csv files
    gwas.all.dir <- "/Users/chaochih/Dropbox/Landrace_Environmental_Association/Analyses/GWAS-GAPIT/GWAS_Results_Full_with_PhysPos"

    #   Read in GWAS significant SNPs data
    gwas.df <- readCsv(filename = gwas.filepath)

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
    d.ss <- scalePsuedo(
        df = gwas.df,
        chr1.whole = chr1.w,
        chr2.whole = chr2.w,
        chr3.whole = chr3.w,
        chr4.whole = chr4.w,
        chr5.whole = chr5.w,
        chr6.whole = chr6.w,
        chr7.whole = chr7.w
    )
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
    
    
    #   Bio1, 6, and 11 in one panel
    par(mfrow = c(3, 1))
    #   All SNPs MAF > 0.01 Bio1
    plot.manhattan(
        x.val = d.all.bio1$scaled.BP[d.all.bio1$MAF > 0.01],
        y.val = -log10(d.all.bio1$P[d.all.bio1$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 20,
        plot.title = "GWAS BIO1 - Annual Mean Temperature",
        ticks = t.chr
    )
    par(new = TRUE)
    #   Significant SNPs Bio1
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio1"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio1"]),
        color.pheno = "#005C9E",
        p.sym = 17,
        plot.title = "GWAS BIO1 - Annual Mean Temperature",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio1
    plot.manhattan(
        x.val = d.all.bio1$scaled.BP[d.all.bio1$MAF < 0.01],
        y.val = -log10(d.all.bio1$P[d.all.bio1$MAF < 0.01]),
        color.pheno = "gold1",
        p.sym = 4,
        plot.title = "GWAS BIO1 - Annual Mean Temperature",
        ticks = t.chr
    )
    legend("topright", legend = c("Sig SNPs", "MAF > 0.01", "MAF < 0.01"), pch = c(17, 20, 4), col = c("#005c9E", adjustcolor(col = "gray20", alpha.f = 0.5), "gold1"), bty = "n", y.intersp = 0.25)
    
    #   All SNPs MAF > 0.01 Bio6
    plot.manhattan(
        x.val = d.all.bio6$scaled.BP[d.all.bio6$MAF > 0.01],
        y.val = -log10(d.all.bio6$P[d.all.bio6$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 20,
        plot.title = "GWAS BIO6 - Min Temperature of Coldest Month",
        ticks = t.chr
    )
    par(new = TRUE)
    #   Significant SNPs Bio6
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio6"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio6"]),
        color.pheno = "#005C9E",
        p.sym = 17,
        plot.title = "GWAS BIO6 - Min Temperature of Coldest Month",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio6
    plot.manhattan(
        x.val = d.all.bio6$scaled.BP[d.all.bio6$MAF < 0.01],
        y.val = -log10(d.all.bio6$P[d.all.bio6$MAF < 0.01]),
        color.pheno = "gold1",
        p.sym = 4,
        plot.title = "GWAS BIO6 - Min Temperature of Coldest Month",
        ticks = t.chr
    )
    
    #   All SNPs MAF > 0.01 Bio11
    plot.manhattan(
        x.val = d.all.bio11$scaled.BP[d.all.bio11$MAF > 0.01],
        y.val = -log10(d.all.bio11$P[d.all.bio11$MAF > 0.01]),
        color.pheno = adjustcolor(col = "gray20", alpha.f = 0.5),
        p.sym = 20,
        plot.title = "GWAS BIO11 - Mean Temperature of Coldest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   Significant SNPs Bio11
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio11"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio11"]),
        color.pheno = "#005C9E",
        p.sym = 17,
        plot.title = "GWAS BIO11 - Mean Temperature of Coldest Quarter",
        ticks = t.chr
    )
    par(new = TRUE)
    #   All SNPs MAF < 0.01 Bio11
    plot.manhattan(
        x.val = d.all.bio11$scaled.BP[d.all.bio11$MAF < 0.01],
        y.val = -log10(d.all.bio11$P[d.all.bio11$MAF < 0.01]),
        color.pheno = "gold1",
        p.sym = 4,
        plot.title = "GWAS BIO11 - Mean Temperature of Coldest Quarter",
        ticks = t.chr
    )
    
    ############################## Start Here next day #################################
    
    
    par(mfrow = c(4,1))
    #   Bio2, 3, 4, and 5
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio2"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio2"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO2 - Mean Diurnal Range (Mean of monthly (max temp - min temp))"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio3"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio3"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO3 - Isothermality (BIO2/BIO7) (* 100)"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio4"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio4"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO4 - Temperature Seasonality (standard deviation *100)"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio5"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio5"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO5 - Max Temperature of Warmest Month"
    )

    #   Bio7, 8, 9, and 10
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio7"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio7"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO7 - Temperature Annual Range (BIO5-BIO6)"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio8"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio8"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO8 - Mean Temperature of Wettest Quarter"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio9"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio9"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO9 - Mean Temperature of Driest Quarter"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio10"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio10"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO10 - Mean Temperature of Warmest Quarter"
    )

    #   Bio12, 13, 14, and 15
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio12"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio12"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO12 - Annual Precipitation"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio13"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio13"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO13 - Precipitation of Wettest Month"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio14"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio14"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO14 - Precipitation of Driest Month"
    )
    arrows(
        x0 = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"],
        y1 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"]),
        y0 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"] + (10 * 0.0001636497)),
        length = 0.1, lwd = 2
    )
    text(
        x = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"],
        y = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio14"] + (10 * 0.0001636497) + 0.05),
        labels = "11_20784\n(BIO14, BIO17)"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio15"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio15"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO15 - Precipitation Seasonality (Coefficient of Variation)"
    )


    #   Bio16, 17, 18, and 19
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio16"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio16"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO16 - Precipitation of Wettest Quarter"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio17"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio17"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO17 - Precipitation of Driest Quarter"
    )
    arrows(
        x0 = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"],
        y1 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"]),
        y0 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"] + (10 * 0.0001636497)),
        length = 0.1, lwd = 2
    )
    text(
        x = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"],
        y = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"] + (10 * 0.0001636497) + 0.05),
        labels = "11_20784\n(BIO14, BIO17)"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio18"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio18"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO18 - Precipitation of Warmest Quarter"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "bio19"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "bio19"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO19 - Precipitation of Coldest Quarter"
    )

    #   Latitude, longitude, and altitude
    par(mfrow = c(3, 1))
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "Latitude"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "Latitude"]),
        color.pheno = "orange",
        p.sym = 20,
        plot.title = "GWAS - Latitude"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "Longitude"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "Longitude"]),
        color.pheno = "orange",
        p.sym = 20,
        plot.title = "GWAS - Longitude"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "altitude"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "altitude"]),
        color.pheno = "orange",
        p.sym = 20,
        plot.title = "GWAS - Altitude"
    )

    #   IC1, 2, and 3
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "IC1"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "IC1"]),
        color.pheno = "limegreen",
        p.sym = 20,
        plot.title = "GWAS IC1 - Independent components applied to BIO1-19"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "IC2"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "IC2"]),
        color.pheno = "limegreen",
        p.sym = 20,
        plot.title = "GWAS IC2 - Independent components applied to BIO1-19"
    )
    plot.manhattan(
        x.val = d.ss$scaled.BP[d.ss$PHENO == "IC3"],
        y.val = -log10(d.ss$P[d.ss$PHENO == "IC3"]),
        color.pheno = "limegreen",
        p.sym = 20,
        plot.title = "GWAS IC3 - Independent components applied to BIO1-19"
    )
}

main()
