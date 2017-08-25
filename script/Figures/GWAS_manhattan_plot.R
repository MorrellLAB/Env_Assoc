#!/usr/bin/env Rscript

#   Function to read in Csv file used for plotting
readCsv <- function(filename) {
    data.file <- read.csv(
        file = filename,
        header = TRUE,
        na.strings = "NA"
    )
    return(data.file)
}

#   Plot function
plot.manhattan <- function(x.val, y.val, color.pheno, p.sym, plot.title) {
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
    #   Set path to file
    gwas.filepath <- "/Users/chaochih/Dropbox/Landrace_Environmental_Association/Analyses/GWAS-GAPIT/compiled.5e_4.0.01.v2_physPos.csv"
    gwas.data <- readCsv(filename = gwas.filepath)
    #   In Chr_2016 column, replace "chr" with nothing and replace "H" with nothing
    gwas.data$Chr_2016 <- sub("chr", "", gwas.data$Chr_2016)
    gwas.data$Chr_2016 <- sub("H", "", gwas.data$Chr_2016)
    #   Pull out columns of interest
    gwas.df <- data.frame(
        SNP = gwas.data$SNP,
        PHENO = gwas.data$Phenotype,
        CHR = gwas.data$Chr_2016,
        BP = gwas.data$PhysPos_2016,
        P = gwas.data$P.value
    )

    #   Ignore the following, need to shorten chunks of hardcoded code
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
    chr1.whole <- chr1H_part1 + chr1H_part2
    chr2.whole <- chr2H_part1 + chr2H_part2
    chr3.whole <- chr3H_part1 + chr3H_part2
    chr4.whole <- chr4H_part1 + chr4H_part2
    chr5.whole <- chr5H_part1 + chr5H_part2
    chr6.whole <- chr6H_part1 + chr6H_part2
    chr7.whole <- chr7H_part1 + chr7H_part2
    #   Define ticks
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

    #   Subset chromosomes
    gwas.df.chr1 <- subset(gwas.df, subset = gwas.df$CHR == 1)
    gwas.df.chr2 <- subset(gwas.df, subset = gwas.df$CHR == 2)
    gwas.df.chr3 <- subset(gwas.df, subset = gwas.df$CHR == 3)
    gwas.df.chr4 <- subset(gwas.df, subset = gwas.df$CHR == 4)
    gwas.df.chr5 <- subset(gwas.df, subset = gwas.df$CHR == 5)
    gwas.df.chr6 <- subset(gwas.df, subset = gwas.df$CHR == 6)
    gwas.df.chr7 <- subset(gwas.df, subset = gwas.df$CHR == 7)
    #   Add new column of scaled positions for plotting
    gwas.df.chr1["scaled.BP"] <- gwas.df.chr1$BP
    gwas.df.chr2["scaled.BP"] <- gwas.df.chr2$BP + chr1.whole
    gwas.df.chr3["scaled.BP"] <- gwas.df.chr3$BP + chr1.whole + chr2.whole
    gwas.df.chr4["scaled.BP"] <- gwas.df.chr4$BP + chr1.whole + chr2.whole + chr3.whole
    gwas.df.chr5["scaled.BP"] <- gwas.df.chr5$BP + chr1.whole + chr2.whole + chr3.whole + chr4.whole
    gwas.df.chr6["scaled.BP"] <- gwas.df.chr6$BP + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole
    gwas.df.chr7["scaled.BP"] <- gwas.df.chr7$BP + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole + chr6.whole
    #   Combine data frames
    d.all <- rbind(
        gwas.df.chr1[order(gwas.df.chr1$scaled.BP), ],
        gwas.df.chr2[order(gwas.df.chr2$scaled.BP), ],
        gwas.df.chr3[order(gwas.df.chr3$scaled.BP), ],
        gwas.df.chr4[order(gwas.df.chr4$scaled.BP), ],
        gwas.df.chr5[order(gwas.df.chr5$scaled.BP), ],
        gwas.df.chr6[order(gwas.df.chr6$scaled.BP), ],
        gwas.df.chr7[order(gwas.df.chr7$scaled.BP), ]
    )

    #   Define chr5 inverted region 1
    c5.inv1.start <- 126746171 + chr1.whole + chr2.whole + chr3.whole + chr4.whole
    c5.inv1.end <- 305528375 + chr1.whole + chr2.whole + chr3.whole + chr4.whole
    c5.inv2.start <- 598975647 + chr1.whole + chr2.whole + chr3.whole + chr4.whole
    c5.inv2.end <- 609073831 + chr1.whole + chr2.whole + chr3.whole + chr4.whole
    c2.inv.start <- 267303750 + chr1.whole
    c2.inv.end <- 508786535 + chr1.whole

    #   Generate Plots
    #   Bio1, 6, and 11
    par(mfrow = c(3, 1))
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio1"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio1"]),
        color.pheno = "#005C9E",
        p.sym = 17,
        plot.title = "GWAS BIO1 - Annual Mean Temperature"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio6"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio6"]),
        color.pheno = "#005C9E",
        p.sym = 17,
        plot.title = "GWAS BIO6 - Min Temperature of Coldest Month"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio11"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio11"]),
        color.pheno = "#005C9E",
        p.sym = 17,
        plot.title = "GWAS BIO11 - Mean Temperature of Coldest Quarter"
    )
    
    par(mfrow = c(4,1))
    #   Bio2, 3, 4, and 5
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio2"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio2"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO2 - Mean Diurnal Range (Mean of monthly (max temp - min temp))"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio3"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio3"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO3 - Isothermality (BIO2/BIO7) (* 100)"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio4"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio4"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO4 - Temperature Seasonality (standard deviation *100)"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio5"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio5"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO5 - Max Temperature of Warmest Month"
    )
    
    #   Bio7, 8, 9, and 10
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio7"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio7"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO7 - Temperature Annual Range (BIO5-BIO6)"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio8"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio8"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO8 - Mean Temperature of Wettest Quarter"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio9"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio9"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO9 - Mean Temperature of Driest Quarter"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio10"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio10"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO10 - Mean Temperature of Warmest Quarter"
    )
    
    #   Bio12, 13, 14, and 15
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio12"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio12"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO12 - Annual Precipitation"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio13"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio13"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO13 - Precipitation of Wettest Month"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio14"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio14"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO14 - Precipitation of Driest Month"
    )
    arrows(
        x0 = d.all$scaled.BP[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"],
        y1 = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"]),
        y0 = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"] + (10 * 0.0001636497)),
        length = 0.1, lwd = 2
    )
    text(
        x = d.all$scaled.BP[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"],
        y = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"] + (10 * 0.0001636497) + 0.05),
        labels = "11_20784\n(BIO14, BIO17)"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio15"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio15"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO15 - Precipitation Seasonality (Coefficient of Variation)"
    )
    
    
    #   Bio16, 17, 18, and 19
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio16"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio16"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO16 - Precipitation of Wettest Quarter"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio17"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio17"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO17 - Precipitation of Driest Quarter"
    )
    arrows(
        x0 = d.all$scaled.BP[d.all$SNP == "11_20784" & d.all$PHENO == "bio17"],
        y1 = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio17"]),
        y0 = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio17"] + (10 * 0.0001636497)),
        length = 0.1, lwd = 2
    )
    text(
        x = d.all$scaled.BP[d.all$SNP == "11_20784" & d.all$PHENO == "bio17"],
        y = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio17"] + (10 * 0.0001636497) + 0.05),
        labels = "11_20784\n(BIO14, BIO17)"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio18"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio18"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO18 - Precipitation of Warmest Quarter"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "bio19"],
        y.val = -log10(d.all$P[d.all$PHENO == "bio19"]),
        color.pheno = "#BD1550",
        p.sym = 20,
        plot.title = "GWAS BIO19 - Precipitation of Coldest Quarter"
    )
    
    #   Latitude, longitude, and altitude
    par(mfrow = c(3, 1))
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "Latitude"],
        y.val = -log10(d.all$P[d.all$PHENO == "Latitude"]),
        color.pheno = "orange",
        p.sym = 20,
        plot.title = "GWAS - Latitude"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "Longitude"],
        y.val = -log10(d.all$P[d.all$PHENO == "Longitude"]),
        color.pheno = "orange",
        p.sym = 20,
        plot.title = "GWAS - Longitude"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "altitude"],
        y.val = -log10(d.all$P[d.all$PHENO == "altitude"]),
        color.pheno = "orange",
        p.sym = 20,
        plot.title = "GWAS - Altitude"
    )
    
    #   IC1, 2, and 3
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "IC1"],
        y.val = -log10(d.all$P[d.all$PHENO == "IC1"]),
        color.pheno = "limegreen",
        p.sym = 20,
        plot.title = "GWAS IC1 - Independent components applied to BIO1-19"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "IC2"],
        y.val = -log10(d.all$P[d.all$PHENO == "IC2"]),
        color.pheno = "limegreen",
        p.sym = 20,
        plot.title = "GWAS IC2 - Independent components applied to BIO1-19"
    )
    plot.manhattan(
        x.val = d.all$scaled.BP[d.all$PHENO == "IC3"],
        y.val = -log10(d.all$P[d.all$PHENO == "IC3"]),
        color.pheno = "limegreen",
        p.sym = 20,
        plot.title = "GWAS IC3 - Independent components applied to BIO1-19"
    )
}

main()
