#!/usr/bin/env Rscript

library(data.table)

#   Function to read in MAF file with two columns: SNP and Maf
readMaf <- function(filename) {
    data.file <- fread(
        input = filename,
        sep = "\t",
        header = TRUE
    )
    return(data.file)
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

fixMaf <- function(full.df, new_maf.df) {
    merged.df <- merge(
        x = full.df,
        y = new_maf.df,
        by.x = "SNP",
        by.y = "maf_snp",
        all.x = TRUE, # Rows that do not have a match will remain in dataframe
        all.y = TRUE # Rows that do not have a match with x will not be added to x
    )
    return(merged.df)
}

scalePsuedo <- function(df, maf.data, chr1.whole, chr2.whole, chr3.whole, chr4.whole, chr5.whole, chr6.whole, chr7.whole) {
    #   Add column of correct MAF values before scaling
    #   Calls on function defined above
    df.merged <- fixMaf(full.df = df, new_maf.df = maf.data)
    
    #   Subset chromosomes
    df.chr1 <- subset(df.merged, subset = df.merged$CHR == 1)
    df.chr2 <- subset(df.merged, subset = df.merged$CHR == 2)
    df.chr3 <- subset(df.merged, subset = df.merged$CHR == 3)
    df.chr4 <- subset(df.merged, subset = df.merged$CHR == 4)
    df.chr5 <- subset(df.merged, subset = df.merged$CHR == 5)
    df.chr6 <- subset(df.merged, subset = df.merged$CHR == 6)
    df.chr7 <- subset(df.merged, subset = df.merged$CHR == 7)

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
plot.manhattan <- function(df, maf.threshold, plot.title, ticks) {
    #   Extract all rows where Chromosome is an odd number
    df.odd <- df[!as.numeric(as.character(df$CHR)) %% 2 == 0, ]
    all.x.gth.odd <- df.odd$scaled.BP[df.odd$correct_maf >= maf.threshold]
    all.y.gth.odd <- -log10(df.odd$P[df.odd$correct_maf >= maf.threshold])
    df.even <- df[as.numeric(as.character(df$CHR)) %% 2 == 0, ]
    all.x.gth.even <- df.even$scaled.BP[df.even$correct_maf >= maf.threshold]
    all.y.gth.even <- -log10(df.even$P[df.even$correct_maf >= maf.threshold])
    #   Subset all chromosomes less than maf.threshold
    all.x.lth <- df$scaled.BP[df$correct_maf < maf.threshold]
    all.y.lth <- -log10(df$P[df$correct_maf < maf.threshold])
    
    #   Plot all SNPs greater than threshold
    p <- plot(
        x = all.x.gth.even,
        y = all.y.gth.even,
        col = adjustcolor(col = "gray70", alpha.f = 0.9),
        xlim = c(0, 4569031868),
        ylim = c(0, 12),
        pch = 1,
        cex = 0.7,
        cex.lab = 1.4,
        cex.main = 1.6,
        cex.axis = 1.2,
        xlab = "Chromosome",
        ylab = expression(-log[10](italic(p))),
        xaxt = "n",
        main = plot.title,
        las = 1
    )
    par(new = TRUE)
    #   Make all SNPs greater than threshold that are odd chromosomes black
    plot(
        x = all.x.gth.odd,
        y = all.y.gth.odd,
        col = adjustcolor(col = "gray10", alpha.f = 0.6),
        xlim = c(0, 4569031868),
        ylim = c(0, 12),
        pch = 1,
        cex = 0.7,
        cex.lab = 1.4,
        cex.main = 1.6,
        cex.axis = 1.2,
        xlab = "Chromosome",
        ylab = expression(-log[10](italic(p))),
        xaxt = "n",
        main = plot.title,
        las = 1
    )
    par(new = TRUE)
    #   Plot all SNPs less than threshold
    plot(
        x = all.x.lth,
        y = all.y.lth,
        col = adjustcolor(col = "dodgerblue", alpha.f = 0.9),
        xlim = c(0, 4569031868),
        ylim = c(0, 12),
        pch = 20,
        cex = 1,
        cex.lab = 1.4,
        cex.main = 1.6,
        cex.axis = 1.2,
        xlab = "Chromosome",
        ylab = expression(-log[10](italic(p))),
        xaxt = "n",
        main = plot.title,
        las = 1
    )
    axis(side = 1, at = c(0, ticks[2, ]), labels = FALSE, tick = TRUE)
    axis(side = 1, at = ticks[1, ], labels = c("1", "2", "3", "4", "5", "6", "7"), tick = FALSE, cex.axis = 1.2)
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
    ##################################################
    ##########   User provided arguments    ##########
    ##################################################
    
    #   Full filepath to updated (correct) MAF values for each SNP
    #   Updated Jan 2018
    maf.fp <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Contributors/Fumi/Fumi_Env_GWAS/GWAS.landraces/trans_myGD.v2.sorted.maf"
    #   Directory path to full dataset csv files
    gwas.all.dir <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Analyses/GWAS-GAPIT/GWAS_Results_Full_with_PhysPos"
    out.dir <- "/Users/chaochih/Dropbox/Projects/Landrace_Environmental_Association/Manuscript/EnvAssoc_Manuscript_2017-09-28/Figures_Tables/candiate_figs/GWAS All SNPs"
    
    #   Read in list of all GWAS data files
    #   These files contain the full dataset including the significant SNPs for each bioclim variable
    full.d.fp <- list.files(path = gwas.all.dir, pattern = ".csv", full.names = TRUE)
    #   Apply function to read in all files containing GAPIT GWAS data
    full.d <- lapply(
        X = full.d.fp,
        FUN = readGAPIT
    )
    
    #   Read in updated maf values
    maf.df <- readMaf(filename = maf.fp)
    #   Rename columns
    colnames(maf.df) <- c("maf_snp", "correct_maf")

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

    #   Scale physical positions for each list for each phenotype value
    #   Output variable naming scheme: d.all.pheno
    #   i.e. d.all.altitude, d.all.bio1
    for (i in 1:25) {
        
        assign(
            paste0("d.all.", as.character(unique(full.d[[i]]$PHENO))),
            scalePsuedo(
                df = full.d[[i]],
                maf.data = maf.df,
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

    ##########################################
    ##########   Generate Plots     ##########
    ##########################################
    
    #   MAF 0.01 filter was used by Fumi to get list of significant SNPs from GWAS analysis
    #   See: ReadMe.GWAS.landraces in ~/Dropbox/Landrace_Environmental_Association/GWAS-GAPIT
    
    ########## Bio1, 6, and 11 in one panel - Cold Tolerance ##########
    pdf(file = paste0(out.dir, "/GWAS_bio1_bio6_bio11_Manhattan_Plots.pdf"), width = 12, height = 14)
        par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Bio1
        plot.manhattan(
        df = d.all.bio1,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO1 - Annual Mean Temperature",
        ticks = t.chr
        )
        #   Legend only in first plot panel of 3
        #   Commented out because there already exists a screenshot of the legend
        #   The concern was legend would be covering some of the significant SNPs
        #   This makes it easier to control how the final publication figure looks and where to put the legend
        #legend("topright", legend = c("Sig SNPs Threshold", "MAF < 0.01"), lty = c(3, 0), pch = c(NA, 20), col = c("black", adjustcolor(col = "red", alpha.f = 0.9)), bty = "n", y.intersp = 0.5)
        
        #   Bio6
        plot.manhattan(
            df = d.all.bio6,
            maf.threshold = 0.01,
            plot.title = "GWAS BIO6 - Min Temperature of Coldest Month",
            ticks = t.chr
        )
        
        #   Bio11
        plot.manhattan(
            df = d.all.bio11,
            maf.threshold = 0.01,
            plot.title = "GWAS BIO11 - Mean Temperature of Coldest Quarter",
            ticks = t.chr
        )
    dev.off()
    
    
    ########## Bio9, 14, and 17 in one panel - Drought Tolerance ##########
    pdf(file = paste0(out.dir, "/GWAS_bio9_bio14_bio17_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Bio9
    plot.manhattan(
        df = d.all.bio9,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO9 - Mean Temperature of Driest Quarter",
        ticks = t.chr
    )
    
    #   Bio14
    plot.manhattan(
        df = d.all.bio14,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO14 - Precipitation of Driest Month",
        ticks = t.chr
    )
    #   Need the following code when adding SNP names
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
    
    #   Bio17
    plot.manhattan(
        df = d.all.bio17,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO17 - Precipitation of Driest Quarter",
        ticks = t.chr
    )
    #   Need the following code when adding SNP names
    # arrows(
    #     x0 = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"],
    #     y1 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"]),
    #     y0 = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"] + (10 * 0.0001636497)),
    #     length = 0.1, lwd = 2
    # )
    # text(
    #     x = d.ss$scaled.BP[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"],
    #     y = -log10(d.ss$P[d.ss$SNP == "11_20784" & d.ss$PHENO == "bio17"] + (10 * 0.0001636497) + 0.05),
    #     labels = "11_20784\n(BIO14, BIO17)"
    # )
    dev.off()
    
    ########## Bio2, 3, and 4 in one panel ##########
    pdf(file = paste0(out.dir, "/GWAS_bio2_bio3_bio4_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Bio2
    plot.manhattan(
        df = d.all.bio2,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO2 - Mean Diurnal Range (Mean of monthly (max temp - min temp))",
        ticks = t.chr
    )
    
    #   Bio3
    plot.manhattan(
        df = d.all.bio3,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO3 - Isothermality (BIO2/BIO7) (* 100)",
        ticks = t.chr
    )
    
    #   Bio4
    plot.manhattan(
        df = d.all.bio4,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO4 - Temperature Seasonality (standard deviation *100)",
        ticks = t.chr
    )
    dev.off()
    
    ########## Bio5, 7, and 8 in one panel ##########
    pdf(file = paste0(out.dir, "/GWAS_bio5_bio7_bio8_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Bio5
    plot.manhattan(
        df = d.all.bio5,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO5 - Max Temperature of Warmest Month",
        ticks = t.chr
    )
    
    #   Bio7
    plot.manhattan(
        df = d.all.bio7,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO7 - Temperature Annual Range (BIO5-BIO6)",
        ticks = t.chr
    )
    
    #   Bio8
    plot.manhattan(
        df = d.all.bio8,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO8 - Mean Temperature of Wettest Quarter",
        ticks = t.chr
    )
    dev.off()

    
    ########## Bio 10, 12, and 13 in one panel ##########
    pdf(file = paste0(out.dir, "/GWAS_bio10_bio12_bio13_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Bio10
    plot.manhattan(
        df = d.all.bio10,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO10 - Mean Temperature of Warmest Quarter",
        ticks = t.chr
    )
    
    #   Bio12
    plot.manhattan(
        df = d.all.bio12,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO12 - Annual Precipitation",
        ticks = t.chr
    )
    
    #   Bio13
    plot.manhattan(
        df = d.all.bio13,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO13 - Precipitation of Wettest Month",
        ticks = t.chr
    )
    dev.off()

    
    ########## Bio15, 16, and 18 in one panel ##########
    pdf(file = paste0(out.dir, "/GWAS_bio15_bio16_bio18_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Bio15
    plot.manhattan(
        df = d.all.bio15,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO15 - Precipitation Seasonality (Coefficient of Variation)",
        ticks = t.chr
    )
    
    #   Bio16
    plot.manhattan(
        df = d.all.bio16,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO16 - Precipitation of Wettest Quarter",
        ticks = t.chr
    )
    
    #   Bio18
    plot.manhattan(
        df = d.all.bio18,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO18 - Precipitation of Warmest Quarter",
        ticks = t.chr
    )
    dev.off()


    ########## Bio19 in one panel ##########
    pdf(file = paste0(out.dir, "/GWAS_bio19_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Bio19
    plot.manhattan(
        df = d.all.bio19,
        maf.threshold = 0.01,
        plot.title = "GWAS BIO19 - Precipitation of Coldest Quarter",
        ticks = t.chr
    )
    dev.off()

    
    ########## Latitude, longitude, and altitude in one panel ##########
    pdf(file = paste0(out.dir, "/GWAS_Latitude_Longitude_Altitude_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   Latitude
    plot.manhattan(
        df = d.all.Latitude,
        maf.threshold = 0.01,
        plot.title = "GWAS - Latitude",
        ticks = t.chr
    )
    
    #   Longitude
    plot.manhattan(
        df = d.all.Longitude,
        maf.threshold = 0.01,
        plot.title = "GWAS - Longitude",
        ticks = t.chr
    )
    
    #   Altitude
    plot.manhattan(
        df = d.all.altitude,
        maf.threshold = 0.01,
        plot.title = "GWAS - Altitude",
        ticks = t.chr
    )
    dev.off()

    
    ########## IC1, 2, and 3 in one panel ##########
    pdf(file = paste0(out.dir, "/GWAS_IC1_IC2_IC3_Manhattan_Plots.pdf"), width = 12, height = 14)
    par(mfrow = c(3, 1), mar = c(5, 5, 5, 2))
    #   IC1
    plot.manhattan(
        df = d.all.IC1,
        maf.threshold = 0.01,
        plot.title = "GWAS IC1 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    
    #   IC2
    plot.manhattan(
        d.all.IC2,
        maf.threshold = 0.01,
        plot.title = "GWAS IC2 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    
    #   IC3
    plot.manhattan(
        df = d.all.IC3,
        maf.threshold = 0.01,
        plot.title = "GWAS IC3 - Independent components applied to BIO1-19",
        ticks = t.chr
    )
    dev.off()
}

main()
