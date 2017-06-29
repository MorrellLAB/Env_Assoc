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

#   Set path to file
gwas.filepath <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/GWAS-GAPIT/compiled.5e_4.0.01.v2_physPos.csv"
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
chr1.t <- c(chr1.whole/2, chr1.whole)
chr2.t <- c(chr2.whole/2 + chr1.whole, chr1.whole + chr2.whole)
chr3.t <- c(chr3.whole/2 + chr1.whole + chr2.whole, chr1.whole + chr2.whole + chr3.whole)
chr4.t <- c(chr4.whole/2 + chr1.whole + chr2.whole + chr3.whole, chr4.whole + chr1.whole + chr2.whole + chr3.whole)
chr5.t <- c(chr5.whole/2 + chr1.whole + chr2.whole + chr3.whole + chr4.whole, chr5.whole + chr1.whole + chr2.whole + chr3.whole + chr4.whole)
chr6.t <- c(chr6.whole/2 + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole, chr6.whole + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole)
chr7.t <- c(chr7.whole/2 + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole + chr6.whole, chr7.whole + chr1.whole + chr2.whole + chr3.whole + chr4.whole + chr5.whole + chr6.whole)
ticks <- as.vector(cbind(chr1.t, chr2.t, chr3.t, chr4.t, chr5.t, chr6.t, chr7.t))
#   Define labels
chr1.l <- c(chr1.whole/2, chr1.whole)
chr2.l <- c(chr2.whole/2 + chr1.whole, chr1.whole + chr2.whole)
chr3.l <- c(chr3.whole/2 + chr1.whole + chr2.whole, chr1.whole + chr2.whole + chr3.whole)
chr4.l <- c(chr4.whole/2 + chr1.whole + chr2.whole + chr3.whole, chr4.whole + chr1.whole + chr2.whole + chr3.whole)
chr5.l <- c(chr5.whole/2 + chr1.whole + chr2.whole + chr3.whole + chr4.whole, chr5.whole + chr1.whole + chr2.whole + chr3.whole + chr4.whole)
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
#   Plot function
plot.manhattan <- function(x.val, y.val, color.pheno, p.sym) {
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
        xaxt = "n"
    )
}

plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "Latitude"], y.val = -log10(d.all$P[d.all$PHENO == "Latitude"]), color.pheno = "orange", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "Longitude"], y.val = -log10(d.all$P[d.all$PHENO == "Longitude"]), color.pheno = "orange", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "altitude"], y.val = -log10(d.all$P[d.all$PHENO == "altitude"]), color.pheno = "orange", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "IC1"], y.val = -log10(d.all$P[d.all$PHENO == "IC1"]), color.pheno = "limegreen", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "IC2"], y.val = -log10(d.all$P[d.all$PHENO == "IC2"]), color.pheno = "limegreen", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "IC3"], y.val = -log10(d.all$P[d.all$PHENO == "IC3"]), color.pheno = "limegreen", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio2"], y.val = -log10(d.all$P[d.all$PHENO == "bio2"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio3"], y.val = -log10(d.all$P[d.all$PHENO == "bio3"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio4"], y.val = -log10(d.all$P[d.all$PHENO == "bio4"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio5"], y.val = -log10(d.all$P[d.all$PHENO == "bio5"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio7"], y.val = -log10(d.all$P[d.all$PHENO == "bio7"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio8"], y.val = -log10(d.all$P[d.all$PHENO == "bio8"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio9"], y.val = -log10(d.all$P[d.all$PHENO == "bio9"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio10"], y.val = -log10(d.all$P[d.all$PHENO == "bio10"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio12"], y.val = -log10(d.all$P[d.all$PHENO == "bio12"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio13"], y.val = -log10(d.all$P[d.all$PHENO == "bio13"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio14"], y.val = -log10(d.all$P[d.all$PHENO == "bio14"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio16"], y.val = -log10(d.all$P[d.all$PHENO == "bio16"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio17"], y.val = -log10(d.all$P[d.all$PHENO == "bio17"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio15"], y.val = -log10(d.all$P[d.all$PHENO == "bio15"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio19"], y.val = -log10(d.all$P[d.all$PHENO == "bio19"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio18"], y.val = -log10(d.all$P[d.all$PHENO == "bio18"]), color.pheno = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio1"], y.val = -log10(d.all$P[d.all$PHENO == "bio1"]), color.pheno = "#005C9E", p.sym = 17)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio6"], y.val = -log10(d.all$P[d.all$PHENO == "bio6"]), color.pheno = "#005C9E", p.sym = 17)
par(new = TRUE)
plot.manhattan(x.val = d.all$scaled.BP[d.all$PHENO == "bio11"], y.val = -log10(d.all$P[d.all$PHENO == "bio11"]), color.pheno = "#005C9E", p.sym = 17)
legend("topright", bty = 'n', legend = c("Bio1, 6, 11", "Bio 2-17 except 1, 6, and 11", "Altitude, Latitude, Longitude", "IC1, IC2, IC3"), pch = c(17, 20, 20, 20), col = c("#005C9E", "#BD1550", "orange", "limegreen"), cex = 0.75)
axis(side = 1, at = c(0, ticks[2, ]), labels = FALSE, tick = TRUE)
axis(side = 1, at = ticks[1, ], labels = c("1", "2", "3", "4", "5", "6", "7"), tick = FALSE)
abline(v = c5.inv1.start, col = "gray10", lty = 3)
abline(v = c5.inv1.end, col = "gray10", lty = 3)
abline(v = c5.inv2.start, col = "red", lty = 3)
abline(v = c5.inv2.end, col = "red", lty = 3)
abline(v = c2.inv.start, col = "blue", lty = 3)
abline(v = c2.inv.end, col = "blue", lty = 3)
arrows(x0 = d.all$scaled.BP[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"], 
       y1 = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"]),
       y0 = -log10(d.all$P[d.all$SNP == "11_20784" & d.all$PHENO == "bio14"] + (10 * 0.0001636497)),
       length = 0.1, lwd = 2
)
