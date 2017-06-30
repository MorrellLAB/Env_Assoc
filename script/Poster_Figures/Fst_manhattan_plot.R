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

plot.manhattan <- function(x.val, y.val, color.cat, p.sym) {
    p <- plot(
        x = x.val,
        y = y.val,
        col = color.cat,
        xlim = c(0, 4569031868),
        ylim = c(0, 0.8),
        pch = p.sym,
        cex = 1.5,
        xlab = "Chromosome",
        ylab = expression('F'[ST]), # subscript
        xaxt = "n"
    )
}

#   Set path to file
#   Elevation outliers
e.1000vsAll.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_1000vsAll_outliers_physPos.txt"
e.2500vs1000.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_2500vs1000_outliers_physPos.txt"
e.5000vs1000.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_5000_vs1000_outliers_physPos.txt"
e.5000vs2500.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_5000_vs2500_outliers_physPos.txt"
#   Latitude Outliers
l.80samp_10iter.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Latitude/Outliers/Fst_average_80samp_10iter_outliers_physPos.txt"
l.topvsBotLat.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Latitude/Outliers/Fst_top_vs_bottomLat_outWildRange_outliers_physPos.txt"
l.wrvshighLat80.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Latitude/Outliers/Fst_wild_range_vs_higherLat_80samples_outliers_physPos.txt"
l.wrvshighLat.fp <- "/Users/chaochih/Dropbox/Landrace Environmental Assocation/Analyses/Fst/Results/Latitude/Outliers/Fst_wild_range_vs_higherLat_outliers_physPos.txt"

#   Read in data
#   Elevation outliers
e.df.1000vsAll <- readData(filename = e.1000vsAll.fp)
e.df.2500vs1000 <- readData(filename = e.2500vs1000.fp)
e.df.5000vs1000 <- readData(filename = e.5000vs1000.fp)
e.df.5000vs2500 <- readData(filename = e.5000vs2500.fp)
#   Latitdue outliers
l.df.80samp_10iter <- readData(filename = l.80samp_10iter.fp)
l.df.topvsBotLat <- readData(filename = l.topvsBotLat.fp)
l.df.wrvshighLat80 <- readData(filename = l.wrvshighLat80.fp)
l.df.wrvshighLat <- readData(filename = l.wrvshighLat.fp)

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
#   Combine
ticks <- cbind(chr1.t, chr2.t, chr3.t, chr4.t, chr5.t, chr6.t, chr7.t)

#   In Chr_2016 column, replace "chr" with nothing and replace "H" with nothing
e.df.1000vsAll <- replace(df = e.df.1000vsAll)
e.df.2500vs1000 <- replace(df = e.df.2500vs1000)
e.df.5000vs1000 <- replace(df = e.df.5000vs1000)
e.df.5000vs2500 <- replace(df = e.df.5000vs2500)
l.df.80samp_10iter <- replace(df = l.df.80samp_10iter)
l.df.topvsBotLat <- replace(df = l.df.topvsBotLat)
l.df.wrvshighLat80 <- replace(df = l.df.wrvshighLat80)
l.df.wrvshighLat <- replace(df = l.df.wrvshighLat)

#   Scale physical positions
e.dfs.1000vsAll <- scaledPos(df = e.df.1000vsAll, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)
e.dfs.2500vs1000 <- scaledPos(df = e.df.2500vs1000, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)
e.dfs.5000vs1000 <- scaledPos(df = e.df.5000vs1000, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)
e.dfs.5000vs2500 <- scaledPos(df = e.df.5000vs2500, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)
l.dfs.80samp_10iter <- scaledPos(df = l.df.80samp_10iter, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)
l.dfs.topvsBotLat <- scaledPos(df = l.df.topvsBotLat, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)
l.dfs.wrvshighLat80 <- scaledPos(df = l.df.wrvshighLat80, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)
l.dfs.wrvshighLat <- scaledPos(df = l.df.wrvshighLat, c1.w = chr1.whole, c2.w = chr2.whole, c3.w = chr3.whole, c4.w = chr4.whole, c5.w = chr5.whole, c6.w = chr6.whole, c7.w = chr7.whole)

#   Define inverted regions
c2.inv.start <- 267303750 + chr1.whole
c2.inv.end <- 508786535 + chr1.whole
c5.inv1.start <- 126746171 + chr1.whole + chr2.whole + chr3.whole + chr4.whole
c5.inv1.end <- 305528375 + chr1.whole + chr2.whole + chr3.whole + chr4.whole
c5.inv2.start <- 598975647 + chr1.whole + chr2.whole + chr3.whole + chr4.whole
c5.inv2.end <- 609073831 + chr1.whole + chr2.whole + chr3.whole + chr4.whole

#   Generate Plots
plot.manhattan(x.val = l.dfs.80samp_10iter$scaled.BP, y.val = l.dfs.80samp_10iter$FST, color.cat = "#FD93B4", p.sym = 2)
par(new = TRUE)
plot.manhattan(x.val = l.dfs.topvsBotLat$scaled.BP, y.val = l.dfs.topvsBotLat$FST, color.cat = "#B984F3", p.sym = 2)
par(new = TRUE)
plot.manhattan(x.val = l.dfs.wrvshighLat80$scaled.BP, y.val = l.dfs.wrvshighLat80$FST, color.cat = "cornflowerblue", p.sym = 2)
par(new = TRUE)
plot.manhattan(x.val = l.dfs.wrvshighLat$scaled.BP, y.val = l.dfs.wrvshighLat$FST, color.cat = "#F8ED00", p.sym = 2)
par(new = TRUE)
plot.manhattan(x.val = e.dfs.1000vsAll$scaled.BP, y.val = e.dfs.1000vsAll$FST, color.cat = "#BD1550", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = e.dfs.2500vs1000$scaled.BP, y.val = e.dfs.2500vs1000$FST, color.cat = "#005C9E", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = e.dfs.5000vs1000$scaled.BP, y.val = e.dfs.5000vs1000$FST, color.cat = "orange", p.sym = 20)
par(new = TRUE)
plot.manhattan(x.val = e.dfs.5000vs2500$scaled.BP, y.val = e.dfs.5000vs2500$FST, color.cat = "limegreen", p.sym = 20)
legend("topright",
       bty = 'n', 
       legend = c("Elev: 1000 vs All", "Elev: 2500 vs 1000",
                  "Elev: 5000 vs 1000", "Elev: 5000 vs 2500",
                  "Lat: 80 samp 10 iter", "Lat: top vs bottom", 
                  "Lat: wild range vs higher lat 80 samp", 
                  "Lat: wild range vs higher lat"),
       pch = c(2, 2, 2, 2, 20, 20, 20, 20), 
       col = c("#FD93B4", "#B984F3", "cornflowerblue", "#F8ED00",
               "#BD1550", "#005C9E", "orange", "limegreen"), 
       cex = 1,
       ncol = 3)
axis(side = 1, at = c(0, ticks[2, ]), labels = FALSE, tick = TRUE)
axis(side = 1, at = ticks[1, ], labels = c("1", "2", "3", "4", "5", "6", "7"), tick = FALSE)
abline(v = c5.inv1.start, col = "gray10", lty = 3)
abline(v = c5.inv1.end, col = "gray10", lty = 3)
abline(v = c5.inv2.start, col = "red", lty = 3)
abline(v = c5.inv2.end, col = "red", lty = 3)
abline(v = c2.inv.start, col = "blue", lty = 3)
abline(v = c2.inv.end, col = "blue", lty = 3)
#   11_20361 chr 4H AFB2
arrows(x0 = l.dfs.wrvshighLat80$scaled.BP[l.dfs.wrvshighLat80$SNP == "11_20361"], 
       y1 = l.dfs.wrvshighLat80$FST[l.dfs.wrvshighLat80$SNP == "11_20361"],
       y0 = l.dfs.wrvshighLat80$FST[l.dfs.wrvshighLat80$SNP == "11_20361"] + 0.075,
       length = 0.1, lwd = 2
)
#   11_10496 chr 6H AOC
arrows(x0 = e.dfs.1000vsAll$scaled.BP[e.dfs.1000vsAll$SNP == "11_10496"], 
       y1 = e.dfs.1000vsAll$FST[e.dfs.1000vsAll$SNP == "11_10496"],
       y0 = e.dfs.1000vsAll$FST[e.dfs.1000vsAll$SNP == "11_10496"] - 0.075,
       length = 0.1, lwd = 2
)
#   12_30848 chr 5H CBF3
arrows(x0 = l.dfs.topvsBotLat$scaled.BP[l.dfs.topvsBotLat$SNP == "12_30848"], 
       y1 = l.dfs.topvsBotLat$FST[l.dfs.topvsBotLat$SNP == "12_30848"],
       y0 = l.dfs.topvsBotLat$FST[l.dfs.topvsBotLat$SNP == "12_30848"] + 0.075,
       length = 0.1, lwd = 2
)
#   12_20187 chr 1H PEAMT
arrows(x0 = e.dfs.5000vs2500$scaled.BP[e.dfs.5000vs2500$SNP == "12_20187"], 
       y1 = e.dfs.5000vs2500$FST[e.dfs.5000vs2500$SNP == "12_20187"],
       y0 = e.dfs.5000vs2500$FST[e.dfs.5000vs2500$SNP == "12_20187"] - 0.075,
       length = 0.1, lwd = 2
)
#   SCRI_RS_179411 chr 5H SPP
arrows(x0 = l.dfs.80samp_10iter$scaled.BP[l.dfs.80samp_10iter$SNP == "SCRI_RS_179411"], 
       y1 = l.dfs.80samp_10iter$FST[l.dfs.80samp_10iter$SNP == "SCRI_RS_179411"],
       y0 = l.dfs.80samp_10iter$FST[l.dfs.80samp_10iter$SNP == "SCRI_RS_179411"] - 0.075,
       length = 0.1, lwd = 2
)
