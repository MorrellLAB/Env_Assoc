#!/usr/bin/env Rscript

#   This is a script that generates LD heatmap
#   The concept of this script was drawn from Peter Morrell's script:
    #   https://github.com/pmorrell/Utilities/blob/master/LDheatmap/CAP_heatmap_7H.r
#   Written by Chaochih Liu
#   September 26, 2016

#   Required arguments:
    #   1) genoData.txt: genotype data frame that has SNPs sorted by SNP names
    #   2) physPos.txt: file that includes two columns named - "Query_SNP" and "PhysPos"
    #   3) heatmap plot name: this will be the title of the heatmap figure
    #   4) outPrefix: outfile prefix (do not include file extension)
    #   5) outDir: full path to where output files go
    #   6) includeSnpName: onlty takes in two args - include/exclude

#   To run: ./LDheatmap.R <genoData.txt> <physPos.txt> <Plot Name> <Out File Prefix> <out directory>

library(LDheatmap)
library(genetics)
library(RColorBrewer)
library(grDevices)
require(chopsticks) # used for LDheatmap function

#   Read in genotype data
#   Input data frame should be sorted by SNP names
readGenoFile <- function(filename) {
    #   we want the SNPs as columns & sample names (i.e. WBDC) as rows
    #   the file we read in has sample names (i.e. WBDC) as columns & SNPs as rows
    genoData <- t(read.delim(file = filename,
                             row.names = 1, # include row names
                             header = TRUE, # include column names
                             na.strings = "NN") # 'NN' is missing data, convert to NA
    )
    return(genoData)
}

#   Read in file containing SNP names and physical position info
#   Input file should already be sorted by Query_SNP
readPhysPos <- function(filename) {
    tmp.data.file <- read.delim(
        file = filename, # passed as an argument
        header = TRUE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
    data.file <- data.frame(Query_SNP = tmp.data.file$Query_SNP,
                            PhysPos = tmp.data.file$PhysPos)
    return(data.file)
}

#   Convert genotype format to LDheatmap function compatible format
#   Example: from AA to A/A
makeGeno <- function(genoData) {
    cat("Running makeGenotypes...", sep = "\n")
    compatible.df <- makeGenotypes(data = genoData, sep = "")
    cat("Done running makeGenotypes...", sep = "\n")
    return(compatible.df)
}

#   Remove columns that contain only NA's
findAllNA <- function(df.column) {
    #   This data frame only has 315 columns
    #   Test if NA's are present in columns
    search.all.na <- sum(is.na(df.column))
    finding <- search.all.na == 315
    return(finding)
}

#   Fix compatible genotypes file
#   Columns with too many NA's tend to not convert successfully
#   Remove these columns and keep a log of SNPs removed
#   Function to find incompatible columns
#   Paul Hoffman assisted in writing code for this function
findIncompatible <- function(df.column) {
    #   Remove NA's temporarily for TRUE/FALSE tests
    no.NA.df <- df.column[!is.na(df.column)]
    #   Use grepl to return logical vector
    searches <- grepl(pattern="[ACTG]/[ACTG]", x = no.NA.df) # Get a vector of TRUE/FALSE for being compatible
    search.summary <- sum(!searches) # Invert all TRUES and FALSES to better find failures using sum()
    #   This works because all trues inverted become zero
    #   So, a single fail means that we get a sum of greater than zero
    #   If we didn't invert, a single fail would be n-1 rather than n
    #   And I don't want to test for n-1, testing for 0 is easier
    search.bool <- as.logical(search.summary) # Convert that number
    return(search.bool)
}

#   Generate heatmap for r^2 values
hm.r2 <- function(genoData, PhysPos, plotName, snpName, outName, directory) {
    #   Generate heatmap color palette
    heatmapColors <- brewer.pal(n = 9, name = "YlOrRd")
    #   Output file naming
    outputName <- paste("HM", "r2", sep = "-")
    svg(filename = paste0(directory, "/", outName, "-", outputName, ".svg"))
    #   Run LDheatmap
    heatmap.r2 <- LDheatmap(gdat = genoData,
                            genetic.distances = PhysPos$PhysPos,
                            distances = "physical",
                            LDmeasure = "r",
                            add.map = FALSE,
                            add.key = TRUE,
                            SNP.name = snpName,
                            flip = FALSE,
                            color = heatmapColors,
                            title = plotName)
    dev.off()
    return(heatmap.r2)
}

#   Generate heatmap for D' values
hm.Dprime <- function(genoData, PhysPos, plotName, snpName, outName, directory) {
    #   Generate heatmap color palette
    heatmapColors <- brewer.pal(n = 9, name = "YlOrRd")
    #   Output file naming
    outputName <- paste("HM", "Dprime", sep = "-")
    svg(filename = paste0(directory, "/", outName, "-", outputName, ".svg"))
    #   Generate D' plot
    heatmap.D <- LDheatmap(gdat = genoData,
                           genetic.distances = PhysPos$PhysPos,
                           distances = "physical",
                           LDmeasure = "D'",
                           add.map = FALSE,
                           add.key = TRUE,
                           SNP.name = snpName,
                           flip = FALSE,
                           color = heatmapColors,
                           title = plotName)
    dev.off()
    return(heatmap.D)
}

#   Save info to spreadsheet
outCsv <- function(df, rowNames, outName) {
    write.csv(x = df,
              file = outName,
              row.names = rowNames)
}

#   Write data to output file
outCompatibleFile <- function(df, outName, outDirectory) {
    #   Name output .svg file
    name <- paste0(outDirectory, "/", outName, ".txt")
    write.table(x = df,
                file = name,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = TRUE, # To prevent top row from shifting due to empty first cell
                row.names = TRUE) # Want to keep SNP names
}

#   Save heatmap to .svg file
outFile <- function(outName, directory, heatmap, LDcalc) {
    #   Output file naming
    outputName <- paste("HM", LDcalc, sep = "_")
    #   Name output .svg file
    name <- paste0(directory, "/", outName, "_", outputName, ".txt")
    heatmap.df <- as.data.frame(x = heatmap$LDmatrix, row.names = NULL)
    write.table(x = heatmap.df,
                file = name,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = NA, # To prevent top row from shifting due to empty first cell
                row.names = TRUE) # Want to keep SNP names
}

#   Driver function
main <- function() {
    #   Take command line arguments
    #   Stores arguments into a vector
    args <- commandArgs(trailingOnly = TRUE)
    #   User provided arguments
    genoData <- args[1] # 1) genotype data frame that has SNPs sorted
    physPosData <- args[2] # 2) file that includes two columns - "Query_SNP" and "PhysPos"
    plotName <- args[3] # 3) heatmap plot name
    outPrefix <- args [4] # 4) outfile prefix (do not include file extension)
    outDir <- args[5] # 5) where do our output files go?
    includeSnpName <- args[6] # 6) Do we want to include SNP names in our plot? (include/exclude)

    #   Read in genotype data and SNP_BAC.txt file
    genoFile <- readGenoFile(filename = genoData)
    physPos <- readPhysPos(filename = physPosData)
    X.names <- physPos # We don't need to fix names in this dataset, just use physPos

    #   Create output and plot names
    outname <- paste0(outDir, "/", outPrefix, "_SNP_info")

    #   Convert data to compatible format
    geno.converted <- makeGeno(genoData = genoFile)
    #   Remove columns with only NA values
    cat("Removing empty columns filled with NA's...", sep = "\n")
    results.na <- apply(X = geno.converted, # genotype data we are fixing
                        MARGIN = 2,
                        FUN = findAllNA)
    #   Keep a record of empty columns
    no.data.cols <- names(x = which(x = results.na))
    cat("Saving samples with empty columns to spreadsheet.", sep = "\n")
    tryCatch({
        #   Empty columns will go in 3rd sheet only if they exist
        outCsv(df = no.data.cols,
                rowNames = FALSE,
                outName = paste0(outname, "-empty_cols.csv"))
    }, error = function(e) {
        cat("No empty columns found.", sep = "\n")
        }
    )
    #   Data frame excluding emtpy columns
    results.na.removed <- geno.converted[, !results.na]

    #   Remove incompatible genotype columns
    cat("Removing incompatible columns...", sep = "\n")
    results <- apply(X = results.na.removed, # genotype data we are fixing
                     MARGIN = 2, # Applied over columns
                     FUN = findIncompatible) # function to find incompatible columns
    #   Keep a record of column names that failed
    failed.samples <- names(x = which(x = results))
    cat("Saving failed samples to spreadsheet.", sep = "\n")
    tryCatch({
        outCsv(df = failed.samples,
                rowNames = FALSE,
                outName = paste0(outname, "-failed_snps.csv"))
    }, error = function(e) {
        cat("No failed samples.", sep = "\n")
        }
    )

    #   Remove SNPs with emtpy columns from PhysPos file
    X.names.noEmpty <- X.names[!(X.names$Query_SNP %in% no.data.cols),]
    #   Remove failed SNPs
    X.names.filtered <- X.names.noEmpty[!(X.names.noEmpty$Query_SNP %in% failed.samples),]

    #   Keep compatible data
    pass.samples <- results.na.removed[, !results] # Pull out all rows but specific columns
    cat("Saving compatible data used for heatmap to file.", sep = "\n")
    outCompatibleFile(
        df = pass.samples,
        outName = paste0(outPrefix, "_compatibleSnps"),
        outDirectory = outDir
    )

    #   Do we include SNP names or exclude SNP names?
    if(includeSnpName == "exclude") {
        cat("Exclude SNP names from plots.", sep = "\n")
        snpname <- FALSE
    } else {
        cat("Include SNP names in plots.", sep = "\n")
        snpname <- X.names.filtered$Query_SNP
    }

    #   Message for starting r^2 heatmap
    cat("LDheatmap - r2 starting...", sep = "\n")
    #   r^2 heatmap
    plot.r2 <- hm.r2(genoData = pass.samples,
                     PhysPos = X.names.filtered,
                     plotName = plotName,
                     snpName = snpname,
                     outName = outPrefix,
                     directory = outDir)
    cat("LDheatmap - r2 done.", sep = "\n")

    #   D' heatmap
    cat("LDheatmap - D' starting...", sep = "\n")
    plot.D <- hm.Dprime(genoData = pass.samples,
                        PhysPos = X.names.filtered,
                        plotName = plotName,
                        snpName = snpname,
                        outName = outPrefix,
                        directory = outDir)
    cat("LDheatmap - D' done", sep = "\n")

    cat("Saving files to out directory...", sep = "\n")
    #   Save r2 plot to .svg file
    outFile(outName = outPrefix,
            directory = outDir,
            heatmap = plot.r2,
            LDcalc = "r2")
    #   Save D' plot to .svg file
    outFile(outName = outPrefix,
            directory = outDir,
            heatmap = plot.D,
            LDcalc = "D_prime")
    cat("Done.", sep = "\n")
}

#   Run the program
main()
