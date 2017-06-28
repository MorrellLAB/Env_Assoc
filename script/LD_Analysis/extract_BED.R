#!/usr/bin/env Rscript

#   This is a script that extracts intervals from a VCF file
#   Script written by Chaochih Liu
#   June 27, 2017

#   Usage:
#   ./extract_BED.R [file.vcf] [bp] [outputFile]

#   A function to read in VCF file
readVcf <- function(filename) {
    data.file <- read.table(
        file = filename, # passed as an argument
        header = FALSE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
}

writeOutFile <- function(data.file, outName) {
    write.table(
        x = data.file,
        file = outName,
        sep = "\t",
        na = "NA",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
}

#   Run program
main <- function() {
    #   User provided arguments
    args <- commandArgs(trailingOnly = TRUE)
    vcf.filepath <- args[1]
    bp.size <- args[2] # What is the window size upstream/downstream of our SNP?
    outFile <- args[3]
    
    #   Read in file
    vcf.df <- readVcf(filename = vcf.filepath)
    bed.df <- data.frame(
        Chr = vcf.df$V1,
        Start.pos = vcf.df$V2 - (bp.size + 1), # BED file is 0-based
        End.pos = vcf.df$V2 + bp.size
    )
    #   Save file
    writeOutFile(data.file = bed.df, outName = outFile)
}