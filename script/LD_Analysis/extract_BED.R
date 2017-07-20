#!/usr/bin/env Rscript

#   This is a script that extracts intervals from a VCF file and outputs a BED file
#   Script written by Chaochih Liu
#   June 27, 2017

#   Usage:
#   ./extract_BED.R [file.vcf] [bp] [outputFile]

#   User provided arguments:
#       1) [file.vcf] is the VCF file we want to extract intervals from
#       2) [bp] is the window size upstream/downstream of our SNP
#           i.e. if we want a window size of 100Kb around our SNP,
#           [bp] input would be 50000 to get 50Kb upstream and 50Kb downstream
#           of our SNP.
#       3) [outputFile] is the full filepath, including filename, to our output file

#   A function to read in VCF file
readVcf <- function(filename) {
    data.file <- read.table(
        file = filename, # passed as an argument
        header = FALSE, # First line is a header
        fill = TRUE, # Fill empty fields with NAs
        na.strings = "NA"
    )
}

#   Function to save data to output file
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
    vcf.filepath <- args[1] # 1) VCF file we want to extract interval from
    bp.size <- args[2] # 2) What is the window size upstream/downstream of our SNP?
    outFile <- args[3] # 3) Full filepath to our output file (include filename)

    #   Read in file
    vcf.df <- readVcf(filename = vcf.filepath)

    #   Take VCF file and create interval that is n bp around SNP position
    #   Expected output is BED file
    bed.df <- data.frame(
        Chr = vcf.df$V1,
        Start.pos = vcf.df$V2 - (bp.size + 1), # BED file is 0-based
        End.pos = vcf.df$V2 + bp.size
    )

    #   Save file
    writeOutFile(data.file = bed.df, outName = outFile)
}
