#!/usr/bin/env Rscript

#   This script is a small function that transposes dataframes
#   and outputs a new transposed file
#   Script written by Chaochih Liu
#   June 13, 2016

#   To run the script: ./transpose_data.R [data.txt] [out directory]


#   Read in data and create data frame from it
#   Transpose data
#   Read in file without row.names = 1 to prevent manually shifting first row of cells over one cell
readFile <- function(filename) {
    transposedData <- t(read.delim(file = filename))
    return(transposedData)
}

#   Write file to outfile
writeOutFile <- function(transposedData, filename, outDirectory) {
    inputName <- unlist(strsplit(x = filename, split = ".txt"))
    outputName <- paste(inputName, "transposed.txt", sep = "_")
    write.table(x = transposedData,
                file = outputName,
                quote = FALSE,
                sep = "\t",
                eol = "\n",
                col.names = FALSE)
}

#   Run the script
createOutFile <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    originalData <- args[1]
    outDir <- args[2]
    formatData <- readFile(originalData)
    output <- writeOutFile(
        transposedData = formatData,
        filename = originalData,
        outDirectory = outDir
    )
}

createOutFile() # Run program
