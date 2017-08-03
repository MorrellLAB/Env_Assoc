#!/usr/bin/env Rscript

#   This script installs required packages for a custom script, LDheatmap.R
#   NOTE: use R 3.4.0
#   This script was adapted from ANGSD-wrapper startShiny.R

#   Run this line on the command line before running script: module load R/3.4.0

#   Install required packages for LDheatmap.R
args <- commandArgs(trailingOnly = TRUE)

#   Test if required packages are installed
pkgTest <- function(package) {
    #   Check if package is available
    if(package %in% rownames(installed.packages()) == FALSE) {
        install.packages(package) # If not, install package
    }
}

#   List of required packages
pkgList1 <- c("MASS", "boot", "cluster", "curl", "foreign", "lattice", "Matrix", "mgcv", "nlme", "rpart", "survival")
#   Before install rJava, make sure we are either running R 3.3.2 or R 3.3.3
pkgList2 <- c("LDheatmap", "genetics", "RColorBrewer", "grDevices")

#   Load the packages
batchInstall <- function(pkgList) {
    options(repos = c(CRAN = "http://cran.cnr.berkeley.edu/")) # Set a repo mirror
    for(dep in pkgList) {
        pkgTest(dep) # Test to see if package is installed
    }
    lapply(X = pkgList, FUN = library, character.only = TRUE) # Load the packages that will be used
}

#   Install packages from BioConductor
bioInstall <- function() {
    if("biocLite" %in% rownames(installed.packages()) == FALSE) {
        source("http://bioconductor.org/biocLite.R")
    }
    library(BiocInstaller)
    if("chopsticks" %in% rownames(installed.packages()) == FALSE) {
        biocLite("chopsticks")
    }
    library(chopsticks)
}

#   Driver function
main <- function() {
    LDheatmap <- args[1] # Path to where LDheatmap.R is located
    dir.create(file.path(LDheatmap, ".RLibs", fsep = "/"), showWarnings = FALSE) # Created directory for packages
    .libPaths(paste0(LDheatmap, "/.RLibs")) # Enable this as a place to look for packages
    batchInstall(pkgList1) # Install packages needed for BioConductor
    bioInstall() # Install required packages from BioConductor
    batchInstall(pkgList2) # Install required package, some depend on biocLite and chopsticks
}

main() # Run the program
