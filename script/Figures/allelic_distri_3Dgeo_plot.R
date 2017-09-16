#!/usr/bin/env Rscript

#   This is a script that used for plot the geographic distribution of the alleles with highest or lowest Fst
#   Written by Li Lei
#   Sep8, 2017

#   Usage: ./allelic_3Dgeo_plot.R /Users/lilei/Downloads/northeast_TIF/SRTM_NE_250m.tif /Users/lilei/Downloads/southeast_TIF/SRTM_SE_250m.tif /Users/lilei/Downloads/combined_gps_Land_6152_SNPs_AB.txt SCRI_RS_149971 SCRI_RS_149971_top_bottom_maxFs-Fst=0.740634006t /Users/lilei/Downloads SCRI_RS_149971_top_bottom_maxFst

library(sp)
library(raster)
library(rgdal)

#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)

#   Function to raster the data
#raster_file <- function(filename) {
#  imported_raster <- raster(filename)
#}

# Function to read and merge the raster
# the NE_tif and SE_tif are eleveation data download from the http://gisweb.ciat.cgiar.org/TRMM/SRTM_Resampled_250m/

#merge_raster <- function (NE_tif, SE_tif){
#  NE_raster <- raster(NE_tif)
#  SE_raster <- raster(SE_tif)
#  rm <- merge(NE_raster,SE_raster)
#  return(rm)
#}

#Raster the tif file:
merge_raster <- function (NE_tif){
  NE_raster <- raster(NE_tif)
  #SE_raster <- raster(SE_tif)
  #rm <- merge(NE_raster,SE_raster)
  return(NE_raster)
}

###read files:
#   Function to read in data
readFile <- function(filename) {
  data.file <- read.delim(
    file = filename, # passed as an argument
    header = TRUE, # First line is a header
    sep="\t", #seperate with tab
    fill = TRUE, # Fill empty fields with NAs
    na.strings = "NA"
  )
  return(data.file)
}

###subset the files, here use Fumi's file include 748 samples rather than Ana's file

subsetSamplePlot <- function(targetSNP,fulldata,mergedRaster,title,outDir,outName){
  #subset the target SNP dataset
  targetset <- data.frame(
    Accession = fulldata$Taxa,
    Latitude = fulldata$Latitude,
    Longitude = fulldata$Longitude,
    Altitude = fulldata$altitude,
    #grepl return a bunch of false and one "true" matched to the pattern, then use fulldata to only extract the column with true value
    targetSNP = fulldata[grepl(x = colnames(fulldata), pattern = targetSNP)]
    )
  #split into different allelc types: 0=AA, 2=BB
    target_AA <- targetset[(targetset[5]=="0"),]
    target_BB <- targetset[(targetset[5]=="2"),]
  #plot
    pdf(file = paste0(outDir, "/", outName, ".pdf"))
    plot(mergedRaster) #plot the map with elevation
    points(target_AA$Longitude,target_AA$Latitude,col=adjustcolor("blue",alpha.f = 0.9),pch=1,cex=1.2)
    points(target_BB$Longitude,target_BB$Latitude,col=adjustcolor("deeppink1",alpha.f = 0.5),pch=20,cex=1.5)
    abline(h=30,lty=2,col="grey30",lwd=1.5)
    abline(h=46,lty=2,col="grey30",lwd=1.5)
    abline(h=49,lty=2,col="grey30",lwd=1.5)
    title(main=title)
    dev.off()
  }


###debug
#createdf <- function(dataset, SNPid)
#  targetset <- data.frame(
#    Accession = fullset$Accession_ID,
#    Country = fullset$Country,
#    Latitude = fullset$Latitude,
#    Longitude = fullset$Longitude,
#    targetSNP = fullset[grepl(x = colnames(fullset), pattern = SNPid)]
#)


#   Driver function
main <- function(){
  rasterNE <- args[1] # NE_tif file
  #rasterSE <- args[2] # SE_tif file
  input <- args[2] # fulldataset including GPS and genotype data
  targetSNP <- args[3] #target SNPs
  title <- args[4] #title of the plot
  outDir <- args[5] #output directory
  outName <- args[6] #output name
  mr <- merge_raster(NE_tif=rasterNE, SE_tif=rasterSE)
  fullset <- readFile(input)
  plot_allelic <- subsetSamplePlot(targetSNP= targetSNP,fulldata=fullset,mergedRaster=mr,title=title,outDir,outName)
}

main() # Run the program







