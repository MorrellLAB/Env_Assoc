#!/usr/bin/env Rscript

#   This is a script that used for plot the geographic distribution of samples located >3000m location
#   Written by Li Lei
#   Sep8, 2017

#   Usage: ./allelic_3Dgeo_plot.R /Users/lilei/Downloads/northeast_TIF/SRTM_NE_250m.tif /Users/lilei/Downloads/lat_log_ele_myGD.v2.txt 12_11529 Alt 12_11529_elevationm_maxFs-Fst=0.6493 /Users/lilei/Downloads 12_11529_elevation_maxFst

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

###subset the files:

subsetSamplePlot <- function(fulldata,mergedRaster,title,outDir,outName){
  #subset the target SNP dataset
  targetset <- data.frame(
    Latitude = fulldata$Latitude,
    Longitude = fulldata$Longitude,
    Altitude = fulldata$altitude
  )
  #subset the sample
  above3k <- fulldata[(fulldata$altitude >3000),]
  above3kTar <- data.frame(
    Latitude = above3k$Latitude,
    Longitude = above3k$Longitude,
    Altitude = above3k$altitude
  )
  #plot
  pdf(file = paste0(outDir, "/", outName, ".pdf"),width=16.00, height=8.67)
  plot(mergedRaster,xlim=range(-10,145), ylim=range(0,55),ylab="Latitude",xlab="Longitude") #plot the map with elevation
  points(targetset$Longitude,targetset$Latitude,col=adjustcolor("blue",alpha.f = 0.9),pch=1,cex=1.2)
  points(above3kTar$Longitude,above3kTar$Latitude,col=adjustcolor("deeppink1",alpha.f = 0.6),pch=19,cex=1.2)
  abline(h=30,lty=2,col="grey60",lwd=1.0)
  abline(h=40,lty=2,col="grey60",lwd=1.0)
  #abline(h=30,lty=2,col="grey30",lwd=1.5)
  #abline(h=40,lty=2,col="grey30",lwd=1.5)
  #abline(h=49,lty=2,col="grey30",lwd=1.5)
  title(main=title)
  dev.off()
}


#   Driver function
main <- function(){
  rasterNE <- args[1] # NE_tif file
  #rasterSE <- args[2] # SE_tif file
  input <- args[2] # fulldataset including GPS and genotype data
  title <- args[3] #title of the plot
  outDir <- args[4] #output directory
  outName <- args[5] #output name
  mr <- merge_raster(NE_tif=rasterNE)
  fullset <- readFile(input)
  plot_allelic <- subsetSamplePlot(fulldata=fullset,mergedRaster=mr,title=title,outDir,outName)
}

main() # Run the program






