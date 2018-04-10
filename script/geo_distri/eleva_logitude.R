#!/usr/bin/env Rscript

#   This is a script that used for plot the geographic distribution of the alleles with highest or lowest Fst according to Fumi's suggestion (Latitude versus the elevation)
#   Written by Li Lei
#   April,2018

#   Usage: ./allelic_3Dgeo_plot.R /Users/lilei/Downloads/northeast_TIF/SRTM_NE_250m.tif /Users/lilei/Downloads/southeast_TIF/SRTM_SE_250m.tif /Users/lilei/Downloads/combined_gps_Land_6152_SNPs_AB.txt SCRI_RS_149971 SCRI_RS_149971_top_bottom_maxFs-Fst=0.740634006t /Users/lilei/Downloads SCRI_RS_149971_top_bottom_maxFst

#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)


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

#filenames <- "/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Geo_distri/data/lat_log_ele_myGD.v2.txt"

#mydata <- readFile(filenames)
#head(mydata)

###subset the files:

subsetSamplePlot <- function(targetSNP,fulldata,alleleStatus,title,outDir,outName){
  #subset the target SNP dataset
  targetset <- data.frame(
    Accession = fulldata$Taxa,
    Latitude = fulldata$Latitude,
    Longitude = fulldata$Longitude,
    Altitude = fulldata$altitude,
    #grepl return a bunch of false and one "true" matched to the pattern, then use fulldata to only extract the column with true value
    targetSNP = fulldata[grepl(x = colnames(fulldata), pattern = targetSNP)]
  )
  #split into different allelc types:
  if(alleleStatus == "Ref") {
    target_AA <- targetset[(targetset[5]=="0"),]
    target_BB <- targetset[(targetset[5]=="2"),]
  }
  else{
    target_AA <- targetset[(targetset[5]=="2"),]
    target_BB <- targetset[(targetset[5]=="0"),]
  }
  #plot
  pdf(file = paste0(outDir, "/", outName, ".pdf"),width=16.00, height=8.67)
  plot(target_AA$Longitude,target_AA$Altitude,col=adjustcolor("blue",alpha.f = 0.9),pch=1,cex=1.2,ylab="Elevation",xlab="Longitude",xlim=c(0,140),ylim=c(0,6000)) #plot the map with elevation
  points(target_BB$Longitude,target_BB$Altitude,col=adjustcolor("deeppink1",alpha.f = 0.5),pch=20,cex=1.5)
  #abline(h=49,lty=2,col="grey30",lwd=1.5)
  title(main=title)
  dev.off()
}

#plot_allelic <- subsetSamplePlot(targetSNP= "X12_11529",fulldata=mydata,alleleStatus="Alt",title="12_11529",outDir="/Users/lilei/Downloads/",outName="X12_11529")

#mydata[grepl(x = colnames(mydata), pattern = "12_11529")]

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
  alleleStatus <- args[1] # Ref or Alt
  #rasterSE <- args[2] # SE_tif file
  input <- args[2] # fulldataset including GPS and genotype data
  targetSNP <- args[3] #target SNPs
  title <- args[4] #title of the plot
  outDir <- args[5] #output directory
  outName <- args[6] #output name
  fullset <- readFile(input)
  plot_allelic <- subsetSamplePlot(targetSNP= targetSNP,fulldata=fullset,alleleStatus=alleleStatus,title=title,outDir,outName)
}

main() # Run the program