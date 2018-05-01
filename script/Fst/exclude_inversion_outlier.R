#!/usr/bin/env Rscript

#   This is a script that used for excluding the SNPs in the putative inversion region based on LD definitation;
############################################
#      Inversion regions                  #
# chr2H: 267,303,750 - 508,786,535        #
# chr5H: invt1: 126,746,171 - 305,528,375 #
# chr5H: invt2: 598,975,647-609,073,831   #
###########################################
#   Written by Li Lei
#   05-01-2018

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

#file <- "/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/AdegenetIndividual/FstLoci_out_Longitude_physPos.txt"
#data <- readFile(file)
#head(data)
#inv <- subsetFile(data)
#head(inv)
#outlier <- getoutliers(inv)
#head(outlier)
#quantile(inv$Fst,probs = seq(0, 1, by= 0.005))

# subset files:
subsetFile <- function(delfile) {
subdel <- delfile[(delfile$Chr_2016!="chr2H" & delfile$Chr_2016!="chr5H"),] #get the data without chr2H and 5H
#get the SNP outside the inversions:
#chr2H
chr2H <- delfile[(delfile$Chr_2016 =="chr2H"),]
subdel1 <- chr2H[(chr2H$PhysPos_2016 < 267303750 | chr2H$PhysPos_2016 > 508786535),]
#chr5H
chr5H <- delfile[(delfile$Chr_2016 =="chr5H"),]
subdel2 <- chr5H[(chr5H$PhysPos_2016 <126746171 |  (chr5H$PhysPos_2016 >305528375 & chr5H$PhysPos_2016 <598975647) |  chr5H$PhysPos_2016 >609073831),]
#combined all the files:
noninv <- rbind(subdel,subdel1,subdel2)
#tosse the lines with NA:
noninv <- na.omit(noninv)
return(noninv)
}

#get outliers
getoutliers <- function(noninv) {
###get outliers with 97.5% threshold:
threshold <- quantile(noninv$Fst, probs = seq(0, 1, by= 0.005))[196] #take 97.5% as threhold
outliers <- noninv[(noninv$Fst >threshold),]
return(outliers)
}
#write file:
writeOutFile <- function(data.file, outName) {
  write.table(
    x = data.file,
    file = outName,
    sep = "\t",
    na = "NA",
    row.names = F,
    col.names = T,
    quote = FALSE
  )
}

#   Driver function
main <- function(){
  input <- args[1] # input file
  exinvfile <- args[2] # excluded inversion file
  out <- args[3] #outliers
  fullset <- readFile(input)
  noninv <- subsetFile(delfile=fullset)
  outfile <- getoutliers(noninv=noninv)
  writeOutFile(data.file=noninv, outName=exinvfile)
  writeOutFile(data.file=outfile, outName=out)
}

main() # Run the program
