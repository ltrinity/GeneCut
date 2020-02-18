#F Cros Analysis Script
#Author: Luke Trinity
#Date: December 19, 2019
#import fcros
library(fcros)
#this method runs fcros with parameters for how many significant genes to output,
#which column to use for test data, and two strings for output filenames
runFCros <- function(numberOfSignificantGenesToOutput, test, outFileName1, outFileName2){
  #identify control
  cont <- c("iPSnorm");
  #identify test from input param
  test <- c(test);
  #default settings
  log2.opt <- 0;
  trim.opt <- 0.25;
  # perform fcros()
  af <- fcros(Hosna.copy, cont, test, log2.opt, trim.opt);
  #output all the genes and fcros values 
  #(top is from initial docs, expanded to output values for all genes in dataset)
  topUpDown <- fcrosTopN(af, numberOfSignificantGenesToOutput);
  #thresholds
  alpha1 <- topUpDown$alpha[1];
  alpha2 <- topUpDown$alpha[2];
  id.down  <- matrix(0, 1);
  fval.down  <- matrix(0, 1);
  foldChange <- matrix(0,1);
  id.up <- matrix(0, 1);
  fval.up  <- matrix(0, 1);
  n <- length(af$FC);
  f.value <- af$f.value;
  idown <- 1;
  iup <- 1;
  #modification here to save fcros value
  for (i in 1:n) {
    if (f.value[i] <= alpha1) { id.down[idown] <- i; 
                                idown <- idown + 1;
    }
    if (f.value[i] >= alpha2) { id.up[iup] <- i;
                                iup <- iup + 1;
    }
  }
  #handle some issues with parsing
  data.down <- gsub(" ","", Hosna.copy[id.down[1:(idown-1)],"symbol" ]);
  fval.down <- f.value[id.down[1:(idown-1)]]; 
  ndown <- nrow(data.down);
  data.up <- gsub(" ","", Hosna.copy[id.up[1:(iup-1)],"symbol" ]);
  fval.up <- f.value[id.up[1:(iup-1)]]; 
  nup <- nrow(data.up)
  #tables for output
  upReg <- data.frame(data.up, fval.up)
  downReg <- data.frame(data.down, fval.down)
  write.table(upReg,file=outFileName1,row.names = FALSE)
  write.table(downReg,file=outFileName2, row.names = FALSE)
}
#may need to set working directory **IMPORTANT**
#setwd('C:/Users/luket/Desktop/GeneCut')
#read csv, equivalent to file -> import dataset
Hosna.copy <- read.csv("input/Hosna-copy.csv")
rownames(Hosna.copy) <- make.names(Hosna.copy[,7],unique=TRUE)
colnames(Hosna.copy) <- c("ID","HUVECnorm","iPSnorm","ECDiff","NnDiff","HFNnorm","symbol","gene","ontologybio","ontologycellular","ontologymicro","alignments")
#run fcros based on specifications
runFCros(48817, "NnDiff",'temp/upregulatedNnMaster.csv','temp/downregulatedNnMaster.csv')
runFCros(48817, "ECDiff",'temp/upregulatedECMaster.csv','temp/downregulatedECMaster.csv')