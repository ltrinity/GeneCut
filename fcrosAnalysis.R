#F Cros Analysis Script
#Author: Luke Trinity
#Date: December 19, 2019
library(fcros)
#this method runs fcros with parameters for how many significant genes to output,
#which column to use for test data, and two strings for outpur filenames
runFCros <- function(numberOfSignificantGenesToOutput, test, outFileName1, outFileName2){
  cont <- c("iPSnorm");
  #modularity in test column
  test <- c(test);
  log2.opt <- 0;
  trim.opt <- 0.25;
  # perform fcros()
  af <- fcros(Hosna.copy, cont, test, log2.opt, trim.opt);
  af$idnames
  #output all the genes and fcros values based on input number
  topUpDown <- fcrosTopN(af, numberOfSignificantGenesToOutput);
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
  #modification here to save f cros value
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
Hosna.copy <- read.table("input/Hosna-copy.csv", header=TRUE,fill = TRUE, sep=",")
#load in the data and update the column names
rownames(Hosna.copy) <- make.names(Hosna.copy[,7],unique=TRUE)
colnames(Hosna.copy) <- c("ID","HUVECnorm","iPSnorm","ECDiff","NnDiff","HFNnorm","symbol","gene","ontologybio","ontologycellular","ontologymicro","alignments")
#run fcros based on specifications
runFCros(49400, "NnDiff",'upregulatedNnMaster.csv','downregulatedNnMaster.csv')
runFCros(49400, "ECDiff",'upregulatedECMaster.csv','downregulatedECMaster.csv')