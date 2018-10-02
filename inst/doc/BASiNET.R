## ------------------------------------------------------------------------
mRNA <- system.file("extdata", "sequences2.fasta", package = "BASiNET")
lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
library(BASiNET)
classification(mRNA,lncRNA, save="example")

## ---- out.width = "400px"------------------------------------------------
knitr::include_graphics("2d.png")

## ------------------------------------------------------------------------
mRNApredict <- system.file("extdata", "sequences2-predict.fasta", package = "BASiNET")
lncRNApredict <- system.file("extdata", "sequences-predict.fasta", package = "BASiNET")
modelPredict <- system.file("extdata", "modelPredict.dat", package = "BASiNET")
library(BASiNET)
classification(mRNApredict,lncRNApredict,load=modelPredict)

