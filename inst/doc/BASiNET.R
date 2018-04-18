## ------------------------------------------------------------------------
mRNA <- system.file("extdata", "sequences2.fasta", package = "BASiNET")
lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
library(BASiNET)
classification(mRNA,lncRNA, save="example")

## ---- out.width = "400px"------------------------------------------------
knitr::include_graphics("2d.png")

## ------------------------------------------------------------------------
dataPredict <- system.file("extdata", "predict.fasta", package = "BASiNET")
modelPredict <- system.file("extdata", "modelPredict.dat", package = "BASiNET")
library(BASiNET)
classification(predicting=dataPredict,load=modelPredict)

