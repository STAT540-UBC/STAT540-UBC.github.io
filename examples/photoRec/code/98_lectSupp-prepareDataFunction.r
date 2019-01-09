##########################################################
## helper functions for finding and featuring genes
##########################################################

## handy for isolating and reshaping the data for a small set of
## probesets

## WARNING: assumes prDat and prDes exist in the calling environment!

prepareData <- function(myGenes) {
    miniDat <- t(prDat[myGenes, ])
    miniDat <- data.frame(gExp = as.vector(miniDat),
                          gene = factor(rep(colnames(miniDat), each =
                          nrow(miniDat)), levels = colnames(miniDat)))
    miniDat <- suppressWarnings(data.frame(prDes, miniDat))
    miniDat
}
