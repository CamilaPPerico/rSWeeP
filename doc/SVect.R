## --------------------------------------------------------------------------
library(SVect)
baseMatrix <- orthbase(160000,10)

## --------------------------------------------------------------------------
data(datastring)
cat(datastring, file = "exdna.fas", sep = "\n")

## --------------------------------------------------------------------------
return <- SWeeP("exdna.fas",baseMatrix)
distancia <- dist(return, method = "euclidean")
tree <- hclust(distancia, method="ward.D")
plot(tree, hang = -1, cex = 1)

## ----label='Session information', eval=TRUE, echo=FALSE--------------------
sessionInfo()

