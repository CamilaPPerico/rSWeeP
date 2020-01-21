fas2mat <- function(allSeq){
  allSeq <- as.list(strsplit(as.character(allSeq), ""))
  allSeq<-vapply(allSeq,aa2mat,FUN.VALUE =double(160000))
  t(allSeq)
}
aa2mat <- function(xseq) {
  ## defines the slide window size
  nucol <- 5
  ## creates a matrix that which row advances on number until the the sequence end
  matrixSlide <- dna2list(xseq, nucol)
  ## creates a matrix with the index position of each couple caught by sliding window
  xy <- cbind(aa2num2(matrixSlide[, 1:2]), aa2num2(matrixSlide[, 4:nucol]))
  ## verify and removes any null position on index matrix
  pinot <- (1 - matrix(xy <= 0, ncol = 2))
  inot <- which(matrix(as.numeric(!(pinot[, 1] * pinot[, 2])), ncol = 1) != 0)
  xy[inot, ] <- NA
  ## plot an occurrence matrix of 400x400 using the xy variable as index positions and removes invalid
  ## values
  inds <- c((xy[, 2] - 1) * 400 + (xy[, 1]))
  inds[inds > (160000)] <- NA
  M <- matrix(0, ncol = 400, nrow = 400)
  M[inds] <- 1
  ## turns the occurrence matrix in vector as returning it
  c(M)
}
