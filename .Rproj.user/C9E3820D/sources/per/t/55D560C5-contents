dna2list <- function(sdna, colunSize) {
  ## geting the sizes of the sequence given
  size <- length(sdna)
  ## create a column of integers from 1 to sequence length
  xcol <- matrix(1:(size - colunSize + 1), ncol = 1)
  ## creates a matrix with the column Size of length and the second number in each line being the first
  ## number in the next one until the sequence length be reached
  inds <- kronecker(matrix(1, 1, colunSize), xcol) + kronecker(matrix(1, size - colunSize + 1, 1), matrix(1:colunSize,
  ncol = colunSize)) - 1
  ## uses inds as index to create a matrix of colunSize with the sdna sequence
  mret = matrix(sdna[inds], ncol = colunSize)
  ## returns matrix
  mret
}
