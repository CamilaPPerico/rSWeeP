dna2list <- function(sdna, q){
  n <- length(sdna)
  xcol <- matrix(1:(n-q+1),ncol =1) 
  inds <- kronecker(matrix(1,1,q),xcol) + kronecker(matrix(1,n-q+1,1),matrix(1:q,ncol = q)) -1
  mret = matrix(sdna[inds],ncol = q)
  list(mret)
}
