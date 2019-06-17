fas2mat <- function(xfas){
  S <-list()
  S<-seqinr::getSequence(xfas)
  S<-vapply(S,aa2mat,FUN.VALUE =double(160000))
  t(S)
}
aa2mat <- function(xseq){

  nucol<- 5
  L5 <- dna2list(xseq,nucol);
  L5 <-L5[[1]]
  xy <- cbind(aa2num2(L5[ ,1:2]), aa2num2(L5[ ,4:nucol]))

  pinot<-(1-matrix(xy<=0,ncol=2))
  inot<- which(matrix(as.numeric(!(pinot[ ,1]*pinot[ ,2])), ncol = 1)!=0)
  xy[inot, ] <- NA
  inds <- ij2inds(xy,400)
  inds[inds>(160000)] <- NA
  M <-matrix(0,ncol= 400,nrow = 400)
  M[inds]<-1
  c(M)

}
