aa2num2 <- function(xseq){
  b <- dim(xseq);
  m <- b[2]
  n <- b[1]
  xseq <- matrix(as.numeric(aa2int(xseq)), ncol = length(xseq[1, ]))
  vls<-matrix(as.double(round(xseq, digits=4))-1, ncol = length(xseq[1, ]))
  pot<-kronecker(matrix(1,n,1),matrix(0:(m-1),ncol=2))
  premat <-(kronecker(matrix(1,n,m),20)^pot)*vls
  mret <- matrix(premat[ ,1]+premat[ , 2], ncol=1) +1
  pinot<-(1-matrix(as.numeric(vls<0 | vls>19), ncol= length(xseq[1, ])))
  inot<- which(matrix(as.numeric(!(pinot[ ,1]*pinot[ ,2])), ncol = 1)!=0)
  mret[inot] = -1
  mret
}
