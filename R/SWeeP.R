#' @title A vectorial comparative method to amino acid sequence.
#' @name SWeeP
#' @description The "Spaced Words Projection (SWeeP)" is a method for
#' representing biological sequences using relativety, it uses the
#' spacedwords concept by scanning sequences and generating indices
#' to create a higherdimensional vector that is later projected into a
#' smaller randomly oriented orthonormal base (PIERRI, 2019)
#' @param file A FASTA format file
#' @param orthBase A orthonormal matrix with 160.000 rows
#' @details The SWeeP method was developed to favor the comparison between
#' complete genomic sequences and to assist in phylogenetic analyzes. This
#' method is based on the concept of sparse words, which are used to scan
#' biological sequences and turn them into vectors, favoring the manipulation
#' of the data. It can project a matrix n by m, where n is the number of
#' sequences in the analized \code{file} and m is the number of columns
#' in \code{orthBase} matrix.
#' @return A matrix resulted by the projection of the sequences in \code{file}
#' in the \code{orthBase} matrix
#' @author Roberto Raittz, et. al.
#' @seealso \code{\link[rSWeeP]{orthbase}}
#' @references Pierri,C. R. et al. SWeeP: Representing large biological
#' sequences data sets in compact vectors. (2019).
#' @examples
#' baseMatrix <- orthbase(160000,10)
#' data(datastring)
#' cat(datastring,file = "exdna.fas", sep = "\n")
#' return <- SWeeP("exdna.fas",baseMatrix)
#' distancia <- dist(return, method = "euclidean")
#' tree <- hclust(distancia, method="ward.D")
#' plot(tree, hang = -1, cex = 1)
#' @export
SWeeP <- function(file,orthBase){
  xfas <- ape::read.dna(file, format = "fasta",as.character = TRUE)
  if(is.list(xfas)){
    if (length(xfas) > 2001){
      barra<-round(length(xfas)/2000)
      Bigprojected <- matrix(2,nrow=length(xfas) ,ncol=length(orthBase))
      blocksToConvert <- length(xfas)/barra
      block <- round(length(xfas)/round(length(xfas)/blocksToConvert))
      p<-1
      for (i in c(1:(round(length(xfas)/block)))){
        print(paste0("Executando parte: ",i, ", de: ",barra))
        if (p+block<length(xfas)){
          projected<-fas2mat(xfas[p:(p+block)])%*%as.matrix(orthBase)
          Bigprojected[p:(p+block),]<-projected
          p <- p+block+1
        }else
        {
          projected<-fas2mat(xfas[p:length(xfas)])%*%as.matrix(orthBase)
          Bigprojected[p:length(xfas),]<-projected
        }
        projected <-Bigprojected
      }
    }else
    {
      matrizid <-fas2mat(xfas)
      projected <- list();
      projected <- (matrizid%*%as.matrix(orthBase))
    }
    end_time <- Sys.time()
    projected
  }else
  {
    print("The amino acid sequences are too short ")
  }
}
