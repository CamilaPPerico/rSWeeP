#' @title A vectorial comparative method to amino acid sequence.
#' @name sWeeP
#' @description The "Spaced Words Projection (SWeeP)" is a method for
#' representing biological sequences using compact vectors. SWeeP uses the
#' spacedwords concept by scanning sequences and generating indices
#' to create a higherdimensional matrix of occurrences that is later projected into a
#' smaller randomly oriented orthonormal base (PIERRI, 2019). This way the resulting
#' matrix  will conserve the comparational data but will have a selectable size
#' @param xfas A AAStringSet or a FASTA format file
#' @param baseMatrix A orthonormal matrix with 160.000 coordinates
#' @details The SWeeP method was developed to favor the comparison between
#' complete proteomic sequences and to assist in machine learning analyzes. This
#' method is based on the concept of spaced words, which are used to scan
#' biological sequences and project them into  matrix of occurrences, favoring the manipulation
#' of the data. The sWeeP function can project a matrix n by m, where n is the number of
#' sequences in the analized \code{xfas} and m is the number of columns
#' in \code{baseMatrix} matrix.
#' @return A matrix resulted by the projection of the sequences in \code{xfas}
#' in the \code{baseMatrix} matrix
#' @author Danrley R. Fernandes.
#' @references Pierri,C. R. et al. SWeeP: Representing large biological
#' sequences data sets in compact vectors. Scientific Reports, accepted in
#' December 2019.doi: 10.1038/s41598-019-55627-4.
#' @examples
#' baseMatrix <- orthBase(160000,10)
#' path <- system.file(package = "rSWeeP", "extdata", "exdna.fas")
#' return <- sWeeP(path,baseMatrix)
#' distancia <- dist(return, method = "euclidean")
#' tree <- hclust(distancia, method="ward.D")
#' plot(tree, hang = -1, cex = 1)
#' @rdname sWeeP
#' @export

##assuring that exposed functions of your package interoperate with standard Bioconductor objects by an S4 Generic and S4 methods
setGeneric(
  "sWeeP",
  function(xfas, baseMatrix) standardGeneric("sWeeP"),
  signature = "xfas"
)
#' @rdname sWeeP
setMethod("sWeeP", "character", function(xfas, baseMatrix) {
  aa = Biostrings::readAAStringSet(xfas,"fasta", use.names = FALSE)
  sWeeP(aa, baseMatrix)
})
#' @rdname sWeeP
setMethod("sWeeP", "AAStringSet", function(xfas, baseMatrix) {
  ##validating user input.
  stopifnot(methods::is(baseMatrix, "matrix"))
  print('starting')
  ##if the data sequence has more than than 2000 specimens it's divide to process in pieces of 2000 at time
  if (length(xfas) > 2001) {
    ##define the number of fragmentation needed
    barra <- round(length(xfas)/2000)
    blocksToConvert <- length(xfas)/barra
    block <- round(length(xfas)/round(length(xfas)/blocksToConvert))
    ##create the matrix who will contain the projection
    Bigprojected <- matrix(0,0,dim(baseMatrix)[2])
      ##position defines the start of projection a at position 1
     position <- 1
      ##runs the projection thru all the matrix
    for (i in  1:barra) {
      message("running ", i, " of ", barra)
      ##verifi if the next matrix area isn't outside limits
      if (position + block < length(xfas)) {
        projected <- fas2mat(xfas[position:(position + block)]) %*% baseMatrix
        Bigprojected<-rbind(Bigprojected,projected)
        position <- position + block + 1
      }
      ##if the next matrix area is outside limits performs the projection only until end of the matrix
      else {
        projected <- fas2mat(xfas[position:length(xfas)]) %*% baseMatrix
        Bigprojected<-rbind(Bigprojected,projected)
      }
      projected <- Bigprojected
    }
  }
  ##if the data sequence has less than than 2000 specimens it's run the projection directly
  else {
    matrizid <- fas2mat(xfas)
    projected <- matrizid %*% baseMatrix
  }
  ##returns the projection matrix
  projected
})
