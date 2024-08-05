\name{rSWeeP-package}
\alias{rSWeeP-package}
\alias{rSWeeP}
\docType{package}
\title{
  rSWeeP -- Alignment-free method for vectorising biological sequences
}
\description{
  The ´rSWeeP´ package is an R implementation of the Spaced Words Projection (SWeeP) method (De Pierri, 2019). The main function of this package is to provide a vector representation of biological sequences (nucleotides or amino acids), and thus favor alignment-free phylogenetic studies. Each sequence provided is represented by a compact numerical vector which is easier to analyze. SWeeP uses k-mers counting for representing the sequences in high dimensional vector (HDV) and then projected into a low dimensional vector (LDV) through random projection using an orthonormal base. The LDV represents the biological sequence and is handable for comparative analisys and machine learning implements.

  In addition, the package allows general dimensionality reduction of RNAseq data and matrices. 
  
  More information about 'rSWeeP' can be found at 
  <https://github.com/CamilaPPerico/rSWeeP>
  Tutorials on package use can be found at
  <https://aibialab.github.io/rSWeeP>

}
\details{
  % This section should provide a more detailed overview of how to use the
  % package, including the most important functions.
  The main functions are SWeeP, SWeePlite and orthBase. Additionally are available extractHDV 
}
\author{
Maintainer: Camila P Perico <camilapp94@gmail.com>

Camila P. Perico [cre, aut, cph] 
Roberto T. Raittz  [aut, cph]
Danrley R. Fernandes  [aut]
Mariane G. Kulik [aut]
}
\references{
  % This optional section can contain literature or other references for
  % background information.
  De Pierri, C. R., Voyceik, R., Santos de Mattos, L. G. C., Kulik, M. G., Camargo, J. O., Repula de Oliveira, A. M., ... & Raittz, R. T. (2020). Sweep: representing large biological sequences datasets in compact vectors. Scientific reports, 10(1), 91.

}
\keyword{ aligment-free, dimension reduction, SWeeP, sequence vectorization, RNA sequencing }
\examples{
  \dontrun{
    
    path <- paste (system.file("examples/aaMitochondrial/",package = "rSWeeP"),'/', sep = '')
    sw = SWeePlite(path,mask=c(1,1,0,1,1),psz=1500,seqtype='AA',ncores=4)

    pathmetadata <- system.file(package = "rSWeeP" , "examples" , "metadata_mitochondrial.csv")
    mt = read.csv(pathmetadata,header=TRUE)
    pca_output <- prcomp (sw$proj , scale = FALSE)
    plot(pca_output$x[,1],pca_output$x[,2] , xlab = 'PC-1' , ylab = 'PC-2', pch=20, col=mt$id)
    legend("bottomright",unique(mt$family),col=as.character(c(1:length(unique(mt$family)))),pch=20)
  
  }
}
