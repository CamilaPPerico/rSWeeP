\name{rSWeeP-package}
\alias{rSWeeP-package}
\alias{rSWeeP}
\docType{package}
\title{
  Spaced Words Projection (SWeeP) -- Alignment-free method for vectorising biological sequences
}
\description{
  The "Spaced Words Projection (SWeeP)" is a method for
  representing biological sequences using compact vectors. SWeeP uses k-mers for representing the sequences
  in high dimensional vector (HDV) and then projected into a low dimensional vector (LDV) through a 
  orthonormal base (De Pierri, 2019). The LDV represents the biological sequence and is handable 
  for comparative analisys and machine learning implements.
  
  ´rSWeeP´ fornece funções para vetorização de sequências biológicas

  ´rSWeeP´ provides a vectorial comparative method to biological sequences (amino acid and nucleotide), 
  expression matrices and generic matrices, with dimension reduction and low loss of information.

  Mais informações sobre ´rSWeeP´ podem ser encontradas em
  github
  A more detailed description of what the package does. A length
  of about one to five lines is recommended.


  % Description:

     % ‘ape’ provides functions for reading, writing, manipulating,
     % analysing, and simulating phylogenetic trees and DNA sequences,
     % computing DNA distances, translating into AA sequences, estimating
     % trees with distance-based methods, and a range of methods for
     % comparative analyses and analysis of diversification.
     % Functionalities are also provided for programming new phylogenetic
     % methods.

     % The complete list of functions can be displayed with ‘library(help
     % = ape)’.

     % More information on ‘ape’ can be found at
     % <http://ape-package.ird.fr/>.

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
    # mt=readRDS(pathmetadata)
    vunq = c('Hominidae','Psittacidae','Muridae','Bovidae')
    cl = mt$id
    ncl = length(vunq)

    pca_output <- prcomp (sw$proj , scale = FALSE)
    plot(pca_output$x[,1] , pca_output$x [,2] , xlab = 'PC-1' , ylab = 'PC-2' , pch =20 , col = cl)
    legend("topleft" , vunq , col = as.character(c(1:ncl)) , pch =20)


  
  }
}
