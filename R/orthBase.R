#' @title Generate a orthonormal matrix (\code{lin} , \code{col})
#' @name orthBase
#'
#' @description Generate a orthonormal matrix in a specified
#' size, \code{lin} by \code{col}.
#'
#' @param lin Number of rows in the desired matrix
#' @param col Number of columns in the desired matrix
#' @param mask      reading mask. Use this option or `lin' option. Default c(2,1,2).
#' @param seqtype   type of data: AA for amino acid, NT for nucleotide. Parameter required if a mask is provided. The default is AA
#' @param seed   provide, if necessary, a seed to generate the matrix. The default is 647474747
#'
#' @return An orthonormal matrix (basis) whose dimensions correspond to the given mask
#'         to be used and a desired projection size (length of the output vector). 
#'         The basis must be supplied to the function \link{SWeeP} (see examples).
#'
#'         `orthBase' returns a `list` containing:
#' \itemize{
#'          \item mat: the orthonormal matrix (basis)
#'          \item seed: the random seed (metadata to identify the matrix)
#'          \item version: the rSWeeP version
#'         }
#'
#' @author Camila P. Perico
#'
#' @seealso \code{\link{SWeeP}}
#' @examples
#' 
#' # define the mask - determines the length of input vector (20^4 = 160000)
#' mask <- c(1,1,0,1,1) 
#' 
#' # define the length of output vector
#' col <- 600
#' 
#' # get the basis matrix to projection
#' Mybase <- orthBase(mask = mask, col = col,seqtype='AA')
#' 
#' @import methods
#' @export
orthBase <- function(lin=NULL, col,seqtype='AA',mask = c(2,1,2),seed=NULL) {

    if (max(mask)>1){
        mask = convertMask(mask)
    }

    # COLOCAR PSEUDORANDOMICO - sempre a mesma!
    if(is.null(seed)){
        seed = 647474747
        set.seed(seed) # fixed
    } else{ 
        set.seed(seed)
    }

    if( length(lin) == 0 ){    # if the user give the mask, not the number of lines

        SW.checks('mask',mask)

        if (seqtype=='AA'){
            lin = 20^sum(mask)
        } else if(seqtype=='NT'){
            lin = 4^sum(mask)
        }
    }



    idx = 1:lin
    pslist = c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229)
    nps = length(pslist)
    Mproj =  matrix(stats::runif(nps * col), ncol=col)
    xnorm = sqrt(lin/3)

    mret = (matrix(rep(idx,nps),ncol=nps))%%(matrix(rep(pslist,length(idx)),ncol=nps,byrow=T))
    dt = (1+(mret))%*%Mproj
    bs = ((dt - floor(dt))-0.5)/0.5

    output=NULL
    output$mat = bs/xnorm
    output$seed = seed
    output$version = 'SWeeP v2.9'

    return(output)
}
