#' @title Generate a orthonormal matrix (\code{lin} , \code{col})
#' @name orthBase
#'
#' @description Generate a orthonormal matrix in a specified
#' size, \code{lin} by \code{col}.
#'
#' @param lin Number of rows in the desired matrix
#' @param col Number of columns in the desired matrix
#'
#' @return A orthonormal matrix in a specified size, \code{lin} by \code{col}.
#'
#' @author Danrley R. Fernandes.
#'
#' @seealso \code{\link[rSWeeP]{sWeeP}}, \code{\link[pracma]{orth}}
#' @examples
#' orthBase(160000, 10)
#'
#' lin <- 160000
#' col <- 10
#' orthBase(lin = lin, col = col)
#'
#'
#' @export

orthBase <- function(lin, col) {
  Ro <- pracma::orth(matrix(stats::runif((col + 1) * lin), lin))

  mret <- Ro[, -1]
}
