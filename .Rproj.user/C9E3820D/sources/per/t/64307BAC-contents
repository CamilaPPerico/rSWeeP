aa2num2 <- function(xseq) {

  map <- c(A = 1L, a = 1L, R = 2L, r = 2L, n = 3L, N = 3L, d = 4L, D = 4L, c = 5L, C = 5L, q = 6L, Q = 6L,
           e = 7L, E = 7L, g = 8L, G = 8L, h = 9L, H = 9L, i = 10L, I = 10L, l = 11L, L = 11L, k = 12L, K = 12L,
           m = 13L, M = 13L, f = 14L, F = 14L, p = 15L, P = 15L, s = 16L, S = 16L, t = 17L, T = 17L, w = 18L,
           W = 18L, y = 19L, Y = 19L, v = 20L, V = 20L, b = 21L, B = 21L, z = 22L, Z = 22L, x = 23L, X = 23L,
           `*` = 24L, `-` = 25L)
  # geting the sizes of the sequence given
  xseqSize <- dim(xseq)
  ## converting the characters in numeric
  xseq <- matrix(map[xseq], ncol = length(xseq[1, ]))
  ## recovers a matrix of index positions of each amino acids in the sequence
  vls <- matrix(xseq - 1, ncol = length(xseq[1, ]))
  pot <- kronecker(matrix(1, xseqSize[1], 1), matrix(0:(xseqSize[2] - 1), ncol = 2))
  premat <- (matrix(20, xseqSize[1], xseqSize[2])^pot) * vls
  mret <- matrix(premat[, 1] + premat[, 2], ncol = 1) + 1
  ## removes all invalids amino acids from the index
  pinot <- (1 - matrix(as.numeric(vls < 0 | vls > 19), ncol = length(xseq[1, ])))
  inot <- which(matrix(as.numeric(!(pinot[, 1] * pinot[, 2])), ncol = 1) != 0)
  mret[inot] <- -1
  ## returns the index matrix
  mret
}
