# test-integration-sweep-matrix.R: Projeção de matrizes genéricas e scRNA-seq (matrix, array, integer, dgCMatrix)

test_that("SWeeP e SWeePlite projetam matrizes genéricas e esparsas com total paridade", {
  # Carrega matriz esparsa de teste e projeções de referência
  sp_mat <- as(load_golden_mtx("sweep_matrix_sparse_10x25.mtx"), "CsparseMatrix")
  dense_mat <- as.matrix(sp_mat)
  golden_proj <- load_golden_csv("sweep_matrix_dense_10x25_all_modes.csv")

  psz <- 5
  base <- orthBase(lin = 25, col = psz, seed = 647474747)

  # 1. Projeção clássica densa vs esparsa em count_none
  proj_dense <- SWeeP(dense_mat, orthbase = base, bin = FALSE, norm = "none", verbose = FALSE)$proj
  proj_sparse <- SWeeP(sp_mat, orthbase = base, bin = FALSE, norm = "none", verbose = FALSE)$proj
  expect_equal(proj_dense, as.matrix(proj_sparse), tolerance = 1e-12)
  expect_equal(proj_dense, unname(golden_proj[, 1:5]), tolerance = 1e-12)

  # 2. Equivalência SWeeP == SWeePlite na matriz esparsa
  proj_lite_sparse <- SWeePlite(sp_mat, psz = psz, bin = FALSE, norm = "none", ncores = 1, verbose = FALSE)$proj
  expect_equal(as.matrix(proj_sparse), proj_lite_sparse, tolerance = 1e-12)

  # 3. Modo bin_none
  proj_bin_sweep <- SWeeP(dense_mat, orthbase = base, bin = TRUE, norm = "none", verbose = FALSE)$proj
  proj_bin_lite <- SWeePlite(dense_mat, psz = psz, bin = TRUE, norm = "none", ncores = 1, verbose = FALSE)$proj
  expect_equal(proj_bin_sweep, proj_bin_lite, tolerance = 1e-12)
  expect_equal(proj_bin_sweep, unname(golden_proj[, 6:10]), tolerance = 1e-12)

  # 4. Modo count_log
  proj_log_sweep <- SWeeP(dense_mat, orthbase = base, bin = FALSE, norm = "log", verbose = FALSE)$proj
  proj_log_lite <- SWeePlite(dense_mat, psz = psz, bin = FALSE, norm = "log", ncores = 1, verbose = FALSE)$proj
  expect_equal(proj_log_sweep, proj_log_lite, tolerance = 1e-12)
  expect_equal(proj_log_sweep, unname(golden_proj[, 11:15]), tolerance = 1e-12)

  # 5. Modo bin_logneg
  proj_logneg_sweep <- SWeeP(dense_mat, orthbase = base, bin = TRUE, norm = "logNeg", verbose = FALSE)$proj
  proj_logneg_lite <- SWeePlite(dense_mat, psz = psz, bin = TRUE, norm = "logNeg", ncores = 1, verbose = FALSE)$proj
  expect_equal(proj_logneg_sweep, proj_logneg_lite, tolerance = 1e-12)
  expect_equal(proj_logneg_sweep, unname(golden_proj[, 16:20]), tolerance = 1e-12)
})

test_that("Transposição e RNAseqdata operam de forma equivalente", {
  mat <- matrix(1:20, nrow = 4, ncol = 5)
  base <- orthBase(lin = 4, col = 3, seed = 647474747)

  # Se transpose = TRUE, cada coluna é uma amostra (dimensão 5 x 3)
  res_trans <- SWeeP(mat, orthbase = base, transpose = TRUE, verbose = FALSE)$proj
  res_rna   <- SWeeP(mat, orthbase = base, RNAseqdata = TRUE, verbose = FALSE)$proj

  expect_equal(dim(res_trans), c(5, 3))
  expect_equal(res_trans, res_rna)

  # Transposição explícita antes
  mat_t <- t(mat)
  res_direct <- SWeeP(mat_t, orthbase = base, transpose = FALSE, verbose = FALSE)$proj
  expect_equal(res_trans, res_direct)
})

test_that("Linearidade algébrica da projeção randômica (Invariante Matemático)", {
  # Com norm='none' e bin=FALSE, SWeeP(a*X + b*Y) == a*SWeeP(X) + b*SWeeP(Y)
  set.seed(123)
  X <- matrix(stats::runif(30), nrow = 6, ncol = 5)
  Y <- matrix(stats::runif(30), nrow = 6, ncol = 5)
  base <- orthBase(lin = 5, col = 3, seed = 647474747)

  a <- 2.5
  b <- -1.3

  proj_comb <- SWeeP(a * X + b * Y, orthbase = base, bin = FALSE, norm = "none", verbose = FALSE)$proj
  proj_X    <- SWeeP(X, orthbase = base, bin = FALSE, norm = "none", verbose = FALSE)$proj
  proj_Y    <- SWeeP(Y, orthbase = base, bin = FALSE, norm = "none", verbose = FALSE)$proj

  expect_equal(proj_comb, a * proj_X + b * proj_Y, tolerance = 1e-12)
})

test_that("Despacho S4 suporta array e integer", {
  int_mat <- matrix(1L:15L, nrow = 3, ncol = 5)
  storage.mode(int_mat) <- "integer"
  base <- orthBase(lin = 5, col = 2, seed = 647474747)

  res_int <- SWeeP(int_mat, orthbase = base, verbose = FALSE)$proj
  expect_equal(dim(res_int), c(3, 2))

  arr <- as.array(int_mat)
  res_arr <- SWeeP(arr, orthbase = base, verbose = FALSE)$proj
  expect_equal(res_int, res_arr)
})
