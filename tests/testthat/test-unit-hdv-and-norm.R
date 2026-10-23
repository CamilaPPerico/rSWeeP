# test-unit-hdv-and-norm.R: Etapa 3 do Artigo (Construção e Normalização do HDV)

test_that("readmask e readmaskVEC são consistentes em todos os modos bin e norm", {
  # Sequência sintética gerando k-mers com contagens controladas
  seq_str <- "ACDACDACD"
  seq_num <- rSWeeP:::seq2num(seq_str, "AA")
  mask <- c(1, 0, 1) # k=2 ativos, lenmax = 20^2 = 400
  lenmax <- 400

  # 1. bin=FALSE, norm='none'
  r_vec_count <- rSWeeP:::readmaskVEC(seq_num, mask, "AA", bin = FALSE, norm = "none", lenmax = lenmax)
  r_list_count <- rSWeeP:::readmask(seq_num, mask, "AA", bin = FALSE, norm = "none", lenmax = lenmax)
  expect_equal(which(r_vec_count != 0), r_list_count$idx)
  expect_equal(r_vec_count[r_list_count$idx], r_list_count$count)

  # 2. bin=TRUE, norm='none'
  r_vec_bin <- rSWeeP:::readmaskVEC(seq_num, mask, "AA", bin = TRUE, norm = "none", lenmax = lenmax)
  r_list_bin <- rSWeeP:::readmask(seq_num, mask, "AA", bin = TRUE, norm = "none", lenmax = lenmax)
  expect_equal(which(r_vec_bin != 0), r_list_bin$idx)
  expect_true(all(r_vec_bin[r_list_bin$idx] == 1))
  expect_null(r_list_bin$count) # em bin=TRUE e norm='none', readmask omite count

  # 3. bin=FALSE, norm='log' (log na base 10)
  r_vec_log <- rSWeeP:::readmaskVEC(seq_num, mask, "AA", bin = FALSE, norm = "log", lenmax = lenmax)
  r_list_log <- rSWeeP:::readmask(seq_num, mask, "AA", bin = FALSE, norm = "log", lenmax = lenmax)
  expect_equal(which(r_vec_log != 0), r_list_log$idx)
  expect_equal(r_vec_log[r_list_log$idx], r_list_log$count)

  # 4. bin=FALSE, norm='Neg' (zeros viram -1)
  r_vec_neg <- rSWeeP:::readmaskVEC(seq_num, mask, "AA", bin = FALSE, norm = "Neg", lenmax = lenmax)
  r_list_neg <- rSWeeP:::readmask(seq_num, mask, "AA", bin = FALSE, norm = "Neg", lenmax = lenmax)
  expect_equal(r_vec_neg, r_list_neg$count)
  expect_true(all(r_vec_neg[r_vec_count == 0] == -1))

  # 5. bin=FALSE, norm='logNeg' (zeros viram -1 e contagens viram log10)
  r_vec_logneg <- rSWeeP:::readmaskVEC(seq_num, mask, "AA", bin = FALSE, norm = "logNeg", lenmax = lenmax)
  r_list_logneg <- rSWeeP:::readmask(seq_num, mask, "AA", bin = FALSE, norm = "logNeg", lenmax = lenmax)
  expect_equal(r_vec_logneg, r_list_logneg$count)
  expect_true(all(r_vec_logneg[r_vec_count == 0] == -1))
})

test_that("MakeLOG transforma matrizes corretamente para todos os 8 modos", {
  mat <- matrix(c(0, 1, 10, 100), nrow = 2, byrow = TRUE)

  # bin=FALSE, norm='none' -> inalterado
  expect_equal(rSWeeP:::MakeLOG(mat, bin = FALSE, norm = "none"), mat)

  # bin=TRUE, norm='none' -> binarização (não nulos viram 1)
  expect_equal(rSWeeP:::MakeLOG(mat, bin = TRUE, norm = "none"),
               matrix(c(0, 1, 1, 1), nrow = 2, byrow = TRUE))

  # bin=FALSE, norm='log' -> log10 sobre não nulos
  expect_equal(rSWeeP:::MakeLOG(mat, bin = FALSE, norm = "log"),
               matrix(c(0, 0, 1, 2), nrow = 2, byrow = TRUE))

  # bin=FALSE, norm='Neg' -> zeros viram -1
  expect_equal(rSWeeP:::MakeLOG(mat, bin = FALSE, norm = "Neg"),
               matrix(c(-1, 1, 10, 100), nrow = 2, byrow = TRUE))

  # bin=FALSE, norm='logNeg' -> zeros viram -1 e log10 nos demais
  expect_equal(rSWeeP:::MakeLOG(mat, bin = FALSE, norm = "logNeg"),
               matrix(c(-1, 0, 1, 2), nrow = 2, byrow = TRUE))

  # bin=TRUE, norm='logNeg' -> zeros viram -1 e não zeros viram 1
  expect_equal(rSWeeP:::MakeLOG(mat, bin = TRUE, norm = "logNeg"),
               matrix(c(-1, 1, 1, 1), nrow = 2, byrow = TRUE))
})
