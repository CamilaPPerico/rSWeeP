# test-integration-extracthdv.R: Extração do Vetor de Alta Dimensão (HDV)

test_that("extractHDV opera em memória sobre AAStringSet e DNAStringSet", {
  aa_seqs <- Biostrings::AAStringSet(c(
    seq1 = "ACDEFGHIKLMNPQRSTVWY",
    seq2 = "ACDACDACDACD"
  ))

  # mask c(1,0,1) -> 400 colunas
  res <- extractHDV(aa_seqs, mask = c(1, 0, 1), seqtype = "AA", bin = FALSE, verbose = FALSE)

  expect_type(res, "list")
  expect_equal(dim(res$HDV), c(2, 400))
  expect_equal(res$info$headers, c("seq1", "seq2"))
  expect_equal(res$info$mask, c(1, 0, 1))
  expect_equal(res$info$SequenceType, "AA")
  expect_false(res$info$concatenate)
  expect_equal(length(res$info$saturation), 2)
  expect_true(all(res$info$saturation > 0 & res$info$saturation <= 1))

  # Teste em modo binário
  res_bin <- extractHDV(aa_seqs, mask = c(1, 0, 1), seqtype = "AA", bin = TRUE, verbose = FALSE)
  expect_true(all(res_bin$HDV %in% c(0, 1)))

  # Teste com concatenação
  res_cat <- extractHDV(aa_seqs, mask = c(1, 0, 1), seqtype = "AA", concatenate = TRUE, verbose = FALSE)
  expect_equal(dim(res_cat$HDV), c(1, 400))
  expect_true(res_cat$info$concatenate)
})

test_that("extractHDV a partir de pasta é consistente com os Golden Files", {
  path_mito <- system.file("examples", "aaMitochondrial", package = "rSWeeP")
  if (path_mito == "") path_mito <- "/home/rasert/workspace/rSWeeP/inst/examples/aaMitochondrial"

  # Máscara compacta c(1,1,1) -> 400 posições
  res_compact <- extractHDV(path_mito, mask = c(1, 1, 1), seqtype = "AA", bin = FALSE, extension = ".faa", verbose = FALSE)
  golden_compact <- load_golden_csv("hdv_compact_aa_400.csv")
  expect_equal(res_compact$HDV, unname(golden_compact), tolerance = 1e-12)

  # Máscara canônica c(2,1,2) -> 160.000 posições (Golden MTX esparso)
  golden_mtx <- load_golden_mtx("hdv_mito_14x160000.mtx")
  expect_equal(dim(golden_mtx), c(14, 160000))
  # Saturação de proteomas mitocondriais é tipicamente baixa (~1-5%)
  saturation <- Matrix::rowSums(golden_mtx != 0) / 160000
  expect_true(all(saturation > 0.005 & saturation < 0.10))
})
