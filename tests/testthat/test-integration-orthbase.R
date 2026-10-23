# test-integration-orthbase.R: Geração da base de projeção e validação de invariantes matemáticos

test_that("orthBase gera base determinística consistente com os Golden Files", {
  # Base para AA (400 x 10)
  base_aa <- orthBase(lin = 400, col = 10, seed = 647474747)
  expect_equal(dim(base_aa$mat), c(400, 10))
  expect_equal(base_aa$seed, 647474747)
  expect_equal(base_aa$version, "SWeeP v2.9")

  golden_aa <- load_golden_csv("orthbase_aa_400x10.csv")
  expect_equal(base_aa$mat, unname(golden_aa), tolerance = 1e-12)

  # Base para NT (256 x 10)
  base_nt <- orthBase(lin = 256, col = 10, seed = 647474747)
  expect_equal(dim(base_nt$mat), c(256, 10))
  golden_nt <- load_golden_csv("orthbase_nt_256x10.csv")
  expect_equal(base_nt$mat, unname(golden_nt), tolerance = 1e-12)
})

test_that("orthBase cumpre propriedades de quase-ortogonalidade (Invariante do Artigo)", {
  # Conforme o artigo (De Pierri et al., 2020), as colunas de orthBase possuem norma aproximadamente unitária
  # e produtos internos entre colunas distintas próximos de zero
  base <- orthBase(lin = 400, col = 10, seed = 647474747)$mat

  # 1. Normas das colunas aproximadamente 1
  col_norms_sq <- colSums(base^2)
  expect_true(all(col_norms_sq > 0.8 & col_norms_sq < 1.2))

  # 2. Quase-ortogonalidade: B^T B tem diagonal ~ 1 e fora da diagonal próximo de zero
  gram <- t(base) %*% base
  diag_vals <- diag(gram)
  diag(gram) <- 0
  max_off_diag <- max(abs(gram))

  expect_true(all(diag_vals > 0.8 & diag_vals < 1.2))
  expect_true(max_off_diag < 0.3)
})

test_that("orthBase aceita mask diretamente sem especificar lin", {
  # mask c(1,0,1) para AA -> lin = 20^2 = 400
  base1 <- orthBase(mask = c(1, 0, 1), col = 10, seqtype = "AA", seed = 647474747)
  base2 <- orthBase(lin = 400, col = 10, seed = 647474747)
  expect_equal(base1$mat, base2$mat)

  # Default seed é 647474747
  base_default_seed <- orthBase(lin = 400, col = 10)
  expect_equal(base_default_seed$mat, base2$mat)
})
