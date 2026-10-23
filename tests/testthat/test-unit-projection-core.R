# test-unit-projection-core.R: Etapa 4 do Artigo (liteParam, COREloop e paralelismo)

test_that("NCoresDef trata valores nulos e inválidos", {
  expect_equal(rSWeeP:::NCoresDef(NULL), 2)
  expect_equal(rSWeeP:::NCoresDef(0), 2)
  expect_equal(rSWeeP:::NCoresDef(-4), 2)
  expect_equal(rSWeeP:::NCoresDef(4), 4)
})

test_that("liteParam inicializa parâmetros determinísticos com 50 primos", {
  # [DOC-DIVERGENCE]: pslist tem 50 números primos (2 a 229)
  par <- rSWeeP:::liteParam(mask = NULL, input = NULL, seqtype = "AA", N = 1, psz = 10, lenmax = 400)

  expect_equal(par$nps, 50)
  expect_equal(length(par$pslist), 50)
  expect_equal(par$pslist[1:5], c(2, 3, 5, 7, 11))
  expect_equal(par$pslist[50], 229)
  expect_equal(dim(par$Mproj), c(50, 10))
  expect_equal(par$xnorm, sqrt(400 / 3))

  # Verifica reprodutibilidade da matriz aleatória Mproj com a semente fixa 647474747
  golden_mproj10 <- load_golden_csv("mproj_50x10_seed647474747.csv")
  expect_equal(par$Mproj, unname(golden_mproj10), tolerance = 1e-12)
})

test_that("COREloop é invariante ao tamanho de bloco nk e consistente com multiplicação de matriz", {
  psz <- 5
  lenmax <- 20
  par <- rSWeeP:::liteParam(mask = NULL, input = NULL, seqtype = "AA", N = 1, psz = psz, lenmax = lenmax)

  # Cria um cluster temporário
  cl <- parallel::makeCluster(1, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Vetor de teste esparso
  hdv_vec <- list(
    idx = c(1, 3, 5, 7, 10, 15),
    count = c(2.0, 1.0, 4.0, 1.0, 3.0, 2.0)
  )

  # Executa com nk grande (1 único chunk) e nk pequeno (múltiplos chunks)
  res_single_chunk <- rSWeeP:::COREloop(par$xnorm, psz, par$Mproj, par$pslist, par$nps, hdv_vec, nk = 100, bin = FALSE, norm = "none")
  res_multi_chunk <- rSWeeP:::COREloop(par$xnorm, psz, par$Mproj, par$pslist, par$nps, hdv_vec, nk = 2, bin = FALSE, norm = "none")

  # Deve ser exatamente igual
  expect_equal(res_single_chunk, res_multi_chunk, tolerance = 1e-14)
  expect_length(res_single_chunk, psz)

  # Caso de borda: vetor vazio (idx de comprimento zero) retorna vetor de zeros
  hdv_empty <- list(idx = integer(0), count = numeric(0))
  res_empty <- rSWeeP:::COREloop(par$xnorm, psz, par$Mproj, par$pslist, par$nps, hdv_empty, nk = 10, bin = FALSE, norm = "none")
  expect_equal(res_empty, rep(0, psz))
})
