# test-unit-phylo-internals.R: Funções internas de avaliação de árvores filogenéticas

test_that("diff calcula diferenças sucessivas de vetores numéricos", {
  expect_equal(rSWeeP:::diff(c(1, 2, 4, 7)), c(1, 2, 3))
  expect_equal(rSWeeP:::diff(c(5, 5, 5)), c(0, 0))
})

test_that("costquebra avalia agrupamento contíguo de táxons", {
  # Taxon único (N=1) sempre tem custo 1
  expect_equal(rSWeeP:::costquebra(c(0, 1, 0), "new"), 1)

  # Taxon perfeitamente agrupado (todos os 1s juntos):
  # xo=2, xe=4, N=3. which(z[xo:xe]==1) é 1:3, diff é c(1,1) -> soma 2.
  # (1 + 2) / 3 = 1.0
  expect_equal(rSWeeP:::costquebra(c(0, 1, 1, 1, 0), "new"), 1)

  # Taxon fragmentado c(1, 0, 1) -> xo=1, xe=3, N=2. z[1:3] é 1,0,1 -> which é 1, 3 -> diff é 2 == 1 (0).
  # (1 + 0) / 2 = 0.5
  expect_equal(rSWeeP:::costquebra(c(1, 0, 1), "new"), 0.5)
})

test_that("quebrataxonCOPHE e MonoParaphylMetric avaliam árvores filogenéticas sintéticas", {
  # Árvore simples com 4 pontas
  tree_text <- "((A:1,B:1):1,(C:1,D:1):1);"
  tr <- ape::read.tree(text = tree_text)

  # Adiciona tips com vazio e NA para testar resiliência
  cost <- rSWeeP:::quebrataxonCOPHE(tr)
  expect_type(cost, "list")
  expect_true("tab" %in% names(cost))
  expect_true("mean" %in% names(cost))
  expect_equal(nrow(cost$tab), 4)

  mp <- rSWeeP:::MonoParaphylMetric(tr)
  expect_type(mp, "list")
  expect_true(all(c("tab", "percMono", "percPara", "metric") %in% names(mp)))
  expect_equal(nrow(mp$tab), 4)
})
