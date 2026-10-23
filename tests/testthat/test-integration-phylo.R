# test-integration-phylo.R: Métricas de avaliação filogenética (PCCI e PMPG)

test_that("PCCI calcula consistência cofenética taxonômica na árvore mitocondrial de exemplo", {
  tree_path <- system.file("examples", "tree_Mitochondrial.tree", package = "rSWeeP")
  if (tree_path == "") tree_path <- "/home/rasert/workspace/rSWeeP/inst/examples/tree_Mitochondrial.tree"
  meta_path <- system.file("examples", "metadata_mitochondrial.csv", package = "rSWeeP")
  if (meta_path == "") meta_path <- "/home/rasert/workspace/rSWeeP/inst/examples/metadata_mitochondrial.csv"

  tree <- ape::read.tree(tree_path)
  meta <- utils::read.csv(meta_path, header = TRUE)
  mt <- data.frame(sp = meta$fileName, family = meta$family)

  pcci_res <- PCCI(tree, mt)

  expect_type(pcci_res, "list")
  expect_true(all(c("tab", "mean") %in% names(pcci_res)))
  expect_equal(pcci_res$mean, 1.0, tolerance = 1e-6)

  golden_pcci <- load_golden_csv("pcci_mitochondrial.csv")
  expect_equal(unname(as.matrix(pcci_res$tab)), unname(golden_pcci))
})

test_that("PMPG calcula monofilia e parafilia na árvore mitocondrial de exemplo", {
  tree_path <- system.file("examples", "tree_Mitochondrial.tree", package = "rSWeeP")
  if (tree_path == "") tree_path <- "/home/rasert/workspace/rSWeeP/inst/examples/tree_Mitochondrial.tree"
  meta_path <- system.file("examples", "metadata_mitochondrial.csv", package = "rSWeeP")
  if (meta_path == "") meta_path <- "/home/rasert/workspace/rSWeeP/inst/examples/metadata_mitochondrial.csv"

  tree <- ape::read.tree(tree_path)
  meta <- utils::read.csv(meta_path, header = TRUE)
  mt <- data.frame(sp = meta$fileName, family = meta$family)

  pmpg_res <- PMPG(tree, mt)

  expect_type(pmpg_res, "list")
  # [DOC-DIVERGENCE]: código retorna 'metric' em vez de 'mean' conforme a docstring
  expect_true(all(c("tab", "percMono", "percPara", "metric") %in% names(pmpg_res)))

  expect_equal(pmpg_res$percMono, 100, tolerance = 1e-6)
  expect_equal(pmpg_res$percPara, 0, tolerance = 1e-6)
  expect_equal(pmpg_res$metric, 100, tolerance = 1e-6)

  golden_pmpg <- load_golden_csv("pmpg_mitochondrial.csv")
  expect_equal(unname(as.matrix(pmpg_res$tab)), unname(golden_pmpg))
})
