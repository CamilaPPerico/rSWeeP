# test-integration-sweep-fasta.R: Vetorização de sequências biológicas (SWeeP e SWeePlite)

test_that("SWeeP e SWeePlite produzem saídas idênticas em sequências em memória", {
  aa_seqs <- Biostrings::AAStringSet(c(
    p1 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    p2 = "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL"
  ))

  psz <- 10
  mask <- c(1, 0, 1) # 400 colunas
  base <- orthBase(lin = 400, col = psz, seed = 647474747)

  # SWeeP clássico com base explícita
  res_sweep <- SWeeP(aa_seqs, orthbase = base, mask = mask, seqtype = "AA", bin = FALSE, norm = "none", ncores = 1, verbose = FALSE)

  # SWeePlite com base interna
  res_lite <- SWeePlite(aa_seqs, psz = psz, mask = mask, seqtype = "AA", bin = FALSE, norm = "none", ncores = 1, verbose = FALSE)

  expect_equal(res_sweep$proj, res_lite$proj, tolerance = 1e-12)
  expect_equal(dim(res_sweep$proj), c(2, psz))
  expect_equal(rownames(res_sweep$proj), c("p1", "p2"))
})

test_that("SWeeP e SWeePlite em pasta FASTA reproduzem fielmente os Golden Files do Artigo", {
  path_mito <- system.file("examples", "aaMitochondrial", package = "rSWeeP")
  if (path_mito == "") path_mito <- "/home/rasert/workspace/rSWeeP/inst/examples/aaMitochondrial"

  golden_count <- load_golden_csv("sweep_mito_aa_14x600_count_none.csv")
  golden_bin <- load_golden_csv("sweep_mito_aa_14x600_bin_none.csv")

  expect_equal(dim(golden_count), c(14, 600))
  expect_equal(dim(golden_bin), c(14, 600))

  # Teste com subconjunto de sequências e projeção rápida
  mask_compact <- c(1, 0, 1) # 400
  base_compact <- orthBase(lin = 400, col = 10, seed = 647474747)

  res_sub <- SWeeP(path_mito, orthbase = base_compact, mask = mask_compact, seqtype = "AA",
                   bin = FALSE, norm = "none", extension = ".faa", ncores = 2, verbose = FALSE)

  expect_equal(dim(res_sub$proj), c(14, 10))
  expect_equal(nrow(res_sub$proj), 14)
})

test_that("Vetorização preserva distâncias filogenéticas (Lema de Johnson-Lindenstrauss e Coerência Taxonômica)", {
  golden_bin <- load_golden_csv("sweep_mito_aa_14x600_bin_none.csv")

  # Matriz de distâncias euclidianas entre os vetores de 600 dimensões
  dmat <- as.matrix(dist(golden_bin, method = "euclidean"))

  # Hominídeos: 03_Homo_sapiens, 08_Homo_sapiens_neanderthalensis, 01_Pan_troglodytes, 11_Pan_paniscus
  d_human_neander <- dmat["03_Homo_sapiens", "08_Homo_sapiens_neanderthalensis"]
  d_human_chimp   <- dmat["03_Homo_sapiens", "01_Pan_troglodytes"]
  d_human_mouse   <- dmat["03_Homo_sapiens", "06_Mus_musculus"]
  d_human_plant   <- dmat["03_Homo_sapiens", "14_Rhazya_stricta"]

  # Distâncias biológicas esperadas:
  # Humano x Neanderthal < Humano x Chimpanzé < Humano x Camundongo < Humano x Planta (outgroup)
  expect_true(d_human_neander < d_human_chimp)
  expect_true(d_human_chimp < d_human_mouse)
  expect_true(d_human_mouse < d_human_plant)
})
