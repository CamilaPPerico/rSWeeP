# test-unit-checks.R: Validação de contratos e checagens de parâmetros

test_that("SW.checks valida seqtype corretamente", {
  # Valores válidos não devem lançar erro
  expect_silent(rSWeeP:::SW.checks("seqtype", "AA"))
  expect_silent(rSWeeP:::SW.checks("seqtype", "NT"))

  # Valores inválidos
  expect_error(rSWeeP:::SW.checks("seqtype", "DNA"), "'seqtype' should be any one of the options")
  expect_error(rSWeeP:::SW.checks("seqtype", c("AA", "NT")), "'seqtype' should be any one of the options")
  expect_error(rSWeeP:::SW.checks("seqtype", 123), "'seqtype' should be any one of the options")
})

test_that("SW.checks valida extension corretamente", {
  expect_silent(rSWeeP:::SW.checks("extension", ".faa"))
  expect_silent(rSWeeP:::SW.checks("extension", c(".faa", ".fasta")))
  expect_error(rSWeeP:::SW.checks("extension", 123), "'extension' should be a character")
})

test_that("SW.checks valida psz corretamente", {
  expect_silent(rSWeeP:::SW.checks("psz", 600))
  expect_silent(rSWeeP:::SW.checks("psz", 1369))
  expect_error(rSWeeP:::SW.checks("psz", 600.5), "'psz' should be a integer numeric value")
  expect_error(rSWeeP:::SW.checks("psz", "600"), "'psz' should be a integer numeric value")
})

test_that("SW.checks valida mask corretamente", {
  expect_silent(rSWeeP:::SW.checks("mask", c(2, 1, 2)))
  expect_silent(rSWeeP:::SW.checks("mask", c(5, 5, 5)))
  expect_error(rSWeeP:::SW.checks("mask", "11011"), "'mask' should be a vector")
})

test_that("SW.checks valida flags booleanas", {
  for (flag_name in c("bin", "readEachFASTA", "lowRAMmode", "RNAseqdata", "transpose", "concatenate")) {
    expect_silent(rSWeeP:::SW.checks(flag_name, TRUE))
    expect_silent(rSWeeP:::SW.checks(flag_name, FALSE))
    expect_error(rSWeeP:::SW.checks(flag_name, "TRUE"))
    expect_error(rSWeeP:::SW.checks(flag_name, 1))
  }
})

test_that("SW.checks valida ncores e fastalist", {
  expect_silent(rSWeeP:::SW.checks("ncores", NULL))
  expect_silent(rSWeeP:::SW.checks("ncores", 2))
  expect_error(rSWeeP:::SW.checks("ncores", 0), "'ncores' should be a positive integer")
  expect_error(rSWeeP:::SW.checks("ncores", -1), "'ncores' should be a positive integer")

  expect_silent(rSWeeP:::SW.checks("fastalist", c("a.faa", "b.faa")))
  expect_error(rSWeeP:::SW.checks("fastalist", character(0)), "The folder provided does not seem to contain files")
})

test_that("SW.checks valida norm para os 4 modos", {
  # [DOC-DIVERGENCE]: código aceita 4 modos ('none', 'log', 'Neg', 'logNeg')
  for (m in c("none", "log", "Neg", "logNeg")) {
    expect_silent(rSWeeP:::SW.checks("norm", m))
  }
  expect_error(rSWeeP:::SW.checks("norm", "invalid"), "'norm' must be one of the following options")
})

test_that("SW.checkORTHBASE valida matriz de projeção", {
  orth <- matrix(1:20, nrow = 5, ncol = 4)
  input <- matrix(1:15, nrow = 3, ncol = 5)

  # Válido: nrow(orth) == ncol(input)
  expect_silent(rSWeeP:::SW.checkORTHBASE(orth, input))

  # Dimensões incompatíveis
  input_err <- matrix(1:12, nrow = 3, ncol = 4)
  expect_error(rSWeeP:::SW.checkORTHBASE(orth, input_err), "The number of rows of the orthbase matrix must be equal")

  # Não-matriz
  expect_error(rSWeeP:::SW.checkORTHBASE(list(1, 2, 3), input), "Orthbase must be an array")
})
