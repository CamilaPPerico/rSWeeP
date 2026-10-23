# helper-golden.R: Utilitários permanentes de leitura e validação dos Golden Files

get_fixtures_dir <- function() {
  testthat::test_path("fixtures", "golden")
}

get_manifest_path <- function() {
  file.path(get_fixtures_dir(), "manifest.json")
}

load_golden_manifest <- function() {
  mpath <- get_manifest_path()
  if (!file.exists(mpath)) {
    stop(paste("Manifest not found at:", mpath))
  }
  jsonlite::fromJSON(mpath, simplifyVector = FALSE)
}

verify_golden_checksum <- function(filename) {
  filepath <- file.path(get_fixtures_dir(), filename)
  if (!file.exists(filepath)) {
    stop(paste("Golden file not found:", filepath))
  }
  manifest <- load_golden_manifest()
  entry <- manifest$files[[filename]]
  if (is.null(entry)) {
    stop(paste("Filename not registered in manifest:", filename))
  }
  actual_sha <- digest::digest(file = filepath, algo = "sha256")
  if (!identical(actual_sha, entry$sha256)) {
    stop(sprintf("Checksum mismatch for %s: expected %s, got %s",
                 filename, entry$sha256, actual_sha))
  }
  TRUE
}

load_golden_csv <- function(filename) {
  verify_golden_checksum(filename)
  filepath <- file.path(get_fixtures_dir(), filename)
  # Leitura em float64 preservando cabeçalhos e dimensões
  df <- utils::read.csv(filepath, header = TRUE, row.names = 1, check.names = FALSE)
  as.matrix(df)
}

load_golden_mtx <- function(filename) {
  verify_golden_checksum(filename)
  filepath <- file.path(get_fixtures_dir(), filename)
  Matrix::readMM(filepath)
}
