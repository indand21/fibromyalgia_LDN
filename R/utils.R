MW <- c(NTX = 341.4, BN = 343.4)  # g/mol; naltrexone, 6-beta-naltrexol

ngml_to_nM <- function(x, species) x * 1000 / MW[[species]]
nM_to_ngml <- function(x, species) x * MW[[species]] / 1000
mg_to_nmol <- function(mg, species) mg * 1e6 / MW[[species]]

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

v5_root <- function() getOption("v5.root", normalizePath(".", winslash = "/"))

source_all <- function(dir = file.path(v5_root(), "R")) {
  for (f in list.files(dir, pattern = "\\.R$", full.names = TRUE)) source(f)
  invisible(TRUE)
}
