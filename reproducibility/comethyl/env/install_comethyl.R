options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# CRAN packages required by comethyl but not auto-installed
cran_pkgs <- c(
  "ggdendro",
  "R.devices",
  "rlist"
)

missing_cran <- cran_pkgs[!vapply(
  cran_pkgs,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

bioc_pkgs <- c(
  "GenomeInfoDb",
  "GenomeInfoDbData",
  "GenomicRanges",
  "IRanges",
  "S4Vectors",
  "BiocGenerics",
  "DelayedArray",
  "DelayedMatrixStats",
  "MatrixGenerics",
  "SummarizedExperiment",
  "annotatr",
  "AnnotationDbi",
  "GO.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "bsseq",
  "dmrseq",
  "sva",
  "preprocessCore",
  "HDO.db",
  "impute"
)

missing_bioc <- bioc_pkgs[!vapply(
  bioc_pkgs,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

remotes::install_github(
  "cemordaunt/comethyl@9a21854",
  upgrade = "never",
  dependencies = FALSE
)