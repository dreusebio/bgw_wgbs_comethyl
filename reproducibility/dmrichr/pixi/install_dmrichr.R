options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes",     quietly = TRUE)) install.packages("remotes")

# ── Bioconductor packages not available as conda binaries ──────────────────
bioc_to_install <- c(
  "GenomeInfoDbData",
  "GO.db",
  "HDO.db",          # <-- ADD: required by ChIPseeker
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ChIPseeker"
)

missing_bioc <- bioc_to_install[
  !vapply(bioc_to_install, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# ── CRAN packages not in conda ──────────────────────────────────────────────
cran_to_install <- c("hablar")
missing_cran <- cran_to_install[
  !vapply(cran_to_install, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# ── Final check before installing DMRichR ───────────────────────────────────
required_pkgs <- c(
  "GenomeInfoDbData", "GO.db","HDO.db", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ChIPseeker", "lsmeans", "R2HTML",
  "hablar", "wesanderson", "PCAtools"
)

still_missing <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(still_missing) > 0) {
  stop("Still missing after install attempts: ", paste(still_missing, collapse = ", "))
}

# ── Install DMRichR ──────────────────────────────────────────────────────────
remotes::install_github(
  "ben-laufer/DMRichR",
  upgrade      = "never",
  dependencies = FALSE
)

if (requireNamespace("DMRichR", quietly = TRUE)) {
  message("✓ DMRichR installation complete")
} else {
  stop("DMRichR installation failed")
}