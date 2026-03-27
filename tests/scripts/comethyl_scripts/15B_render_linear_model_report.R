
---

## 15B: render script

```r
#!/usr/bin/env Rscript
message("Starting ✓")

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

required <- c("--rmd", "--me_file", "--trait_file", "--formulas_file", "--module_name", "--output_dir")
missing_required <- required[sapply(required, function(x) is.null(get_arg(x)))]

if (length(missing_required) > 0) {
  stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
}

rmd_file    <- get_arg("--rmd")
me_file     <- get_arg("--me_file")
trait_file  <- get_arg("--trait_file")
formulas_file <- get_arg("--formulas_file")
module_name <- get_arg("--module_name")
output_dir  <- get_arg("--output_dir")
sample_id_col <- get_arg("--sample_id_col", "Sample_ID")
trait_exclude_file <- get_arg("--trait_exclude_file", "")
tech_vars_file <- get_arg("--tech_vars_file", "")
output_file <- get_arg("--output_file", paste0(module_name, "_linear_model_report.pdf"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  input = rmd_file,
  output_file = output_file,
  output_dir = output_dir,
  params = list(
    me_file = me_file,
    trait_file = trait_file,
    formulas_file = formulas_file,
    module_name = module_name,
    output_dir = output_dir,
    sample_id_col = sample_id_col,
    trait_exclude_file = trait_exclude_file,
    tech_vars_file = tech_vars_file
  ),
  envir = new.env(parent = globalenv())
)