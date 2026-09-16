# Extract chromosome-specific posterior means from the completed direct
# two-population D. integrifolia IM fits.
#
# Usage from anywhere in the repository:
#   Rscript intg/species-im/scripts/extract_fit_parameters.R
#
# Output:
#   intg/species-im/fit/species_im_means.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
})

# Locate intg/species-im from this script's path.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

chromosomes <- paste0("ch", 1:9)

# Read each BPP MCMC sample and calculate posterior means for the two tip thetas,
# root theta, root divergence time, and reciprocal migration rates. Column names
# are checked explicitly so a changed population order cannot silently produce
# mislabeled demographic parameters.
params <- map_dfr(chromosomes, function(chrom) {
  f <- file.path("fit", "output", paste0(chrom, ".mcmc.txt"))
  if (!file.exists(f)) stop("Missing completed MCMC file: ", f)

  x <- read_tsv(f, show_col_types = FALSE)

  expected <- c(
    "theta:1:intg_nGL_Nslope",
    "theta:2:intg_CAswGL",
    "theta:3:R",
    "tau:3:R",
    "W:1->2:intg_nGL_Nslope->intg_CAswGL",
    "W:2->1:intg_CAswGL->intg_nGL_Nslope"
  )
  missing <- setdiff(expected, names(x))
  if (length(missing)) {
    stop(
      "Unexpected BPP columns in ", f, ". Missing: ",
      paste(missing, collapse = ", ")
    )
  }

  tibble(
    chromosome = chrom,
    theta_intg_nGL_Nslope = mean(x$`theta:1:intg_nGL_Nslope`),
    theta_intg_CAswGL = mean(x$`theta:2:intg_CAswGL`),
    theta_R = mean(x$`theta:3:R`),
    tau_R = mean(x$`tau:3:R`),
    W_intg_nGL_Nslope_to_intg_CAswGL =
      mean(x$`W:1->2:intg_nGL_Nslope->intg_CAswGL`),
    W_intg_CAswGL_to_intg_nGL_Nslope =
      mean(x$`W:2->1:intg_CAswGL->intg_nGL_Nslope`)
  )
})

write_csv(params, file.path("fit", "species_im_means.csv"))
print(params)
message("Wrote intg/species-im/fit/species_im_means.csv")
