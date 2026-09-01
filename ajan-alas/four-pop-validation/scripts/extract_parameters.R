# Extract chromosome-specific posterior means from Yuttapong's existing
# four-population ajan/alaskensis gdi input table.
#
# Usage:
#   Rscript extract_parameters.R /path/to/est-all_chr-mean.csv
#
# The source table is the one created in Yuttapong's gdi.R workflow from
# ajan-alas-m3-prior2 (or the preferred prior3 equivalent). It contains one
# row per chromosome and BPP parameter names as columns.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# Locate the analysis directory from this script's path.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Read the source file supplied on the command line.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/extract_parameters.R /path/to/est-all_chr-mean.csv")
}
source_file <- normalizePath(args[1], mustWork = TRUE)
x <- read_csv(source_file, show_col_types = FALSE)

# Accept either Yuttapong's 'chr' column or our repository's 'chromosome' name.
if ("chr" %in% names(x)) x <- rename(x, chromosome = chr)
if (!"chromosome" %in% names(x)) stop("Source table needs a chr or chromosome column")

# These are the exact parameters used in Yuttapong's four-population model.
required <- c(
  "chromosome",
  "theta:1:ajan_Interior",
  "theta:2:ajan_Seward",
  "theta:3:alas_Interior",
  "theta:4:alas_Seward",
  "theta:5:R",
  "theta:6:J",
  "theta:7:L",
  "tau:5:R",
  "tau:6:J",
  "tau:7:L",
  "W:4->2:alas_Seward->ajan_Seward",
  "W:2->4:ajan_Seward->alas_Seward",
  "W:3->1:alas_Interior->ajan_Interior",
  "W:1->3:ajan_Interior->alas_Interior"
)
missing <- setdiff(required, names(x))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

# Rename parameters to compact repository-safe names while preserving the values.
out <- x |>
  transmute(
    chromosome,
    theta_ajan_Interior = `theta:1:ajan_Interior`,
    theta_ajan_Seward = `theta:2:ajan_Seward`,
    theta_alas_Interior = `theta:3:alas_Interior`,
    theta_alas_Seward = `theta:4:alas_Seward`,
    theta_R = `theta:5:R`,
    theta_J = `theta:6:J`,
    theta_L = `theta:7:L`,
    tau_R = `tau:5:R`,
    tau_J = `tau:6:J`,
    tau_L = `tau:7:L`,
    W_alas_Seward_to_ajan_Seward = `W:4->2:alas_Seward->ajan_Seward`,
    W_ajan_Seward_to_alas_Seward = `W:2->4:ajan_Seward->alas_Seward`,
    W_alas_Interior_to_ajan_Interior = `W:3->1:alas_Interior->ajan_Interior`,
    W_ajan_Interior_to_alas_Interior = `W:1->3:ajan_Interior->alas_Interior`
  )

# Verify chromosome ordering and write the clean parameter table used by simulation.
if (!identical(as.character(out$chromosome), paste0("ch", 1:9))) {
  stop("Expected chromosome rows ch1 through ch9 in order")
}
dir.create("parameters", showWarnings = FALSE)
write_csv(out, "parameters/four_pop_means.csv")
message("Wrote four-population posterior means to parameters/four_pop_means.csv")
