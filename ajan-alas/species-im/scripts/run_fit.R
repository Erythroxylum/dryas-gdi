# Run the nine chromosome-specific empirical BPP A00 fits for the collapsed
# ajanensis vs alaskensis two-lineage IM model.
#
# Required environment variable:
#   BPP_BIN=/full/path/to/bpp
#
# Usage from the repository root:
#   BPP_BIN=/full/path/to/bpp Rscript ajan-alas/species-im/scripts/run_fit.R
#
# Fits are run sequentially so failures are obvious and log files remain small.

suppressPackageStartupMessages({
  library(readr)
})

# Locate this species-IM analysis directory from the script path.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Read the BPP executable path and make sure it works before starting long fits.
bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) {
  stop("BPP could not be run. Set BPP_BIN to the full path to bpp.")
}

# Confirm that all nine fit controls have already been generated.
controls <- file.path("fit", "controls", paste0("ch", 1:9, ".ctl"))
missing <- controls[!file.exists(controls)]
if (length(missing)) {
  stop(
    "Missing fit controls: ", paste(missing, collapse = ", "),
    ". Run scripts/build_fit_controls.R first."
  )
}

# Keep a plain-text console log for each chromosome in addition to BPP's own
# mcmcfile/outfile outputs specified inside the control file.
dir.create("fit/logs", recursive = TRUE, showWarnings = FALSE)

# Run each chromosome sequentially. Existing output files are not deleted here;
# this protects completed long runs from accidental removal. BPP itself will report
# if an output path needs attention.
for (chrom in paste0("ch", 1:9)) {
  ctl <- file.path("fit", "controls", paste0(chrom, ".ctl"))
  log_file <- file.path("fit", "logs", paste0(chrom, ".log.txt"))

  message("Starting species-level IM fit for ", chrom)
  result <- system2(
    bpp,
    args = ctl,
    stdout = log_file,
    stderr = log_file
  )

  if (!identical(result, 0L)) {
    stop("BPP failed for ", chrom, ". Inspect ", log_file)
  }

  message("Finished ", chrom)
}

message("All nine species-level ajan vs alas IM fits completed")
