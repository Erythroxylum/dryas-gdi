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

# Locate this species-IM analysis directory from the script path so all relative
# control, log, and output paths resolve consistently.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Read the BPP executable path and make sure the binary is callable before
# starting any chromosome analyses.
bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) {
  stop("BPP could not be run. Set BPP_BIN to the full path to bpp.")
}

# Confirm that all nine empirical fit controls have already been generated.
controls <- file.path("fit", "controls", paste0("ch", 1:9, ".ctl"))
missing <- controls[!file.exists(controls)]
if (length(missing)) {
  stop(
    "Missing fit controls: ", paste(missing, collapse = ", "),
    ". Run scripts/build_fit_controls.R first."
  )
}

# Prepare the console-log directory. Each chromosome also writes BPP job output
# under fit/output/ according to the jobname directive in its control file.
dir.create("fit/logs", recursive = TRUE, showWarnings = FALSE)
dir.create("fit/output", recursive = TRUE, showWarnings = FALSE)

# Run each chromosome sequentially. BPP 4.8.x requires --cfile before the control
# filename for an empirical analysis. Passing the filename alone only prints BPP's
# usage screen and exits successfully, which can otherwise look like a completed fit.
for (chrom in paste0("ch", 1:9)) {
  ctl <- file.path("fit", "controls", paste0(chrom, ".ctl"))
  log_file <- file.path("fit", "logs", paste0(chrom, ".log.txt"))

  # Record existing output files so we can verify that this invocation actually
  # produced new BPP analysis output rather than merely returning exit status 0.
  before <- list.files("fit/output", full.names = TRUE)

  message("Starting species-level IM fit for ", chrom)
  result <- system2(
    bpp,
    args = c("--cfile", ctl),
    stdout = log_file,
    stderr = log_file
  )

  if (!identical(result, 0L)) {
    stop("BPP failed for ", chrom, ". Inspect ", log_file)
  }

  # Require at least one chromosome-specific output file after BPP returns. This
  # catches the exact failure mode where BPP prints its help text but does no fit.
  after <- list.files("fit/output", full.names = TRUE)
  chrom_outputs <- after[grepl(paste0("(^|/)", chrom, "([._]|$)"), after)]
  if (!length(chrom_outputs)) {
    stop(
      "BPP returned status 0 for ", chrom,
      " but produced no chromosome-specific output. Inspect ", log_file
    )
  }

  message("Finished ", chrom)
}

message("All nine species-level ajan vs alas IM fits completed")
