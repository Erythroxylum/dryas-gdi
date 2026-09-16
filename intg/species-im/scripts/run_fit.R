# Run nine chromosome-specific empirical BPP A00 fits for the direct
# two-population D. integrifolia isolation-with-migration model.
#
# Required environment variable if BPP is not on PATH:
#   BPP_BIN=/full/path/to/bpp
#
# Optional environment variable:
#   BPP_JOBS=<number of chromosome fits to run simultaneously>
#
# Example from the repository root:
#   BPP_BIN=/full/path/to/bpp BPP_JOBS=9 \
#     Rscript intg/species-im/scripts/run_fit.R
#
# Each chromosome is an independent BPP fit. --no-pin leaves CPU placement to
# the operating system when several BPP processes are running simultaneously.

suppressPackageStartupMessages({
  library(parallel)
})

# Locate intg/species-im from this script's path.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Confirm that BPP is callable.
bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) {
  stop("BPP could not be run. Set BPP_BIN to the full path to bpp.")
}

# Read requested parallelism and cap it at the nine chromosome analyses.
requested_jobs <- Sys.getenv("BPP_JOBS", unset = "2")
n_jobs <- suppressWarnings(as.integer(requested_jobs))
if (is.na(n_jobs) || n_jobs < 1L) stop("BPP_JOBS must be a positive integer")
n_jobs <- min(n_jobs, 9L)

message(
  "Parallel configuration: BPP_JOBS=", n_jobs,
  "; detected logical cores=", parallel::detectCores(logical = TRUE)
)

chromosomes <- paste0("ch", 1:9)
controls <- file.path("fit", "controls", paste0(chromosomes, ".ctl"))
missing <- controls[!file.exists(controls)]
if (length(missing)) {
  stop(
    "Missing fit controls: ", paste(missing, collapse = ", "),
    ". Run scripts/build_fit_controls.R first."
  )
}

dir.create("fit/logs", recursive = TRUE, showWarnings = FALSE)
dir.create("fit/output", recursive = TRUE, showWarnings = FALSE)

# Protect completed or partial output from accidental overwrite.
existing_outputs <- unlist(lapply(chromosomes, function(chrom) {
  list.files(
    "fit/output",
    pattern = paste0("^", chrom, "([._]|$)"),
    full.names = TRUE
  )
}), use.names = FALSE)
if (length(existing_outputs)) {
  stop(
    "Existing chromosome outputs found in fit/output/. Move or remove them before restarting: ",
    paste(existing_outputs, collapse = ", ")
  )
}

# Run one chromosome inside a parallel worker and record its console output in a
# chromosome-specific log file.
run_one_chromosome <- function(chrom, bpp_path, analysis_path) {
  setwd(analysis_path)
  ctl <- file.path("fit", "controls", paste0(chrom, ".ctl"))
  log_file <- file.path("fit", "logs", paste0(chrom, ".log.txt"))

  start_time <- Sys.time()
  result <- system2(
    bpp_path,
    args = c("--no-pin", "--cfile", ctl),
    stdout = log_file,
    stderr = log_file
  )
  end_time <- Sys.time()

  chrom_outputs <- list.files(
    "fit/output",
    pattern = paste0("^", chrom, "([._]|$)"),
    full.names = TRUE
  )

  list(
    chromosome = chrom,
    status = as.integer(result),
    log_file = log_file,
    outputs = chrom_outputs,
    elapsed_hours = as.numeric(difftime(end_time, start_time, units = "hours"))
  )
}

# PSOCK workers make the parallel implementation portable across platforms.
cluster <- parallel::makeCluster(n_jobs)
on.exit(parallel::stopCluster(cluster), add = TRUE)

message(
  "Starting 9 direct integrifolia IM fits with up to ", n_jobs,
  " simultaneous BPP processes"
)
results <- parallel::parLapply(
  cluster,
  chromosomes,
  run_one_chromosome,
  bpp_path = bpp,
  analysis_path = analysis_dir
)

parallel::stopCluster(cluster)
on.exit(NULL, add = FALSE)

# Summarize return status and verify that every successful process produced BPP
# output files rather than merely returning without a fatal shell error.
for (x in results) {
  message(
    x$chromosome, ": status=", x$status,
    ", elapsed=", sprintf("%.2f", x$elapsed_hours), " h",
    ", outputs=", length(x$outputs)
  )
}

failed <- vapply(
  results,
  function(x) x$status != 0L || length(x$outputs) == 0L,
  logical(1)
)
if (any(failed)) {
  bad <- vapply(results[failed], `[[`, character(1), "chromosome")
  stop(
    "One or more BPP fits failed or produced no output: ", paste(bad, collapse = ", "),
    ". Inspect the corresponding files in fit/logs/."
  )
}

message("All nine direct integrifolia two-population IM fits completed successfully")
