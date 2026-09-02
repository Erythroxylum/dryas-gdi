# Run the nine chromosome-specific empirical BPP A00 fits for the collapsed
# ajanensis vs alaskensis two-lineage IM model.
#
# Required environment variable:
#   BPP_BIN=/full/path/to/bpp
#
# Optional environment variable:
#   BPP_JOBS=<number of chromosome fits to run simultaneously>
#
# Usage from the repository root:
#   BPP_BIN=/full/path/to/bpp BPP_JOBS=4 \
#     Rscript ajan-alas/species-im/scripts/run_fit.R
#
# The nine chromosome fits are statistically independent, so this script runs
# several BPP processes at the same time. Each process has its own control file,
# log, and job output. BPP_JOBS is deliberately user-configurable because the
# best value depends on the CPU and on how many resources one BPP process uses.

suppressPackageStartupMessages({
  library(parallel)
})

# Locate this species-IM analysis directory from the script path so all relative
# control, log, and output paths resolve consistently on any computer.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Read the BPP executable path and confirm that the binary is callable before
# starting any long-running chromosome analyses.
bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) {
  stop("BPP could not be run. Set BPP_BIN to the full path to bpp.")
}

# Read the requested parallelism. Default to two simultaneous chromosome fits,
# which is conservative for a workstation. On a machine with ample CPU/RAM,
# set BPP_JOBS as high as 9 to launch all chromosomes at once.
requested_jobs <- Sys.getenv("BPP_JOBS", unset = "2")
n_jobs <- suppressWarnings(as.integer(requested_jobs))
if (is.na(n_jobs) || n_jobs < 1L) {
  stop("BPP_JOBS must be a positive integer")
}
n_jobs <- min(n_jobs, 9L)

# Report available hardware as a guide, but do not automatically consume every
# detected core because each BPP process may itself use substantial CPU and RAM.
detected_cores <- parallel::detectCores(logical = TRUE)
message(
  "Parallel configuration: BPP_JOBS=", n_jobs,
  if (!is.na(detected_cores)) paste0("; detected logical cores=", detected_cores) else ""
)

# Confirm that all nine empirical fit controls have already been generated.
chromosomes <- paste0("ch", 1:9)
controls <- file.path("fit", "controls", paste0(chromosomes, ".ctl"))
missing <- controls[!file.exists(controls)]
if (length(missing)) {
  stop(
    "Missing fit controls: ", paste(missing, collapse = ", "),
    ". Run scripts/build_fit_controls.R first."
  )
}

# Prepare separate directories for BPP console logs and BPP job outputs. Each
# chromosome writes to unique paths, which prevents collisions between workers.
dir.create("fit/logs", recursive = TRUE, showWarnings = FALSE)
dir.create("fit/output", recursive = TRUE, showWarnings = FALSE)

# Refuse to start if chromosome-specific BPP outputs already exist. This protects
# completed or partially completed runs from accidental overwrite. Console logs
# alone are safe and will be replaced when a chromosome is rerun.
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

# Run one chromosome fit. This function is executed inside a parallel worker.
# BPP 4.8.x requires --cfile before an empirical-analysis control filename.
run_one_chromosome <- function(chrom, bpp_path, analysis_path) {
  setwd(analysis_path)

  ctl <- file.path("fit", "controls", paste0(chrom, ".ctl"))
  log_file <- file.path("fit", "logs", paste0(chrom, ".log.txt"))

  start_time <- Sys.time()
  result <- system2(
    bpp_path,
    args = c("--cfile", ctl),
    stdout = log_file,
    stderr = log_file
  )
  end_time <- Sys.time()

  # Find files created under this chromosome's unique job prefix. Requiring real
  # output catches cases where BPP exits without actually conducting the fit.
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

# Create a portable PSOCK cluster. Unlike mclapply(), PSOCK workers work on both
# Unix-like systems and Windows, so the same script can be moved to another
# workstation or compute server without changing the parallel implementation.
cluster <- parallel::makeCluster(n_jobs)
on.exit(parallel::stopCluster(cluster), add = TRUE)

# Launch the nine independent chromosome fits across the worker pool. parLapply
# waits until all assigned fits finish, while each BPP process streams its console
# output directly to its own chromosome log file.
message(
  "Starting 9 species-level IM fits with up to ", n_jobs,
  " simultaneous BPP processes"
)
results <- parallel::parLapply(
  cluster,
  chromosomes,
  run_one_chromosome,
  bpp_path = bpp,
  analysis_path = analysis_dir
)

# Stop the worker processes as soon as all chromosome jobs have returned.
parallel::stopCluster(cluster)
on.exit(NULL, add = FALSE)

# Summarize each chromosome and identify any failed or output-free analyses.
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

# Reaching this message means all nine BPP processes returned successfully and
# produced chromosome-specific files under fit/output/.
message("All nine species-level ajan vs alas IM fits completed successfully")
