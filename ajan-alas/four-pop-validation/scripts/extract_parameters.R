# Extract chromosome-specific posterior means from the final four-population
# ajanensis/alaskensis m3-prior2 formatted posterior table.
#
# The expected formatted input has parameters as rows and chromosomes as columns:
#   variable, ch1, ch2, ..., ch9
# with cells such as:
#   0.00123 (0.00100, 0.00145)
#
# For compatibility, this script also accepts Yuttapong's alternate table layout
# with one chromosome per row and BPP parameter names as columns.
#
# Usage from anywhere:
#   Rscript scripts/extract_parameters.R /path/to/ajan-alas-m3-prior2-s20-p4-fmt.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# Read and resolve the user-supplied path BEFORE changing working directories.
# This avoids accidentally interpreting a repo-relative path from inside the
# analysis subdirectory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript scripts/extract_parameters.R /path/to/ajan-alas-m3-prior2-s20-p4-fmt.csv")
}
source_file <- normalizePath(args[1], mustWork = TRUE)

# Locate the analysis directory from this script's path, then work there so all
# generated outputs have stable repository-relative locations.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Read the source CSV exactly as supplied by Yuttapong's formatted output.
x <- read_csv(source_file, show_col_types = FALSE)

# Define the exact BPP parameters used by the final m3-prior2 four-population
# model and the shorter output names used by our simulation workflow.
parameter_map <- c(
  theta_ajan_Interior = "theta:1:ajan_Interior",
  theta_ajan_Seward = "theta:2:ajan_Seward",
  theta_alas_Interior = "theta:3:alas_Interior",
  theta_alas_Seward = "theta:4:alas_Seward",
  theta_R = "theta:5:R",
  theta_J = "theta:6:J",
  theta_L = "theta:7:L",
  tau_R = "tau:5:R",
  tau_J = "tau:6:J",
  tau_L = "tau:7:L",
  W_alas_Seward_to_ajan_Seward = "W:4->2:alas_Seward->ajan_Seward",
  W_ajan_Seward_to_alas_Seward = "W:2->4:ajan_Seward->alas_Seward",
  W_alas_Interior_to_ajan_Interior = "W:3->1:alas_Interior->ajan_Interior",
  W_ajan_Interior_to_alas_Interior = "W:1->3:ajan_Interior->alas_Interior"
)

# Helper: extract only the posterior mean from a formatted cell such as
# "0.00123 (0.00100, 0.00145)". Scientific notation is also supported.
extract_mean <- function(z) {
  z <- as.character(z)
  value <- str_extract(str_trim(z), "^[+-]?(?:[0-9]*\\.?[0-9]+)(?:[eE][+-]?[0-9]+)?")
  as.numeric(value)
}

# Detect the formatted orientation used by the uploaded m3-prior2 CSV:
# one parameter per row, with variable + ch1...ch9 columns.
if ("variable" %in% names(x) && all(paste0("ch", 1:9) %in% names(x))) {

  # Verify that every required model parameter is present exactly once.
  missing_parameters <- setdiff(unname(parameter_map), x$variable)
  if (length(missing_parameters)) {
    stop("Formatted table is missing parameters: ", paste(missing_parameters, collapse = ", "))
  }
  duplicated_parameters <- unname(parameter_map)[vapply(
    unname(parameter_map),
    function(v) sum(x$variable == v) != 1,
    logical(1)
  )]
  if (length(duplicated_parameters)) {
    stop("Expected exactly one row for each parameter: ", paste(duplicated_parameters, collapse = ", "))
  }

  # Keep only required rows, extract means from the nine chromosome columns,
  # transpose to one row per chromosome, and rename parameters for the workflow.
  selected <- x |>
    filter(variable %in% unname(parameter_map)) |>
    select(variable, all_of(paste0("ch", 1:9))) |>
    mutate(across(starts_with("ch"), extract_mean))

  out <- selected |>
    pivot_longer(cols = starts_with("ch"), names_to = "chromosome", values_to = "value") |>
    mutate(output_name = names(parameter_map)[match(variable, unname(parameter_map))]) |>
    select(chromosome, output_name, value) |>
    pivot_wider(names_from = output_name, values_from = value) |>
    arrange(match(chromosome, paste0("ch", 1:9)))

# Otherwise accept the alternate orientation: one chromosome per row with
# BPP parameter names as individual columns.
} else {
  if ("chr" %in% names(x)) x <- rename(x, chromosome = chr)
  if (!"chromosome" %in% names(x)) {
    stop(
      "Could not recognize table layout. Expected either variable + ch1...ch9 columns, ",
      "or a chr/chromosome column with BPP parameters as columns. Found columns: ",
      paste(names(x), collapse = ", ")
    )
  }

  missing_parameters <- setdiff(unname(parameter_map), names(x))
  if (length(missing_parameters)) {
    stop("Table is missing parameter columns: ", paste(missing_parameters, collapse = ", "))
  }

  # Select and rename the required parameters while extracting means if cells
  # still contain parenthesized credible intervals.
  out <- tibble(chromosome = as.character(x$chromosome))
  for (output_name in names(parameter_map)) {
    source_name <- unname(parameter_map[[output_name]])
    out[[output_name]] <- extract_mean(x[[source_name]])
  }
}

# Verify chromosome ordering and that every extracted posterior mean is numeric.
if (!identical(as.character(out$chromosome), paste0("ch", 1:9))) {
  stop("Expected chromosome rows ch1 through ch9 in order")
}
numeric_columns <- setdiff(names(out), "chromosome")
if (any(vapply(out[numeric_columns], function(z) any(is.na(z)), logical(1)))) {
  bad <- numeric_columns[vapply(out[numeric_columns], function(z) any(is.na(z)), logical(1))]
  stop("Could not parse numeric posterior means for: ", paste(bad, collapse = ", "))
}

# Write the standardized chromosome-by-parameter table used to generate BPP
# validation simulations.
dir.create("parameters", showWarnings = FALSE)
write_csv(out, "parameters/four_pop_means.csv")
message("Wrote four-population posterior means to parameters/four_pop_means.csv")
