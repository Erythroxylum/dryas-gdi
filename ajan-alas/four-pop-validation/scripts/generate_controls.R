# Generate BPP simulation controls for the four-population validation analysis.
#
# This script:
#   1. Loads chromosome-specific posterior means from the fitted four-pop model.
#   2. Loads the exact four-population topology and migration edges.
#   3. Generates aab/abb three-sequence sampling schemes for the two within-species tests.
#   4. Writes 36 controls: 9 chromosomes x 2 comparisons x 2 directions.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# Find this analysis directory from the script location.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Define required inputs.
param_file <- "parameters/four_pop_means.csv"
model_file <- "model/four_pop_model.R"
if (!file.exists(param_file)) stop("Missing ", param_file, ". Run scripts/extract_parameters.R first.")
if (!file.exists(model_file)) stop("Missing ", model_file)

# Load the exact model objects into this Rscript process.
source(model_file, local = globalenv())

# Validate that the model definition is complete.
required_objects <- c(
  "population_order",
  "comparison_specs",
  "migration_pairs",
  "migration_parameter_names",
  "make_species_tree",
  "required_parameter_columns"
)
missing_objects <- required_objects[
  !vapply(required_objects, function(x) exists(x, envir = globalenv(), inherits = FALSE), logical(1))
]
if (length(missing_objects)) stop("Model file is missing: ", paste(missing_objects, collapse = ", "))

# Read and validate chromosome-specific posterior means.
p <- read_csv(param_file, show_col_types = FALSE)
missing_cols <- setdiff(required_parameter_columns, names(p))
if (length(missing_cols)) stop("Parameter table is missing: ", paste(missing_cols, collapse = ", "))
if (!identical(as.character(p$chromosome), paste0("ch", 1:9))) stop("Expected chromosome rows ch1 through ch9 in order")

# Create output directories before BPP is called.
dir.create("controls", showWarnings = FALSE)
dir.create("imap", showWarnings = FALSE)
dir.create(file.path("output", "trees"), recursive = TRUE, showWarnings = FALSE)

# Write a three-sequence sampling map for one directional gdi calculation.
write_imap <- function(spec) {
  f <- file.path("imap", paste0(spec$comparison, "_", spec$config, ".imap.txt"))
  if (spec$config == "aab") {
    x <- c(
      paste0("a1 ", spec$pop_a),
      paste0("a2 ", spec$pop_a),
      paste0("b1 ", spec$pop_b)
    )
  } else {
    x <- c(
      paste0("a1 ", spec$pop_a),
      paste0("b1 ", spec$pop_b),
      paste0("b2 ", spec$pop_b)
    )
  }
  writeLines(x, f)
  f
}

# Build one numeric BPP simulation control from one chromosome and one direction.
make_ctl <- function(row, spec) {
  # Assign the three simulated samples to the focal pair while leaving the other
  # two fitted populations unsampled but present in the demographic model.
  counts <- setNames(rep(0L, length(population_order)), population_order)
  if (spec$config == "aab") {
    counts[spec$pop_a] <- 2L
    counts[spec$pop_b] <- 1L
  } else {
    counts[spec$pop_a] <- 1L
    counts[spec$pop_b] <- 2L
  }

  # Create a unique imap and tree output path.
  imap <- write_imap(spec)
  treefile <- file.path(
    "output", "trees",
    paste0(row$chromosome, "_", spec$comparison, "_", spec$config, ".tree.txt")
  )

  # Insert chromosome-specific tau/theta posterior means into the population tree.
  tree <- make_species_tree(row)

  # Retain all four migration rates from the fitted model.
  mig <- vapply(migration_pairs, function(pair) {
    edge_name <- paste0("W_", pair[1], "_to_", pair[2])
    param_name <- unname(migration_parameter_names[edge_name])
    if (is.na(param_name) || !param_name %in% names(row)) stop("Missing migration parameter for ", edge_name)
    sprintf("            %s %s %.10g", pair[1], pair[2], row[[param_name]])
  }, character(1))

  # Assemble the BPP simulation control.
  paste(c(
    "seed = -1",
    paste("treefile =", treefile),
    paste("Imapfile =", imap),
    paste("species&tree =", length(population_order), paste(population_order, collapse = " ")),
    paste("                ", paste(unname(counts), collapse = " ")),
    tree,
    "loci&length = 1000000 50",
    paste("migration =", length(migration_pairs)),
    mig
  ), collapse = "\n")
}

# Generate all chromosome/configuration controls.
for (i in seq_len(nrow(p))) {
  row <- p[i, ]
  for (j in seq_len(nrow(comparison_specs))) {
    spec <- comparison_specs[j, ]
    outfile <- file.path(
      "controls",
      paste0(row$chromosome, "_", spec$comparison, "_", spec$config, ".ctl")
    )
    writeLines(make_ctl(row, spec), outfile)
  }
}

message("Generated 36 four-population validation controls")
