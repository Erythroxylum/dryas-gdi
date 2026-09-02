# Generate BPP simulation control files for the two PoOcPu1 gdi comparisons.
#
# This script:
#   1. Locates the poocpu analysis directory.
#   2. Loads the chromosome-specific posterior means and model definition.
#   3. Validates that all required model objects and parameter columns exist.
#   4. Creates three-sequence aab/abb imap files.
#   5. Writes 36 BPP control files: 9 chromosomes x 2 comparisons x 2 directions.

# Load packages used for reading tables and handling tibbles.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# Determine this script's location so it works regardless of the shell working directory.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Identify the two required inputs: extracted posterior means and the exact reduced model.
param_file <- "parameters/h4d-s47-p9-means.csv"
model_file <- "model/poocpu_model.R"
if (!file.exists(param_file)) {
  stop("Missing ", param_file, ". Run scripts/extract_parameters.R first.")
}
if (!file.exists(model_file)) {
  stop("Missing ", model_file, ".")
}

# Source the model into the global environment. This makes its objects directly
# available to the remainder of this top-level Rscript process.
source(model_file, local = globalenv())

# Confirm that the model file defined every object needed by this script.
required_objects <- c(
  "population_order",
  "comparison_specs",
  "make_species_tree",
  "migration_pairs",
  "migration_parameter_names",
  "required_parameter_columns"
)
missing_objects <- required_objects[
  !vapply(
    required_objects,
    function(x) exists(x, envir = globalenv(), inherits = FALSE),
    logical(1)
  )
]
if (length(missing_objects)) {
  stop("Model file is missing: ", paste(missing_objects, collapse = ", "))
}

# Read the chromosome-specific posterior means and verify that all model parameters
# required for the PoOcPu1 controls are present.
p <- read_csv(param_file, show_col_types = FALSE)
missing_cols <- setdiff(required_parameter_columns, names(p))
if (length(missing_cols)) {
  stop("Parameter table is missing columns: ", paste(missing_cols, collapse = ", "))
}
if (!identical(as.character(p$chromosome), paste0("ch", 1:9))) {
  stop("Expected chromosome rows ch1 through ch9 in order")
}

# Create output directories. BPP will not create parent directories for tree files.
dir.create("controls", showWarnings = FALSE)
dir.create(file.path("output", "trees"), recursive = TRUE, showWarnings = FALSE)
dir.create("imap", showWarnings = FALSE)

# Write a three-sequence sampling map for one comparison/configuration.
# aab = two samples from population A and one from B.
# abb = one sample from A and two from B.
write_imap <- function(comparison, config, pop_a, pop_b) {
  f <- file.path("imap", paste0(comparison, "_", config, ".imap.txt"))

  if (config == "aab") {
    x <- c(
      paste0("a1 ", pop_a),
      paste0("a2 ", pop_a),
      paste0("b1 ", pop_b)
    )
  } else {
    x <- c(
      paste0("a1 ", pop_a),
      paste0("b1 ", pop_b),
      paste0("b2 ", pop_b)
    )
  }

  writeLines(x, f)
  f
}

# Build one BPP simulation control file from one chromosome's posterior means
# and one row of the comparison specification table.
make_ctl <- function(row, spec) {
  # Set all five PoOcPu1 populations to zero samples, then place the three focal
  # sequences into the two populations being tested.
  counts <- setNames(rep(0L, length(population_order)), population_order)
  if (spec$config == "aab") {
    counts[spec$pop_a] <- 2L
    counts[spec$pop_b] <- 1L
  } else {
    counts[spec$pop_a] <- 1L
    counts[spec$pop_b] <- 2L
  }

  # Create the matching imap file and unique output tree path.
  imap <- write_imap(spec$comparison, spec$config, spec$pop_a, spec$pop_b)
  treefile <- file.path(
    "output",
    "trees",
    paste0(row$chromosome, "_", spec$comparison, "_", spec$config, ".tree.txt")
  )

  # Convert the chromosome-specific tau/theta means into the numeric BPP tree.
  tree <- make_species_tree(row)

  # Write every migration edge retained in the reduced PoOcPu1 model, mapping
  # reduced-model edge names back to the corresponding H4d posterior parameter.
  mig <- vapply(migration_pairs, function(pair) {
    edge_name <- paste0("W_", pair[1], "_to_", pair[2])
    param_name <- unname(migration_parameter_names[edge_name])

    if (is.na(param_name) || !nzchar(param_name)) {
      stop("No parameter mapping for migration edge: ", edge_name)
    }
    if (!param_name %in% names(row)) {
      stop("Missing parameter column: ", param_name)
    }

    sprintf("            %s %s %.10g", pair[1], pair[2], row[[param_name]])
  }, character(1))

  # Assemble the final BPP control file. One million simulated gene trees are
  # generated per chromosome/configuration, matching the integrifolia workflow.
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

# Loop over nine chromosomes and four directional tests, producing 36 controls.
for (i in seq_len(nrow(p))) {
  row <- p[i, ]

  for (j in seq_len(nrow(comparison_specs))) {
    spec <- comparison_specs[j, ]
    out <- file.path(
      "controls",
      paste0(row$chromosome, "_", spec$comparison, "_", spec$config, ".ctl")
    )
    writeLines(make_ctl(row, spec), out)
  }
}

message("Generated 36 PoOcPu1 gdi controls in poocpu/controls/")
