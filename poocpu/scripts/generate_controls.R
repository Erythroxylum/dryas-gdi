suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

param_file <- "parameters/h4d-s47-p9-means.csv"
model_file <- "model/poocpu_model.R"
if (!file.exists(param_file)) stop("Missing ", param_file, ". Add chromosome-specific H4d posterior means first.")
if (!file.exists(model_file)) stop("Missing ", model_file, ". Add the exact PoOcPu1 demographic model specification first.")
source(model_file)

required_objects <- c("population_order", "comparison_specs", "make_species_tree", "migration_pairs")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = FALSE)]
if (length(missing_objects)) stop("Model file is missing: ", paste(missing_objects, collapse = ", "))

p <- read_csv(param_file, show_col_types = FALSE)
dir.create("controls", showWarnings = FALSE)
dir.create(file.path("output", "trees"), recursive = TRUE, showWarnings = FALSE)
dir.create("imap", showWarnings = FALSE)

write_imap <- function(comparison, config, pop_a, pop_b) {
  f <- file.path("imap", paste0(comparison, "_", config, ".imap.txt"))
  if (config == "aab") {
    x <- c(paste0("a1 ", pop_a), paste0("a2 ", pop_a), paste0("b1 ", pop_b))
  } else {
    x <- c(paste0("a1 ", pop_a), paste0("b1 ", pop_b), paste0("b2 ", pop_b))
  }
  writeLines(x, f)
  f
}

make_ctl <- function(row, spec) {
  counts <- setNames(rep(0L, length(population_order)), population_order)
  if (spec$config == "aab") {
    counts[spec$pop_a] <- 2L; counts[spec$pop_b] <- 1L
  } else {
    counts[spec$pop_a] <- 1L; counts[spec$pop_b] <- 2L
  }
  imap <- write_imap(spec$comparison, spec$config, spec$pop_a, spec$pop_b)
  treefile <- file.path("output", "trees", paste0(row$chromosome, "_", spec$comparison, "_", spec$config, ".tree.txt"))
  tree <- make_species_tree(row)
  mig <- vapply(migration_pairs, function(pair) {
    nm <- paste0("W_", pair[1], "_to_", pair[2])
    if (!nm %in% names(row)) stop("Missing parameter column: ", nm)
    sprintf("            %s %s %.10g", pair[1], pair[2], row[[nm]])
  }, character(1))
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

for (i in seq_len(nrow(p))) {
  row <- p[i, ]
  for (j in seq_len(nrow(comparison_specs))) {
    spec <- comparison_specs[j, ]
    out <- file.path("controls", paste0(row$chromosome, "_", spec$comparison, "_", spec$config, ".ctl"))
    writeLines(make_ctl(row, spec), out)
  }
}
message("Generated PoOcPu1 gdi controls in poocpu/controls/")
