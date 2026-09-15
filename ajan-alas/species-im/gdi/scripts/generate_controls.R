# Generate chromosome-specific BPP gene-tree simulation controls for the
# species-level D. ajanensis vs D. alaskensis gdi analysis.
#
# This script uses posterior mean demographic parameters estimated by the
# empirical two-species IM fits in ../fit/species_im_means.csv. For every
# chromosome it generates two reciprocal three-sequence simulations:
#
#   aab = 2 ajan + 1 alas  -> gdi for D. ajanensis
#   abb = 1 ajan + 2 alas  -> gdi for D. alaskensis
#
# Each simulation retains the chromosome-specific tip/root theta values,
# divergence time tau_R, and bidirectional migration rates estimated from BPP.
#
# Usage from anywhere in the repository:
#   Rscript ajan-alas/species-im/gdi/scripts/generate_controls.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# Locate the gdi analysis directory from this script's own path so the workflow
# is portable across computers and does not depend on the current shell folder.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
gdi_dir <- dirname(dirname(script_path))
setwd(gdi_dir)

# Read the chromosome-specific posterior means from the completed empirical
# species-level isolation-with-migration fits.
parameter_file <- file.path("..", "fit", "species_im_means.csv")
if (!file.exists(parameter_file)) {
  stop(
    "Missing parameter file: ", parameter_file,
    ". Run the empirical-fit parameter extraction first."
  )
}
p <- read_csv(parameter_file, show_col_types = FALSE)

# Validate the exact columns needed to reconstruct each fitted demographic model.
required_columns <- c(
  "chromosome",
  "theta_ajan",
  "theta_alas",
  "theta_R",
  "tau_R",
  "W_ajan_to_alas",
  "W_alas_to_ajan"
)
missing_columns <- setdiff(required_columns, names(p))
if (length(missing_columns)) {
  stop("Parameter file is missing: ", paste(missing_columns, collapse = ", "))
}

# Require one parameter row for each of the nine nuclear chromosomes.
chromosomes <- paste0("ch", 1:9)
if (!setequal(p$chromosome, chromosomes) || nrow(p) != 9L) {
  stop("Expected exactly one parameter row for each chromosome ch1-ch9")
}
p <- p |>
  mutate(chromosome = factor(chromosome, levels = chromosomes)) |>
  arrange(chromosome) |>
  mutate(chromosome = as.character(chromosome))

# Create analysis directories before writing controls or simulation outputs.
dir.create("controls", recursive = TRUE, showWarnings = FALSE)
dir.create("imap", recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("output", "trees"), recursive = TRUE, showWarnings = FALSE)

# Write the two small imap files used to label simulated sequences. The names
# before the population labels become the individual labels in BPP gene trees.
aab_imap <- c(
  "ajan1\tajan",
  "ajan2\tajan",
  "alas1\talas"
)
abb_imap <- c(
  "ajan1\tajan",
  "alas1\talas",
  "alas2\talas"
)
writeLines(aab_imap, file.path("imap", "ajan_alas_aab.imap.txt"))
writeLines(abb_imap, file.path("imap", "ajan_alas_abb.imap.txt"))

# Define reciprocal sampling configurations. Only the sequence counts and imap
# differ; all fitted demographic parameters remain identical within a chromosome.
comparison_specs <- list(
  aab = list(
    counts = "2 1",
    imap = "imap/ajan_alas_aab.imap.txt"
  ),
  abb = list(
    counts = "1 2",
    imap = "imap/ajan_alas_abb.imap.txt"
  )
)

# Build one complete BPP simulation control from a chromosome parameter row.
make_ctl <- function(row, config) {
  spec <- comparison_specs[[config]]
  treefile <- sprintf("output/trees/%s_%s.tree.txt", row$chromosome, config)

  # Reconstruct the fitted two-species MSC-M tree. BPP expects tip theta after
  # each population name, root tau after R:, and root theta after #.
  tree <- sprintf(
    "(ajan #%.10g, alas #%.10g)R:%.10g #%.10g;",
    row$theta_ajan,
    row$theta_alas,
    row$tau_R,
    row$theta_R
  )

  # Retain both directional migration rates from the empirical chromosome fit.
  migration <- c(
    sprintf("            ajan alas %.10g", row$W_ajan_to_alas),
    sprintf("            alas ajan %.10g", row$W_alas_to_ajan)
  )

  paste(c(
    paste0("# Species-level ajan vs alas gdi simulation: ", row$chromosome, " ", config),
    "# Parameters are chromosome-specific posterior means from the empirical IM fit.",
    "# aab estimates gdi_ajan; abb estimates gdi_alas.",
    "seed = -1",
    paste("treefile =", treefile),
    paste("Imapfile =", spec$imap),
    "species&tree = 2 ajan alas",
    paste("                ", spec$counts),
    tree,
    "loci&length = 1000000 50",
    "migration = 2",
    migration
  ), collapse = "\n")
}

# Generate 18 controls: two reciprocal simulations for each of nine chromosomes.
for (i in seq_len(nrow(p))) {
  row <- p[i, ]
  for (config in names(comparison_specs)) {
    outfile <- file.path(
      "controls",
      paste0(row$chromosome, "_", config, ".ctl")
    )
    writeLines(make_ctl(row, config), outfile)
  }
}

message("Generated 18 species-level ajan vs alas gdi controls in gdi/controls/")
message("Generated reciprocal aab/abb imap files in gdi/imap/")
message("Simulation tree directory: gdi/output/trees/")
