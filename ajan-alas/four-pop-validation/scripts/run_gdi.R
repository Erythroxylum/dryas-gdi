# Run BPP gene-tree simulations and calculate Eq. 13 gdi for the four-population
# ajanensis/alaskensis validation analysis.
#
# The fitted four-population demographic model is retained in every simulation,
# but only three focal sequences are sampled at a time. This makes the scoring
# directly comparable to the intg and PoOcPu1 workflows in this repository.

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
})

# Locate this analysis directory from the script path.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Load the model specification and BPP executable settings.
source("model/four_pop_model.R", local = globalenv())
bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
keep_trees <- tolower(Sys.getenv("KEEP_TREES", unset = "false")) %in% c("1", "true", "yes", "y")

# Generate controls automatically if they are not present yet.
expected <- unlist(lapply(paste0("ch", 1:9), function(chrom) {
  paste0(
    "controls/", chrom, "_",
    rep(c("ajan_I_vs_S", "alas_I_vs_S"), each = 2), "_",
    rep(c("aab", "abb"), 2), ".ctl"
  )
}))
if (!all(file.exists(expected))) {
  status <- system2("Rscript", "scripts/generate_controls.R")
  if (status != 0) stop("Control generation failed")
}

# Confirm that BPP can be executed before starting the 36-run batch.
bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) stop("BPP could not be run. Set BPP_BIN to the full path to bpp.")

# Score one simulated gene-tree file using Kornai et al. Eq. 13.
# A success requires the two focal-population samples to form the first cherry
# and their coalescence time (the equal terminal branch length) to be < tau.
gdi_from_gtree <- function(fin, focal_population, tau) {
  lines <- read_lines(fin, progress = FALSE) |>
    str_replace_all("\\^[A-Za-z0-9_]+", "")

  number <- "[0-9.eE+-]+"
  pattern <- paste0(
    "\\(", focal_population, ":(", number, "),",
    focal_population, ":(", number, ")\\)"
  )
  m <- str_match(lines, pattern)
  t1_left <- suppressWarnings(as.numeric(m[, 2]))
  t1_right <- suppressWarnings(as.numeric(m[, 3]))
  cherry <- !is.na(t1_left)

  # BPP ultrametric trees should give equal terminal lengths for the focal cherry.
  unequal <- cherry & abs(t1_left - t1_right) > 1e-10
  if (any(unequal, na.rm = TRUE)) stop("Unexpected unequal focal terminal branches in ", fin)

  success <- cherry & t1_left < tau

  tibble(
    n_gene_trees = length(lines),
    n_focal_pair_first = sum(cherry),
    n_focal_pair_first_before_tau = sum(success),
    gdi = sum(success) / length(lines)
  )
}

# Pull the named divergence-time cutoff from a generated numeric control file.
extract_tau <- function(ctl_file, node) {
  txt <- paste(read_lines(ctl_file), collapse = " ")
  x <- str_match(txt, paste0(node, ":([0-9.eE+-]+)"))[, 2]
  if (is.na(x)) stop("Could not find tau_", node, " in ", ctl_file)
  as.numeric(x)
}

# Run one chromosome/configuration and immediately score its simulated trees.
run_one <- function(chrom, comparison, config, pop_a, pop_b, focal_population, tau_node) {
  ctl <- file.path("controls", paste0(chrom, "_", comparison, "_", config, ".ctl"))
  treefile <- file.path("output", "trees", paste0(chrom, "_", comparison, "_", config, ".tree.txt"))
  dir.create(dirname(treefile), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(treefile)) unlink(treefile)

  tau <- extract_tau(ctl, tau_node)
  message("Running ", chrom, " ", comparison, " ", config,
          " focal=", focal_population, " tau_", tau_node, "=", tau)

  status <- system2(bpp, c("--quiet", "--simulate", ctl))
  if (status != 0) stop("BPP failed for ", chrom, " ", comparison, " ", config)
  if (!file.exists(treefile)) stop("Expected tree file was not created: ", treefile)

  score <- gdi_from_gtree(treefile, focal_population, tau)
  if (!keep_trees) unlink(treefile)

  tibble(
    chromosome = chrom,
    comparison = comparison,
    configuration = config,
    focal_population = focal_population,
    tau_node = tau_node,
    tau = tau
  ) |>
    bind_cols(score)
}

# Run all four directional tests for all nine chromosomes.
chromosomes <- paste0("ch", 1:9)
raw <- map_dfr(chromosomes, function(chrom) {
  pmap_dfr(comparison_specs, function(comparison, config, pop_a, pop_b, focal_population, tau_node) {
    run_one(chrom, comparison, config, pop_a, pop_b, focal_population, tau_node)
  })
})

# Save both a detailed long table and a one-row-per-chromosome summary.
dir.create("output", showWarnings = FALSE)
write_csv(raw, "output/gdi_four_pop_long.csv")
wide <- raw |>
  mutate(column = paste0("gdi_", comparison, "__", focal_population)) |>
  select(chromosome, column, gdi) |>
  pivot_wider(names_from = column, values_from = gdi)
write_csv(wide, "output/gdi_four_pop.csv")
print(wide)
