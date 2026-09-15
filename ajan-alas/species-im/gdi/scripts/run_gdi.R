# Run reciprocal BPP gene-tree simulations and calculate chromosome-specific gdi
# for the direct species-level D. ajanensis vs D. alaskensis comparison.
#
# The demographic parameters come from empirical chromosome-specific two-species
# IM fits. For each chromosome:
#
#   aab = 2 ajan + 1 alas -> gdi_ajan
#   abb = 1 ajan + 2 alas -> gdi_alas
#
# gdi is the proportion of simulated gene trees in which the two sequences from
# the focal population coalesce with each other before the species divergence
# time tau_R.
#
# Required environment variable if bpp is not on PATH:
#   BPP_BIN=/full/path/to/bpp
#
# Optional:
#   KEEP_TREES=true
# By default the 1,000,000-tree simulation files are deleted after scoring.
#
# Usage from anywhere in the repository:
#   BPP_BIN=/full/path/to/bpp \
#     Rscript ajan-alas/species-im/gdi/scripts/run_gdi.R

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

# Locate the gdi directory from this script's path so all relative paths work
# regardless of the directory from which Rscript is called.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine path to run_gdi.R")
script_path <- normalizePath(sub("^--file=", "", file_arg))
gdi_dir <- dirname(dirname(script_path))
setwd(gdi_dir)

# Read runtime options. KEEP_TREES is false by default because each simulation
# creates a very large one-million-tree text file and only its gdi score is needed.
bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
keep_trees <- tolower(Sys.getenv("KEEP_TREES", unset = "false")) %in%
  c("1", "true", "yes", "y")

# Ensure all 18 chromosome/configuration controls exist. If not, regenerate them
# automatically from the empirical species-level posterior means.
chromosomes <- paste0("ch", 1:9)
expected_controls <- unlist(lapply(chromosomes, function(chrom) {
  c(
    file.path("controls", paste0(chrom, "_aab.ctl")),
    file.path("controls", paste0(chrom, "_abb.ctl"))
  )
}))
if (!all(file.exists(expected_controls))) {
  message("One or more gdi controls are missing; generating controls first")
  status <- system2("Rscript", file.path("scripts", "generate_controls.R"))
  if (status != 0) stop("Control generation failed")
}

# Confirm that the requested BPP executable is callable before beginning the
# simulations. This avoids discovering a PATH problem after the workflow starts.
bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) {
  stop("BPP executable could not be run. Set BPP_BIN to the full path to bpp.")
}

# Score one BPP gene-tree file. BPP writes labels such as ajan^ajan1; remove the
# simulated individual suffix first so the focal cherry can be detected using
# only the population name. A successful gdi event requires that the two focal
# sequences form the first coalescence and that this coalescence occurs < tau_R.
gdi_from_gtree <- function(fin, focal_population, tau) {
  lines <- read_lines(fin, progress = FALSE) |>
    str_replace_all("\\^[A-Za-z0-9_]+", "")

  number <- "[0-9.eE+-]+"
  pattern <- paste0(
    "\\(", focal_population, ":(", number, "),",
    focal_population, ":(", number, ")\\)"
  )

  m <- str_match(lines, pattern)
  t_left <- suppressWarnings(as.numeric(m[, 2]))
  t_right <- suppressWarnings(as.numeric(m[, 3]))
  focal_pair_first <- !is.na(t_left)

  # In a rooted three-tip gene tree, the two branches descending from a focal
  # cherry should have identical lengths. Flag malformed/unexpected tree output.
  unequal <- focal_pair_first & abs(t_left - t_right) > 1e-10
  if (any(unequal, na.rm = TRUE)) {
    stop("Unexpected unequal terminal branch lengths in focal cherries: ", fin)
  }

  success <- focal_pair_first & t_left < tau

  tibble(
    n_gene_trees = length(lines),
    n_focal_pair_first = sum(focal_pair_first),
    n_focal_pair_first_before_tau = sum(success),
    gdi = sum(success) / length(lines)
  )
}

# Extract tau_R directly from the generated simulation control. This guarantees
# that the scoring threshold is exactly the divergence time supplied to BPP.
extract_tau_R <- function(ctl_file) {
  txt <- paste(read_lines(ctl_file), collapse = " ")
  x <- str_match(txt, "\\)R:([0-9.eE+-]+)")[, 2]
  if (is.na(x)) stop("Could not find tau_R in ", ctl_file)
  as.numeric(x)
}

# Define the reciprocal configurations and which lineage each one estimates.
comparison_specs <- tibble::tribble(
  ~configuration, ~focal_population,
  "aab",          "ajan",
  "abb",          "alas"
)

# Run and score one chromosome/configuration simulation. Tree files are removed
# after scoring unless KEEP_TREES=true was explicitly requested.
run_one <- function(chrom, configuration, focal_population) {
  ctl <- file.path("controls", paste0(chrom, "_", configuration, ".ctl"))
  treefile <- file.path(
    "output", "trees",
    paste0(chrom, "_", configuration, ".tree.txt")
  )
  dir.create(dirname(treefile), recursive = TRUE, showWarnings = FALSE)

  # Remove a stale tree file so a failed simulation cannot accidentally be scored
  # as if it came from the current run.
  if (file.exists(treefile)) unlink(treefile)
  tau <- extract_tau_R(ctl)

  message(
    "Running ", chrom, " ", configuration,
    " (focal = ", focal_population, ", tau_R = ", tau, ")"
  )

  status <- system2(bpp, c("--quiet", "--simulate", ctl))
  if (status != 0) stop("BPP failed for ", chrom, " ", configuration)
  if (!file.exists(treefile)) {
    stop("BPP did not create expected tree file: ", treefile)
  }

  score <- gdi_from_gtree(treefile, focal_population, tau)
  if (!keep_trees) unlink(treefile)

  tibble(
    chromosome = chrom,
    configuration = configuration,
    focal_population = focal_population,
    tau_R = tau
  ) |>
    bind_cols(score)
}

# Run both reciprocal simulations for each chromosome. Simulations are executed
# sequentially because each one is fast relative to the empirical BPP fits and
# temporarily writes a large tree file to disk.
raw <- map_dfr(chromosomes, function(chrom) {
  pmap_dfr(
    comparison_specs,
    function(configuration, focal_population) {
      run_one(chrom, configuration, focal_population)
    }
  )
})

# Save a long-format audit table containing simulation counts, tau, and gdi.
dir.create("output", recursive = TRUE, showWarnings = FALSE)
write_csv(raw, file.path("output", "gdi_ajan_alas_long.csv"))

# Produce the concise chromosome-by-chromosome table used for interpretation and
# plotting. Each row contains directional gdi for both nominal species.
wide <- raw |>
  mutate(column = case_when(
    focal_population == "ajan" ~ "gdi_ajan",
    focal_population == "alas" ~ "gdi_alas",
    TRUE ~ NA_character_
  )) |>
  select(chromosome, column, gdi) |>
  pivot_wider(names_from = column, values_from = gdi) |>
  arrange(factor(chromosome, levels = chromosomes))

write_csv(wide, file.path("output", "gdi_ajan_alas.csv"))

# Summarize chromosome-level variation without replacing the chromosome-specific
# results. These statistics are convenient for manuscript reporting.
summary_table <- wide |>
  summarise(
    mean_gdi_ajan = mean(gdi_ajan),
    min_gdi_ajan = min(gdi_ajan),
    max_gdi_ajan = max(gdi_ajan),
    mean_gdi_alas = mean(gdi_alas),
    min_gdi_alas = min(gdi_alas),
    max_gdi_alas = max(gdi_alas)
  )
write_csv(summary_table, file.path("output", "gdi_ajan_alas_summary.csv"))

print(wide)
print(summary_table)
message("Species-level ajan vs alas gdi analysis completed successfully")
