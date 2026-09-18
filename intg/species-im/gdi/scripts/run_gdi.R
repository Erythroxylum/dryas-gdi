# Run reciprocal BPP simulations and calculate chromosome-specific gdi for the
# two geographic D. integrifolia populations under two fitted IM parameter sets:
#
#   prior2  = W prior mean 2
#   prior20 = W prior mean 20
#
# Required if bpp is not on PATH:
#   BPP_BIN=/full/path/to/bpp
#
# Optional:
#   KEEP_TREES=true
#
# Usage:
#   BPP_BIN=/full/path/to/bpp Rscript intg/species-im/gdi/scripts/run_gdi.R

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
gdi_dir <- dirname(dirname(script_path))
setwd(gdi_dir)

bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
keep_trees <- tolower(Sys.getenv("KEEP_TREES", unset = "false")) %in%
  c("1", "true", "yes", "y")

scenarios <- c("prior2", "prior20")
chromosomes <- paste0("ch", 1:9)

expected_controls <- unlist(lapply(scenarios, function(scenario) {
  unlist(lapply(chromosomes, function(chrom) {
    c(
      file.path(scenario, "controls", paste0(chrom, "_aab.ctl")),
      file.path(scenario, "controls", paste0(chrom, "_abb.ctl"))
    )
  }))
}))
if (!all(file.exists(expected_controls))) {
  message("One or more controls are missing; generating them first")
  rc <- system2("Rscript", file.path("scripts", "generate_controls.R"))
  if (rc != 0) stop("Control generation failed")
}

bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
rc <- attr(bpp_test, "status")
if (!is.null(rc) && rc != 0) {
  stop("BPP executable could not be run. Set BPP_BIN to the full path to bpp.")
}

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

extract_tau_R <- function(ctl_file) {
  txt <- paste(read_lines(ctl_file), collapse = " ")
  x <- str_match(txt, "\\)R:([0-9.eE+-]+)")[, 2]
  if (is.na(x)) stop("Could not find tau_R in ", ctl_file)
  as.numeric(x)
}

comparison_specs <- tibble::tribble(
  ~configuration, ~focal_population,
  "aab",          "intg_nGL_Nslope",
  "abb",          "intg_CAswGL"
)

run_one <- function(scenario, chrom, configuration, focal_population) {
  scenario_dir <- file.path(scenario)
  ctl <- file.path(scenario_dir, "controls", paste0(chrom, "_", configuration, ".ctl"))
  treefile <- file.path(
    scenario_dir, "output", "trees",
    paste0(chrom, "_", configuration, ".tree.txt")
  )
  dir.create(dirname(treefile), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(treefile)) unlink(treefile)
  tau <- extract_tau_R(ctl)

  message(
    "Running ", scenario, " ", chrom, " ", configuration,
    " (focal = ", focal_population, ", tau_R = ", tau, ")"
  )

  oldwd <- getwd()
  setwd(scenario_dir)
  on.exit(setwd(oldwd), add = TRUE)
  ctl_local <- file.path("controls", paste0(chrom, "_", configuration, ".ctl"))
  rc <- system2(bpp, c("--quiet", "--simulate", ctl_local))
  setwd(oldwd)
  on.exit(NULL, add = FALSE)

  if (rc != 0) stop("BPP failed for ", scenario, " ", chrom, " ", configuration)
  if (!file.exists(treefile)) {
    stop("BPP did not create expected tree file: ", treefile)
  }

  score <- gdi_from_gtree(treefile, focal_population, tau)
  if (!keep_trees) unlink(treefile)

  tibble(
    scenario = scenario,
    chromosome = chrom,
    configuration = configuration,
    focal_population = focal_population,
    tau_R = tau
  ) |>
    bind_cols(score)
}

raw <- map_dfr(scenarios, function(scenario) {
  map_dfr(chromosomes, function(chrom) {
    pmap_dfr(
      comparison_specs,
      function(configuration, focal_population) {
        run_one(scenario, chrom, configuration, focal_population)
      }
    )
  })
})

dir.create("output", recursive = TRUE, showWarnings = FALSE)
write_csv(raw, file.path("output", "gdi_integrifolia_prior_sensitivity_long.csv"))

wide <- raw |>
  mutate(column = case_when(
    focal_population == "intg_nGL_Nslope" ~ "gdi_nGL_Nslope",
    focal_population == "intg_CAswGL" ~ "gdi_CAswGL",
    TRUE ~ NA_character_
  )) |>
  select(scenario, chromosome, column, gdi) |>
  pivot_wider(names_from = column, values_from = gdi) |>
  arrange(
    factor(scenario, levels = scenarios),
    factor(chromosome, levels = chromosomes)
  )

write_csv(wide, file.path("output", "gdi_integrifolia_prior_sensitivity.csv"))

summary_table <- wide |>
  group_by(scenario) |>
  summarise(
    mean_gdi_nGL_Nslope = mean(gdi_nGL_Nslope),
    min_gdi_nGL_Nslope = min(gdi_nGL_Nslope),
    max_gdi_nGL_Nslope = max(gdi_nGL_Nslope),
    mean_gdi_CAswGL = mean(gdi_CAswGL),
    min_gdi_CAswGL = min(gdi_CAswGL),
    max_gdi_CAswGL = max(gdi_CAswGL),
    .groups = "drop"
  )

write_csv(summary_table, file.path("output", "gdi_integrifolia_prior_sensitivity_summary.csv"))

print(wide)
print(summary_table)
message("Direct integrifolia gdi prior-sensitivity analysis completed successfully")
