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
if (length(file_arg) != 1) stop("Could not determine path to run_gdi.R")
script_path <- normalizePath(sub("^--file=", "", file_arg))
intg_dir <- dirname(dirname(script_path))
setwd(intg_dir)

bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
keep_trees <- tolower(Sys.getenv("KEEP_TREES", unset = "false")) %in% c("1", "true", "yes", "y")

expected_controls <- unlist(lapply(paste0("ch", 1:9), function(chrom) {
  c(
    file.path("controls", paste0(chrom, "_intg_aab.ctl")),
    file.path("controls", paste0(chrom, "_intg_abb.ctl")),
    file.path("controls", paste0(chrom, "_ca_hook_aab.ctl")),
    file.path("controls", paste0(chrom, "_ca_hook_abb.ctl"))
  )
}))
if (!all(file.exists(expected_controls))) {
  status <- system2("Rscript", file.path("scripts", "generate_controls.R"))
  if (status != 0) stop("Control generation failed")
}

bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) stop("BPP executable could not be run. Set BPP_BIN to the full path to bpp.")

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
  focal_pair_first <- !is.na(t1_left)

  unequal <- focal_pair_first & abs(t1_left - t1_right) > 1e-10
  if (any(unequal, na.rm = TRUE)) {
    stop("Unexpected unequal terminal branch lengths in focal cherries: ", fin)
  }

  success <- focal_pair_first & t1_left < tau

  tibble(
    n_gene_trees = length(lines),
    n_focal_pair_first = sum(focal_pair_first),
    n_focal_pair_first_before_tau = sum(success),
    gdi = sum(success) / length(lines)
  )
}

extract_tau <- function(ctl_file, node) {
  txt <- paste(read_lines(ctl_file), collapse = " ")
  pattern <- paste0(node, ":([0-9.eE+-]+)")
  x <- str_match(txt, pattern)[, 2]
  if (is.na(x)) stop("Could not find tau for node ", node, " in ", ctl_file)
  as.numeric(x)
}

comparison_specs <- tibble::tribble(
  ~comparison, ~config, ~focal_population, ~tau_node,
  "intg_N_vs_CA", "aab", "intg_nGL_Nslope", "IIH",
  "intg_N_vs_CA", "abb", "intg_CAswGL",      "IIH",
  "CA_vs_hook",   "aab", "intg_CAswGL",      "IH",
  "CA_vs_hook",   "abb", "hook",              "IH"
)

control_prefix <- function(comparison) {
  if (comparison == "intg_N_vs_CA") "intg" else "ca_hook"
}

run_one <- function(chrom, comparison, config, focal_population, tau_node) {
  prefix <- control_prefix(comparison)
  ctl <- file.path("controls", paste0(chrom, "_", prefix, "_", config, ".ctl"))
  treefile <- file.path("output", "trees", paste0(chrom, "_", prefix, "_", config, ".tree.txt"))
  dir.create(dirname(treefile), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(treefile)) unlink(treefile)
  tau <- extract_tau(ctl, tau_node)

  message(
    "Running ", chrom, " ", comparison, " ", config,
    " (focal = ", focal_population, ", tau_", tau_node, " = ", tau, ")"
  )

  status <- system2(bpp, c("--quiet", "--simulate", ctl))
  if (status != 0) stop("BPP failed for ", chrom, " ", comparison, " ", config)
  if (!file.exists(treefile)) stop("BPP did not create expected tree file: ", treefile)

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

chromosomes <- paste0("ch", 1:9)
raw <- map_dfr(chromosomes, function(chrom) {
  pmap_dfr(comparison_specs, function(comparison, config, focal_population, tau_node) {
    run_one(chrom, comparison, config, focal_population, tau_node)
  })
})

dir.create("output", showWarnings = FALSE)
write_csv(raw, file.path("output", "gdi_intg_long.csv"))

wide <- raw |>
  mutate(column = case_when(
    comparison == "intg_N_vs_CA" & focal_population == "intg_nGL_Nslope" ~ "gdi_intg_N_vs_CA__intg_nGL_Nslope",
    comparison == "intg_N_vs_CA" & focal_population == "intg_CAswGL" ~ "gdi_intg_N_vs_CA__intg_CAswGL",
    comparison == "CA_vs_hook" & focal_population == "intg_CAswGL" ~ "gdi_CA_vs_hook__intg_CAswGL",
    comparison == "CA_vs_hook" & focal_population == "hook" ~ "gdi_CA_vs_hook__hook",
    TRUE ~ NA_character_
  )) |>
  select(chromosome, column, gdi) |>
  pivot_wider(names_from = column, values_from = gdi) |>
  arrange(factor(chromosome, levels = chromosomes))

write_csv(wide, file.path("output", "gdi_intg.csv"))
print(wide)
