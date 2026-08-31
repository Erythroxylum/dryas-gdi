suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(tibble)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine path to run_gdi.R")
script_path <- normalizePath(sub("^--file=", "", file_arg))
intg_dir <- dirname(dirname(script_path))
setwd(intg_dir)

bpp <- Sys.getenv("BPP_BIN", unset = "bpp")
keep_trees <- tolower(Sys.getenv("KEEP_TREES", unset = "false")) %in% c("1", "true", "yes", "y")

if (!all(file.exists(file.path("controls", paste0("ch", rep(1:9, each = 2), "_", rep(c("aab", "abb"), 9), ".ctl"))))) {
  system2("Rscript", file.path("scripts", "generate_controls.R"))
}

bpp_test <- suppressWarnings(system2(bpp, "--version", stdout = TRUE, stderr = TRUE))
status <- attr(bpp_test, "status")
if (!is.null(status) && status != 0) stop("BPP executable could not be run. Set BPP_BIN to the full path to bpp.")

gdi_from_gtree <- function(fin, focal_population, tau) {
  lines <- read_lines(fin, progress = FALSE) |> str_replace_all("\\^[A-Za-z0-9_]+", "")
  number <- "[0-9.eE+-]+"
  pattern <- paste0("\\(", focal_population, ":(", number, "),", focal_population, ":(", number, ")\\)")
  m <- str_match(lines, pattern)
  t1 <- suppressWarnings(as.numeric(m[, 2]))
  focal_pair_first <- !is.na(t1)
  success <- focal_pair_first & t1 < tau
  tibble(
    n_gene_trees = length(lines),
    n_focal_pair_first = sum(focal_pair_first),
    n_focal_pair_first_before_tau = sum(success),
    gdi = sum(success) / length(lines)
  )
}

extract_tau <- function(ctl_file, node = "IIH") {
  txt <- paste(read_lines(ctl_file), collapse = " ")
  pattern <- paste0(node, ":([0-9.eE+-]+)")
  x <- str_match(txt, pattern)[, 2]
  if (is.na(x)) stop("Could not find tau for node ", node, " in ", ctl_file)
  as.numeric(x)
}

run_one <- function(chrom, config) {
  ctl <- file.path("controls", paste0(chrom, "_", config, ".ctl"))
  treefile <- file.path("output", "trees", paste0(chrom, "_", config, ".tree.txt"))
  dir.create(dirname(treefile), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(treefile)) unlink(treefile)
  focal <- if (config == "aab") "intg_nGL_Nslope" else "intg_CAswGL"
  tau <- extract_tau(ctl, "IIH")
  message("Running ", chrom, " ", config, " (focal = ", focal, ", tau_IIH = ", tau, ")")
  status <- system2(bpp, c("--quiet", "--simulate", ctl))
  if (status != 0) stop("BPP failed for ", chrom, " ", config)
  if (!file.exists(treefile)) stop("BPP did not create expected tree file: ", treefile)
  score <- gdi_from_gtree(treefile, focal, tau)
  if (!keep_trees) unlink(treefile)
  tibble(chromosome = chrom, configuration = config, focal_population = focal, tau_IIH = tau) |> bind_cols(score)
}

chromosomes <- paste0("ch", 1:9)
raw <- map_dfr(chromosomes, function(chrom) bind_rows(run_one(chrom, "aab"), run_one(chrom, "abb")))
dir.create("output", showWarnings = FALSE)
write_csv(raw, file.path("output", "gdi_intg_long.csv"))
wide <- raw |>
  select(chromosome, focal_population, gdi) |>
  tidyr::pivot_wider(names_from = focal_population, values_from = gdi) |>
  rename(gdi_intg_nGL_Nslope = intg_nGL_Nslope, gdi_intg_CAswGL = intg_CAswGL)
write_csv(wide, file.path("output", "gdi_intg.csv"))
print(wide)
