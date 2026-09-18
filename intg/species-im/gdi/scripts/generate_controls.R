# Generate reciprocal three-sequence BPP simulation controls for the direct
# two-population D. integrifolia gdi analysis under two fitted parameterizations:
#
#   prior2  = original IM fit with W ~ Gamma(2,1), prior mean 2
#   prior20 = sensitivity fit with W ~ Gamma(2,0.1), prior mean 20
#
# Each scenario uses chromosome-specific posterior means from its own empirical
# IM fit. For every chromosome:
#   aab = 2 intg_nGL_Nslope + 1 intg_CAswGL
#   abb = 1 intg_nGL_Nslope + 2 intg_CAswGL
#
# Usage:
#   Rscript intg/species-im/gdi/scripts/generate_controls.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
gdi_dir <- dirname(dirname(script_path))
setwd(gdi_dir)

parameter_files <- c(
  prior2 = file.path("..", "fit", "species_im_means.csv"),
  prior20 = file.path("..", "sensitivity", "prior20", "species_im_means.csv")
)

required_columns <- c(
  "chromosome",
  "theta_intg_nGL_Nslope",
  "theta_intg_CAswGL",
  "theta_R",
  "tau_R",
  "W_intg_nGL_Nslope_to_intg_CAswGL",
  "W_intg_CAswGL_to_intg_nGL_Nslope"
)

chromosomes <- paste0("ch", 1:9)

dir.create("imap", recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    "N1\tintg_nGL_Nslope",
    "N2\tintg_nGL_Nslope",
    "CA1\tintg_CAswGL"
  ),
  file.path("imap", "intg_aab.imap.txt")
)
writeLines(
  c(
    "N1\tintg_nGL_Nslope",
    "CA1\tintg_CAswGL",
    "CA2\tintg_CAswGL"
  ),
  file.path("imap", "intg_abb.imap.txt")
)

comparison_specs <- list(
  aab = list(counts = "2 1", imap = "../../imap/intg_aab.imap.txt"),
  abb = list(counts = "1 2", imap = "../../imap/intg_abb.imap.txt")
)

for (scenario in names(parameter_files)) {
  parameter_file <- parameter_files[[scenario]]
  if (!file.exists(parameter_file)) {
    stop("Missing parameter file for ", scenario, ": ", parameter_file)
  }

  p <- read_csv(parameter_file, show_col_types = FALSE)
  missing <- setdiff(required_columns, names(p))
  if (length(missing)) {
    stop("Parameter file ", parameter_file, " is missing: ", paste(missing, collapse = ", "))
  }
  if (!setequal(p$chromosome, chromosomes) || nrow(p) != 9L) {
    stop("Expected exactly one row for each chromosome ch1-ch9 in ", parameter_file)
  }

  p <- p |>
    mutate(chromosome = factor(chromosome, levels = chromosomes)) |>
    arrange(chromosome) |>
    mutate(chromosome = as.character(chromosome))

  scenario_dir <- file.path(scenario)
  dir.create(file.path(scenario_dir, "controls"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(scenario_dir, "output", "trees"), recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(p))) {
    row <- p[i, ]

    for (config in names(comparison_specs)) {
      spec <- comparison_specs[[config]]
      treefile <- sprintf(
        "output/trees/%s_%s.tree.txt",
        row$chromosome, config
      )

      tree <- sprintf(
        "(intg_nGL_Nslope #%.10g, intg_CAswGL #%.10g)R:%.10g #%.10g;",
        row$theta_intg_nGL_Nslope,
        row$theta_intg_CAswGL,
        row$tau_R,
        row$theta_R
      )

      migration <- c(
        sprintf(
          "            intg_nGL_Nslope intg_CAswGL %.10g",
          row$W_intg_nGL_Nslope_to_intg_CAswGL
        ),
        sprintf(
          "            intg_CAswGL intg_nGL_Nslope %.10g",
          row$W_intg_CAswGL_to_intg_nGL_Nslope
        )
      )

      ctl <- paste(c(
        paste0("# Direct integrifolia gdi simulation: ", scenario, " ", row$chromosome, " ", config),
        paste0("# Parameter source: ", parameter_file),
        "# aab estimates gdi for intg_nGL_Nslope; abb estimates gdi for intg_CAswGL.",
        "seed = -1",
        paste("treefile =", treefile),
        paste("Imapfile =", spec$imap),
        "species&tree = 2 intg_nGL_Nslope intg_CAswGL",
        paste("                ", spec$counts),
        tree,
        "loci&length = 1000000 50",
        "migration = 2",
        migration
      ), collapse = "\n")

      outfile <- file.path(
        scenario_dir, "controls",
        paste0(row$chromosome, "_", config, ".ctl")
      )
      writeLines(ctl, outfile)
    }
  }

  message("Generated 18 gdi controls for ", scenario)
}

message("Generated direct integrifolia gdi controls for prior2 and prior20")
