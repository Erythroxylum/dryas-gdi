suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
intg_dir <- dirname(dirname(script_path))
setwd(intg_dir)

p <- read_csv("parameters/aih-prior3-s16-p4-means.csv", show_col_types = FALSE)
dir.create("controls", showWarnings = FALSE)
dir.create(file.path("output", "trees"), recursive = TRUE, showWarnings = FALSE)

migration_pairs <- list(
  c("ajan_alas", "intg_nGL_Nslope"), c("ajan_alas", "intg_CAswGL"),
  c("ajan_alas", "hook"), c("ajan_alas", "IIH"), c("ajan_alas", "IH"),
  c("intg_nGL_Nslope", "ajan_alas"), c("intg_nGL_Nslope", "intg_CAswGL"),
  c("intg_nGL_Nslope", "hook"), c("intg_nGL_Nslope", "IH"),
  c("intg_CAswGL", "ajan_alas"), c("intg_CAswGL", "intg_nGL_Nslope"),
  c("intg_CAswGL", "hook"), c("hook", "ajan_alas"),
  c("hook", "intg_nGL_Nslope"), c("hook", "intg_CAswGL"),
  c("IIH", "ajan_alas"), c("IH", "ajan_alas"), c("IH", "intg_nGL_Nslope")
)

comparison_specs <- list(
  intg = list(
    aab = list(counts = "0 2 1 0", imap = "imap/aab.imap.txt"),
    abb = list(counts = "0 1 2 0", imap = "imap/abb.imap.txt")
  ),
  ca_hook = list(
    aab = list(counts = "0 0 2 1", imap = "imap/ca_hook_aab.imap.txt"),
    abb = list(counts = "0 0 1 2", imap = "imap/ca_hook_abb.imap.txt")
  )
)

make_ctl <- function(row, comparison, config) {
  spec <- comparison_specs[[comparison]][[config]]
  treefile <- sprintf("output/trees/%s_%s_%s.tree.txt", row$chromosome, comparison, config)

  tree <- sprintf(
    "(ajan_alas #%g, (intg_nGL_Nslope #%g, (intg_CAswGL #%g, hook #%g)IH:%g #%g)IIH:%g #%g)R:%g #%g;",
    row$theta_ajan_alas, row$theta_intg_nGL_Nslope, row$theta_intg_CAswGL,
    row$theta_hook, row$tau_IH, row$theta_IH, row$tau_IIH, row$theta_IIH,
    row$tau_R, row$theta_R
  )

  mig <- vapply(migration_pairs, function(pair) {
    nm <- paste0("W_", pair[1], "_to_", pair[2])
    sprintf("            %s %s %.10g", pair[1], pair[2], row[[nm]])
  }, character(1))

  paste(c(
    "seed = -1",
    paste("treefile =", treefile),
    paste("Imapfile =", spec$imap),
    "species&tree = 4 ajan_alas intg_nGL_Nslope intg_CAswGL hook",
    paste("                ", spec$counts),
    tree,
    "loci&length = 1000000 50",
    "migration = 18",
    mig
  ), collapse = "\n")
}

for (i in seq_len(nrow(p))) {
  row <- p[i, ]
  for (comparison in names(comparison_specs)) {
    for (config in c("aab", "abb")) {
      outfile <- file.path("controls", paste0(row$chromosome, "_", comparison, "_", config, ".ctl"))
      writeLines(make_ctl(row, comparison, config), outfile)
    }
  }
}

message("Generated 36 controls in intg/controls/ and ensured intg/output/trees/ exists")
