suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript poocpu/scripts/extract_parameters.R /path/to/h4d-s47-p9-fmt.csv")
}

infile <- normalizePath(args[1], mustWork = TRUE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

dir.create("parameters", showWarnings = FALSE)
outfile <- file.path("parameters", "h4d-s47-p9-means.csv")

x <- read_csv(infile, show_col_types = FALSE)
chromosomes <- paste0("ch", 1:9)
if (!all(c("variable", chromosomes) %in% names(x))) {
  stop("Expected columns variable, ch1, ..., ch9 in H4d formatted CSV")
}

# Extract the posterior mean, i.e. the number before the parenthesized interval.
mean_value <- function(z) {
  suppressWarnings(as.numeric(sub("\\s*\\(.*$", "", trimws(z))))
}

# Map the reduced PoOcPu1 model onto the corresponding H4d nodes/populations.
# Important mappings from Dawson's PoOcPu1 model diagram:
#   Po      = H4d octo_Pobeda (tip 7)
#   Pu      = H4d punc        (tip 9)
#   RU_SJ   = H4d octo_RU_SJ  (tip 5)
#   R       = H4d S            (node 13), NOT the global H4d root R
#   AOcPu   = H4d node 16
#   Octo    = H4d node 17
#   Punc    = H4d node 19
# Migration edges retained in PoOcPu1 are H4d 4->5, 5->4, and 7->19.
parameter_map <- c(
  tau_Octo = "tau:17:Octo",
  tau_Punc = "tau:19:Punc",
  tau_AOcPu = "tau:16:AOcPu",
  tau_R = "tau:13:S",

  theta_octo_EU = "theta:4:octo_EU",
  theta_octo_Carp_MK = "theta:6:octo_Carp_MK",
  theta_RU_SJ = "theta:5:octo_RU_SJ",
  theta_Pu = "theta:9:punc",
  theta_Po = "theta:7:octo_Pobeda",
  theta_Octo = "theta:17:Octo",
  theta_Punc = "theta:19:Punc",
  theta_AOcPu = "theta:16:AOcPu",
  theta_R = "theta:13:S",

  W_EU_to_RUSJ = "W:4->5:octo_EU->octo_RU_SJ",
  W_RUSJ_to_EU = "W:5->4:octo_RU_SJ->octo_EU",
  W_Po_to_Punc = "W:7->19:octo_Pobeda->Punc"
)

missing_vars <- setdiff(unname(parameter_map), x$variable)
if (length(missing_vars)) {
  stop("Could not find these H4d variables: ", paste(missing_vars, collapse = ", "))
}

out <- tibble(chromosome = chromosomes)
for (new_name in names(parameter_map)) {
  old_name <- parameter_map[[new_name]]
  vals <- x[x$variable == old_name, chromosomes, drop = FALSE]
  if (nrow(vals) != 1) stop("Expected exactly one row for ", old_name)
  out[[new_name]] <- vapply(vals[1, ], mean_value, numeric(1))
}

if (any(!is.finite(as.matrix(out[, -1])))) {
  stop("At least one parameter mean could not be parsed")
}

write_csv(out, outfile)
message("Wrote ", outfile)
print(out, n = Inf, width = Inf)
