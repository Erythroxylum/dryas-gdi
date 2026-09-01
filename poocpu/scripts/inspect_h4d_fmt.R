suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript poocpu/scripts/inspect_h4d_fmt.R /path/to/h4d-s47-p9-fmt.csv")

x <- read_csv(args[1], show_col_types = FALSE)
cat("Rows:", nrow(x), " Columns:", ncol(x), "\n\n")
cat("Column names:\n")
cat(paste(names(x), collapse = "\n"), "\n\n")
cat("First 20 rows:\n")
print(utils::head(x, 20), width = Inf)
