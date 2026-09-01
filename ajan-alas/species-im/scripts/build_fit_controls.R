# Build empirical BPP A00 controls for a two-lineage ajanensis vs alaskensis IM model.
#
# This does NOT invent a new dataset. Instead, it takes the existing four-population
# ajan/alaskensis empirical controls and imap, collapses Interior + Seward within
# each taxon, and preserves the original sequence files and MCMC settings.
#
# Usage:
#   Rscript scripts/build_fit_controls.R /path/to/four_pop_ctl_dir /path/to/four_pop_imap.txt
#
# The source control directory should contain one final empirical BPP control for
# each chromosome ch1-ch9. Prefer the W ~ G(2,0.1) / prior3 controls because
# Yuttapong found their root-age estimates more stable than prior2.

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# Locate this species-IM analysis directory from the script location.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Read the two required local source paths.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Usage: Rscript scripts/build_fit_controls.R ",
    "/path/to/four_pop_ctl_dir /path/to/four_pop_imap.txt"
  )
}
source_ctl_dir <- normalizePath(args[1], mustWork = TRUE)
source_imap <- normalizePath(args[2], mustWork = TRUE)

# Read the original imap as two whitespace-delimited columns: sample and population.
imap <- read_table2(
  source_imap,
  col_names = c("sample", "population"),
  col_types = cols(.default = col_character())
)

# Collapse geographic populations into the two species-level lineages.
collapsed <- imap |>
  mutate(population = case_when(
    population %in% c("ajan_Interior", "ajan_Seward") ~ "ajan",
    population %in% c("alas_Interior", "alas_Seward") ~ "alas",
    TRUE ~ NA_character_
  ))
if (any(is.na(collapsed$population))) {
  bad <- unique(imap$population[is.na(collapsed$population)])
  stop("Unexpected populations in source imap: ", paste(bad, collapse = ", "))
}

# Count samples automatically rather than assuming the source design.
counts <- table(factor(collapsed$population, levels = c("ajan", "alas")))
if (any(counts == 0)) stop("Both ajan and alas must have sampled sequences")

# Write the species-level imap used by every chromosome.
dir.create("fit/imap", recursive = TRUE, showWarnings = FALSE)
new_imap <- normalizePath("fit/imap/ajan_alas.imap.txt", mustWork = FALSE)
write.table(
  collapsed,
  file = new_imap,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = " "
)
new_imap <- normalizePath(new_imap, mustWork = TRUE)

# Helper: identify the unique source control corresponding to one chromosome.
# We intentionally fail if several candidate files match, because silently choosing
# the wrong BPP model would invalidate the fit.
find_ctl <- function(chrom) {
  files <- list.files(source_ctl_dir, pattern = "\\.ctl$", recursive = TRUE, full.names = TRUE)
  hit <- files[str_detect(basename(files), regex(paste0("(^|[^0-9])", chrom, "([^0-9]|$)"), ignore_case = TRUE))]
  if (length(hit) != 1) {
    stop(
      "Expected exactly one source control for ", chrom, "; found ", length(hit),
      if (length(hit)) paste0(": ", paste(hit, collapse = ", ")) else ""
    )
  }
  normalizePath(hit)
}

# Helper: make relative input-file paths in copied controls absolute so the new
# controls can run from ajan-alas/species-im/fit/controls without losing access
# to the original chromosome data.
absolutize_directive <- function(lines, key, source_dir) {
  idx <- grep(paste0("^\\s*", key, "\\s*="), lines, ignore.case = TRUE)
  if (!length(idx)) return(lines)
  for (i in idx) {
    value <- str_trim(sub("^[^=]+=[[:space:]]*", "", lines[i]))
    if (nzchar(value) && !str_detect(value, "^/")) {
      value <- normalizePath(file.path(source_dir, value), mustWork = FALSE)
      lines[i] <- paste(key, "=", value)
    }
  }
  lines
}

# Helper: replace the three-line species&tree block in a BPP control.
replace_species_block <- function(lines) {
  i <- grep("^\\s*species&tree\\s*=", lines, ignore.case = TRUE)
  if (length(i) != 1) stop("Expected exactly one species&tree block")
  if (i + 2 > length(lines)) stop("Incomplete species&tree block")

  replacement <- c(
    "species&tree = 2 ajan alas",
    paste("                ", paste(as.integer(counts), collapse = " ")),
    "(ajan, alas)R;"
  )
  append(lines[-c(i:(i + 2))], replacement, after = i - 1)
}

# Helper: replace the fitted four-pop migration graph with the two-direction IM graph.
# Existing migration priors elsewhere in the source control are intentionally retained.
replace_migration_block <- function(lines) {
  i <- grep("^\\s*migration\\s*=", lines, ignore.case = TRUE)
  if (length(i) != 1) stop("Expected exactly one migration block")

  n_old <- suppressWarnings(as.integer(str_extract(lines[i], "[0-9]+\\s*$")))
  if (is.na(n_old)) stop("Could not parse existing migration count")
  end <- i + n_old
  if (end > length(lines)) stop("Existing migration block is incomplete")

  replacement <- c(
    "migration = 2",
    "            ajan alas",
    "            alas ajan"
  )
  append(lines[-c(i:end)], replacement, after = i - 1)
}

# Create one two-population empirical control for each chromosome while preserving
# the source model's sequence file, mutation model, priors, burnin, sampling, etc.
dir.create("fit/controls", recursive = TRUE, showWarnings = FALSE)
dir.create("fit/output", recursive = TRUE, showWarnings = FALSE)

for (chrom in paste0("ch", 1:9)) {
  source_ctl <- find_ctl(chrom)
  source_dir <- dirname(source_ctl)
  lines <- read_lines(source_ctl)

  # Preserve access to chromosome data after the control is copied.
  for (key in c("seqfile", "locusratefile", "heredityfile")) {
    lines <- absolutize_directive(lines, key, source_dir)
  }

  # Point the fit to the collapsed species-level sample map.
  imap_idx <- grep("^\\s*Imapfile\\s*=", lines, ignore.case = TRUE)
  if (length(imap_idx) != 1) stop("Expected exactly one Imapfile directive in ", source_ctl)
  lines[imap_idx] <- paste("Imapfile =", new_imap)

  # Replace only the biological model components: species tree and migration graph.
  lines <- replace_species_block(lines)
  lines <- replace_migration_block(lines)

  # Keep each chromosome's BPP outputs separate and inside this repository workflow.
  mcmc_idx <- grep("^\\s*mcmcfile\\s*=", lines, ignore.case = TRUE)
  if (length(mcmc_idx) == 1) {
    lines[mcmc_idx] <- paste("mcmcfile =", normalizePath(file.path("fit/output", paste0(chrom, ".mcmc.txt")), mustWork = FALSE))
  }
  outfile_idx <- grep("^\\s*outfile\\s*=", lines, ignore.case = TRUE)
  if (length(outfile_idx) == 1) {
    lines[outfile_idx] <- paste("outfile =", normalizePath(file.path("fit/output", paste0(chrom, ".out.txt")), mustWork = FALSE))
  }

  # Record provenance at the top of every generated control.
  lines <- c(
    paste0("# Generated species-level ajan vs alas IM control from: ", source_ctl),
    "# Interior and Seward samples are collapsed within each species; bidirectional migration is fitted.",
    lines
  )

  write_lines(lines, file.path("fit/controls", paste0(chrom, ".ctl")))
}

message("Built 9 empirical two-lineage IM controls in fit/controls/")
message("Collapsed sample counts: ajan=", counts[["ajan"]], ", alas=", counts[["alas"]])
