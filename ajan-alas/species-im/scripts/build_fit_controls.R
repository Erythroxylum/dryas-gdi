# Build empirical BPP A00 controls for the species-level ajanensis vs alaskensis
# isolation-with-migration (IM) model.
#
# This analysis fits a NEW demographic model to the empirical chromosome data:
#
#   (ajan, alas)R;
#
# Interior and Seward samples are collapsed within each species, and BPP estimates
# species-level theta, tau_R, and bidirectional migration directly from the data.
# We do not average parameters from the fitted four-population model.
#
# Usage from the repository root:
#   Rscript ajan-alas/species-im/scripts/build_fit_controls.R \
#     /path/to/ajan-alas-s20/multilocus \
#     /path/to/ajan-alas-s20-p4.imap.txt \
#     /path/to/bpp-a00-ajan-alas-s20-p4.ctl

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# Resolve user-supplied paths before changing working directories.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Usage: Rscript ajan-alas/species-im/scripts/build_fit_controls.R ",
    "/path/to/multilocus /path/to/imap.txt /path/to/template.ctl"
  )
}
data_dir <- normalizePath(args[1], mustWork = TRUE)
source_imap <- normalizePath(args[2], mustWork = TRUE)
template_ctl <- normalizePath(args[3], mustWork = TRUE)

# Locate the species-IM analysis directory from this script's location.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Verify that all nine chromosome sequence files are present.
chromosomes <- paste0("ch", 1:9)
seqfiles <- file.path(data_dir, paste0(chromosomes, ".txt"))
missing_seqfiles <- seqfiles[!file.exists(seqfiles)]
if (length(missing_seqfiles)) {
  stop("Missing chromosome data files: ", paste(missing_seqfiles, collapse = ", "))
}
seqfiles <- normalizePath(seqfiles)
names(seqfiles) <- chromosomes

# Count loci directly from each BPP multilocus PHYLIP file. Each locus begins
# with a header line containing exactly two integers: number of sequences and
# alignment length. Chromosomes contain different numbers of retained loci, so
# nloci must be set separately rather than copied as 2000 from the template.
count_loci <- function(path) {
  x <- readLines(path, warn = FALSE)
  n <- sum(grepl("^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]*$", x))
  if (n < 1L) stop("Could not identify any locus headers in ", path)
  n
}
nloci_by_chrom <- vapply(seqfiles, count_loci, integer(1))
message(
  "Detected loci per chromosome: ",
  paste(paste0(names(nloci_by_chrom), "=", nloci_by_chrom), collapse = ", ")
)

# Read the original four-population sample map.
imap <- read_table(
  source_imap,
  col_names = c("sample", "population"),
  col_types = cols(.default = col_character())
)

# Collapse Interior and Seward populations within each taxon.
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

# Derive species-level sample counts automatically.
counts <- table(factor(collapsed$population, levels = c("ajan", "alas")))
if (any(counts == 0)) stop("Both ajan and alas must have sampled individuals")

# Write one collapsed imap used by every chromosome.
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

# Read the empirical BPP control used as a template for priors and MCMC settings.
template_lines <- read_lines(template_ctl)

# Replace exactly one simple key=value directive.
replace_directive <- function(lines, key, value, required = TRUE) {
  idx <- grep(paste0("^\\s*", key, "\\s*="), lines, ignore.case = TRUE)
  if (!length(idx)) {
    if (required) stop("Template control is missing directive: ", key)
    return(lines)
  }
  if (length(idx) != 1) stop("Expected exactly one ", key, " directive")
  lines[idx] <- paste(key, "=", value)
  lines
}

# Replace the four-population species block with the two-lineage model.
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

# Replace the original migration graph with bidirectional species-level migration.
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

# Create one empirical A00 control per chromosome.
dir.create("fit/controls", recursive = TRUE, showWarnings = FALSE)
dir.create("fit/output", recursive = TRUE, showWarnings = FALSE)

for (chrom in chromosomes) {
  lines <- template_lines

  # Use the chromosome-specific sequence file, collapsed imap, and actual locus count.
  lines <- replace_directive(lines, "seqfile", seqfiles[[chrom]])
  lines <- replace_directive(lines, "Imapfile", new_imap)
  lines <- replace_directive(lines, "nloci", as.character(nloci_by_chrom[[chrom]]))

  # Give every chromosome a unique BPP output prefix.
  job_prefix <- normalizePath(file.path("fit/output", chrom), mustWork = FALSE)
  lines <- replace_directive(lines, "jobname", job_prefix)

  # Replace model-dependent directives for the collapsed two-species analysis.
  lines <- replace_species_block(lines)
  lines <- replace_migration_block(lines)
  lines <- replace_directive(lines, "phase", "1 1")
  lines <- replace_directive(lines, "wprior", "2 0.01")

  # Record provenance and key model choices at the top of each generated control.
  lines <- c(
    paste0("# Species-level ajan vs alas IM fit for ", chrom),
    paste0("# Empirical sequence file: ", seqfiles[[chrom]]),
    paste0("# Empirical loci detected: ", nloci_by_chrom[[chrom]]),
    paste0("# Source sample map: ", source_imap),
    paste0("# Source BPP template: ", template_ctl),
    "# Interior and Seward are collapsed within species; migration is fitted in both directions.",
    "# Phase = 1 1 because the collapsed model has two unphased species.",
    "# Migration prior: W ~ Gamma(2, 0.01), matching m3-prior2.",
    lines
  )

  write_lines(lines, file.path("fit/controls", paste0(chrom, ".ctl")))
}

message("Built 9 empirical species-level ajan vs alas IM controls in fit/controls/")
message("Collapsed sample counts: ajan=", counts[["ajan"]], ", alas=", counts[["alas"]])
message("Migration prior: W ~ Gamma(2, 0.01) [m3-prior2]")
