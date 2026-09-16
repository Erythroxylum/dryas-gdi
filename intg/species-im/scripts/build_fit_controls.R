# Build empirical BPP A00 controls for a direct two-population D. integrifolia
# isolation-with-migration (IM) model.
#
# Biological question
# -------------------
# Seven of nine chromosome-level CASTER trees support a monophyletic
# D. integrifolia composed of two strongly structured geographic populations,
# whereas the broader AIH BPP MSC-M model places intg_CAswGL with hookeriana.
# This analysis therefore fits the alternative species hypothesis directly:
#
#   (intg_nGL_Nslope, intg_CAswGL)R;
#
# with bidirectional migration between the two integrifolia populations.
# The fitted chromosome-specific theta, tau_R, and W values can subsequently be
# used for reciprocal aab/abb gdi simulations.
#
# IMPORTANT: the original AIH empirical chromosome files contain additional
# populations. This script subsets each multilocus PHYLIP file to samples mapped
# to intg_nGL_Nslope or intg_CAswGL before constructing the two-population model.
#
# Usage from the repository root:
#   Rscript intg/species-im/scripts/build_fit_controls.R \
#     /path/to/aih-s16/multilocus \
#     /path/to/aih-s16-p4.imap.txt \
#     /path/to/bpp-a00-aih-prior3-s16-p4.ctl

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# Resolve user-supplied paths before changing working directories.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Usage: Rscript intg/species-im/scripts/build_fit_controls.R ",
    "/path/to/multilocus /path/to/imap.txt /path/to/template.ctl"
  )
}
data_dir <- normalizePath(args[1], mustWork = TRUE)
source_imap <- normalizePath(args[2], mustWork = TRUE)
template_ctl <- normalizePath(args[3], mustWork = TRUE)

# Locate intg/species-im from this script's own path so generated files are
# portable and do not depend on the directory from which Rscript is called.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

chromosomes <- paste0("ch", 1:9)
source_seqfiles <- file.path(data_dir, paste0(chromosomes, ".txt"))
missing_seqfiles <- source_seqfiles[!file.exists(source_seqfiles)]
if (length(missing_seqfiles)) {
  stop("Missing chromosome data files: ", paste(missing_seqfiles, collapse = ", "))
}
source_seqfiles <- normalizePath(source_seqfiles)
names(source_seqfiles) <- chromosomes

# Read the original AIH sample map and retain only the two geographic populations
# currently assigned to D. integrifolia.
imap <- read_table(
  source_imap,
  col_names = c("sample", "population"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

target_populations <- c("intg_nGL_Nslope", "intg_CAswGL")
target_imap <- imap |>
  filter(population %in% target_populations)

if (!all(target_populations %in% target_imap$population)) {
  stop(
    "Source imap must contain both populations: ",
    paste(target_populations, collapse = ", ")
  )
}
if (anyDuplicated(target_imap$sample)) {
  stop("Duplicate sample IDs found in the target integrifolia imap")
}

target_samples <- target_imap$sample
counts <- table(factor(target_imap$population, levels = target_populations))
if (any(counts == 0)) stop("Both integrifolia populations must contain samples")

# Write the filtered two-population imap used by every chromosome.
dir.create("fit/imap", recursive = TRUE, showWarnings = FALSE)
new_imap <- file.path("fit", "imap", "integrifolia_two_pop.imap.txt")
write.table(
  target_imap,
  file = new_imap,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = " "
)
new_imap <- normalizePath(new_imap, mustWork = TRUE)

# Subset a BPP multilocus PHYLIP file to the target samples. The chromosome files
# used in this project are sequential multilocus PHYLIP: each locus starts with a
# line containing "nseq nsites", followed by nseq one-line sequence records.
# Empty lines between loci are tolerated. Loci lacking one or more target samples
# are discarded so every retained locus has the same expected individual set.
subset_multilocus <- function(infile, outfile, keep_samples) {
  x <- readLines(infile, warn = FALSE)
  header_pattern <- "^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]*$"
  out <- character(0)
  i <- 1L
  n_loci_seen <- 0L
  n_loci_kept <- 0L

  while (i <= length(x)) {
    if (!grepl(header_pattern, x[i])) {
      i <- i + 1L
      next
    }

    header <- str_split(str_trim(x[i]), "[[:space:]]+", simplify = TRUE)
    nseq <- as.integer(header[1])
    nsites <- as.integer(header[2])
    n_loci_seen <- n_loci_seen + 1L

    # Collect the next nseq non-empty sequence records.
    records <- character(0)
    j <- i + 1L
    while (j <= length(x) && length(records) < nseq) {
      if (nzchar(str_trim(x[j]))) records <- c(records, x[j])
      j <- j + 1L
    }
    if (length(records) != nseq) {
      stop("Incomplete locus after line ", i, " in ", infile)
    }

    sample_ids <- vapply(
      str_split(str_trim(records), "[[:space:]]+"),
      function(z) z[1],
      character(1)
    )
    keep <- sample_ids %in% keep_samples
    kept_ids <- sample_ids[keep]

    # Retain only loci containing every requested target sample exactly once.
    if (length(kept_ids) == length(keep_samples) &&
        setequal(kept_ids, keep_samples) &&
        !anyDuplicated(kept_ids)) {
      out <- c(
        out,
        paste(length(keep_samples), nsites),
        records[keep],
        ""
      )
      n_loci_kept <- n_loci_kept + 1L
    }

    i <- j
  }

  if (n_loci_kept < 1L) {
    stop("No complete target loci retained from ", infile)
  }

  writeLines(out, outfile)
  c(seen = n_loci_seen, kept = n_loci_kept)
}

# Produce chromosome-specific two-population empirical sequence files. Keeping
# these under fit/data makes the model completely explicit and avoids relying on
# BPP to ignore sequences from populations absent from the new imap.
dir.create("fit/data", recursive = TRUE, showWarnings = FALSE)
subset_stats <- vector("list", length(chromosomes))
names(subset_stats) <- chromosomes
seqfiles <- character(length(chromosomes))
names(seqfiles) <- chromosomes

for (chrom in chromosomes) {
  outfile <- file.path("fit", "data", paste0(chrom, ".txt"))
  stats <- subset_multilocus(source_seqfiles[[chrom]], outfile, target_samples)
  subset_stats[[chrom]] <- stats
  seqfiles[[chrom]] <- normalizePath(outfile, mustWork = TRUE)
}

nloci_by_chrom <- vapply(subset_stats, function(z) unname(z[["kept"]]), integer(1))
message(
  "Retained complete integrifolia loci: ",
  paste(paste0(chromosomes, "=", nloci_by_chrom), collapse = ", ")
)

# Read the fitted AIH prior3 control as the template. Priors, MCMC settings, and
# other empirical-analysis options are retained unless they depend on the number
# of populations. In particular, the W prior is inherited from the prior3 file
# rather than hard-coded here.
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

# Replace the original four-population AIH species block with the direct
# two-population integrifolia species hypothesis.
replace_species_block <- function(lines) {
  i <- grep("^\\s*species&tree\\s*=", lines, ignore.case = TRUE)
  if (length(i) != 1) stop("Expected exactly one species&tree block")
  if (i + 2 > length(lines)) stop("Incomplete species&tree block")

  replacement <- c(
    "species&tree = 2 intg_nGL_Nslope intg_CAswGL",
    paste("                ", paste(as.integer(counts), collapse = " ")),
    "(intg_nGL_Nslope, intg_CAswGL)R;"
  )
  append(lines[-c(i:(i + 2))], replacement, after = i - 1)
}

# Replace the original AIH migration graph with the two reciprocal migration
# directions relevant to the direct integrifolia population comparison.
replace_migration_block <- function(lines) {
  i <- grep("^\\s*migration\\s*=", lines, ignore.case = TRUE)
  if (length(i) != 1) stop("Expected exactly one migration block")
  n_old <- suppressWarnings(as.integer(str_extract(lines[i], "[0-9]+\\s*$")))
  if (is.na(n_old)) stop("Could not parse existing migration count")
  end <- i + n_old
  if (end > length(lines)) stop("Existing migration block is incomplete")

  replacement <- c(
    "migration = 2",
    "            intg_nGL_Nslope intg_CAswGL",
    "            intg_CAswGL intg_nGL_Nslope"
  )
  append(lines[-c(i:end)], replacement, after = i - 1)
}

# Create one empirical A00 control per chromosome.
dir.create("fit/controls", recursive = TRUE, showWarnings = FALSE)
dir.create("fit/output", recursive = TRUE, showWarnings = FALSE)

for (chrom in chromosomes) {
  lines <- template_lines

  lines <- replace_directive(lines, "seqfile", seqfiles[[chrom]])
  lines <- replace_directive(lines, "Imapfile", new_imap)
  lines <- replace_directive(lines, "nloci", as.character(nloci_by_chrom[[chrom]]))

  # Give every chromosome a unique output prefix.
  job_prefix <- normalizePath(file.path("fit", "output", chrom), mustWork = FALSE)
  lines <- replace_directive(lines, "jobname", job_prefix)

  # Update model-dependent directives for the two-population analysis.
  lines <- replace_species_block(lines)
  lines <- replace_migration_block(lines)
  lines <- replace_directive(lines, "phase", "1 1")

  # Record provenance and the hypothesis being tested at the top of each control.
  lines <- c(
    paste0("# Direct two-population D. integrifolia IM fit for ", chrom),
    paste0("# Original empirical sequence file: ", source_seqfiles[[chrom]]),
    paste0("# Filtered empirical sequence file: ", seqfiles[[chrom]]),
    paste0("# Complete target loci retained: ", nloci_by_chrom[[chrom]]),
    paste0("# Source sample map: ", source_imap),
    paste0("# Source BPP template: ", template_ctl),
    "# Hypothesis: monophyletic integrifolia with two strongly structured populations.",
    "# Migration is fitted in both directions between the two populations.",
    "# Priors, including wprior, are inherited from the aih-prior3 template.",
    lines
  )

  write_lines(lines, file.path("fit", "controls", paste0(chrom, ".ctl")))
}

message("Built 9 direct integrifolia two-population IM controls in fit/controls/")
message(
  "Sample counts: intg_nGL_Nslope=", counts[["intg_nGL_Nslope"]],
  ", intg_CAswGL=", counts[["intg_CAswGL"]]
)
message("Priors inherited from the supplied aih-prior3 BPP template")
