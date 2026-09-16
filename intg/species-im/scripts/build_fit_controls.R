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
# populations. This script subsets each multilocus PHYLIP file to sequences whose
# BPP labels end in ^<sample>, where <sample> occurs in one of the two focal
# integrifolia populations in the source imap.
#
# Usage from the repository root:
#   Rscript intg/species-im/scripts/build_fit_controls.R \
#     /path/to/aih-s16/multilocus \
#     /path/to/aih-s16-p4.imap.txt \
#     /path/to/bpp-a00-aih-s16-p4-short.ctl

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

# Return the biological sample ID from a BPP sequence label. Empirical labels in
# these files have the form locus_identifier^sample_ID, for example:
#   0_OY992832.1_8020_8705^ajan_BABY6724
# The imap contains only the suffix (ajan_BABY6724), so matching the complete
# sequence label directly to the imap would incorrectly retain zero sequences.
extract_sample_id <- function(sequence_label) {
  ifelse(
    grepl("\\^", sequence_label),
    sub("^.*\\^", "", sequence_label),
    sequence_label
  )
}

# Subset a BPP multilocus PHYLIP file to the focal populations. Each locus starts
# with "nseq nsites", followed by nseq one-line sequence records. Unlike the
# previous implementation, a locus does NOT need to contain all eight target
# individuals: missing samples vary among loci in these GBS data. We retain every
# locus containing at least one sequence from EACH focal population and rewrite
# its nseq header to the number of retained sequences.
subset_multilocus <- function(infile, outfile, target_map) {
  x <- readLines(infile, warn = FALSE)
  header_pattern <- "^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]*$"
  out <- character(0)
  i <- 1L
  n_loci_seen <- 0L
  n_loci_kept <- 0L
  retained_sequences <- 0L

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

    sequence_labels <- vapply(
      str_split(str_trim(records), "[[:space:]]+"),
      function(z) z[1],
      character(1)
    )
    sample_ids <- extract_sample_id(sequence_labels)
    keep <- sample_ids %in% target_map$sample

    kept_records <- records[keep]
    kept_samples <- sample_ids[keep]
    kept_pops <- target_map$population[match(kept_samples, target_map$sample)]

    # A two-population IM locus is informative/valid only if both populations are
    # represented. We do not require every individual to be present at every locus.
    if (length(kept_records) >= 2L && all(target_populations %in% kept_pops)) {
      out <- c(
        out,
        paste(length(kept_records), nsites),
        kept_records,
        ""
      )
      n_loci_kept <- n_loci_kept + 1L
      retained_sequences <- retained_sequences + length(kept_records)
    }

    i <- j
  }

  if (n_loci_kept < 1L) {
    stop("No loci containing both focal populations retained from ", infile)
  }

  writeLines(out, outfile)
  c(
    seen = n_loci_seen,
    kept = n_loci_kept,
    sequences = retained_sequences
  )
}

# Produce chromosome-specific two-population empirical sequence files. Keeping
# these under fit/data makes the model completely explicit and prevents nonfocal
# AIH sequences from entering the direct integrifolia fit.
dir.create("fit/data", recursive = TRUE, showWarnings = FALSE)
subset_stats <- vector("list", length(chromosomes))
names(subset_stats) <- chromosomes
seqfiles <- character(length(chromosomes))
names(seqfiles) <- chromosomes

for (chrom in chromosomes) {
  outfile <- file.path("fit", "data", paste0(chrom, ".txt"))
  stats <- subset_multilocus(source_seqfiles[[chrom]], outfile, target_imap)
  subset_stats[[chrom]] <- stats
  seqfiles[[chrom]] <- normalizePath(outfile, mustWork = TRUE)
}

nloci_by_chrom <- vapply(subset_stats, function(z) unname(z[["kept"]]), integer(1))
message(
  "Retained integrifolia loci containing both populations: ",
  paste(paste0(chromosomes, "=", nloci_by_chrom), collapse = ", ")
)

# Read the supplied AIH control as the template. Priors, MCMC settings, and other
# empirical-analysis options are retained unless they depend on population count.
# In particular, wprior is inherited from the template rather than hard-coded.
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

# Replace the original AIH migration graph with reciprocal migration between the
# two geographic integrifolia populations.
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
    paste0("# Loci containing both focal populations: ", nloci_by_chrom[[chrom]]),
    paste0("# Source sample map: ", source_imap),
    paste0("# Source BPP template: ", template_ctl),
    "# Hypothesis: monophyletic integrifolia with two strongly structured populations.",
    "# Migration is fitted in both directions between the two populations.",
    "# Priors, including wprior, are inherited from the supplied AIH template.",
    lines
  )

  write_lines(lines, file.path("fit", "controls", paste0(chrom, ".ctl")))
}

message("Built 9 direct integrifolia two-population IM controls in fit/controls/")
message(
  "Sample counts: intg_nGL_Nslope=", counts[["intg_nGL_Nslope"]],
  ", intg_CAswGL=", counts[["intg_CAswGL"]]
)
message("Priors inherited from the supplied AIH BPP template")
