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
#
# The data directory must contain ch1.txt through ch9.txt. The source control is
# used as a template for the empirical-data and MCMC settings. We explicitly set
# W ~ Gamma(2, 0.01), matching the final m3-prior2 four-population analysis that
# supplied the validation parameters used elsewhere in this repository.

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
})

# Resolve all user-supplied paths before changing the working directory. This
# allows both absolute and repository-relative paths to work predictably.
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

# Locate this species-IM analysis directory from the script path. Generated
# controls, sample maps, and outputs are kept inside this workflow directory.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
analysis_dir <- dirname(dirname(script_path))
setwd(analysis_dir)

# Verify that all nine empirical chromosome data files are present before we
# generate any controls.
chromosomes <- paste0("ch", 1:9)
seqfiles <- file.path(data_dir, paste0(chromosomes, ".txt"))
missing_seqfiles <- seqfiles[!file.exists(seqfiles)]
if (length(missing_seqfiles)) {
  stop("Missing chromosome data files: ", paste(missing_seqfiles, collapse = ", "))
}
seqfiles <- normalizePath(seqfiles)
names(seqfiles) <- chromosomes

# Read the original four-population sample map. The s20 dataset contains 20
# individuals assigned to ajan_Interior, ajan_Seward, alas_Interior, or
# alas_Seward.
imap <- read_table2(
  source_imap,
  col_names = c("sample", "population"),
  col_types = cols(.default = col_character())
)

# Collapse geography within each taxon so the empirical fit estimates the two
# species-level lineages requested for gdi: ajan and alas.
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

# Derive species-level sample counts from the map rather than hard-coding them.
# For the current s20 data these should be 9 ajan and 11 alas individuals.
counts <- table(factor(collapsed$population, levels = c("ajan", "alas")))
if (any(counts == 0)) stop("Both ajan and alas must have sampled individuals")

# Write one collapsed sample map used by all nine chromosome fits.
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

# Read the empirical BPP control that supplied the original data/MCMC settings.
# We preserve settings such as phase, nloci, cleandata, theta/tau priors,
# finetune, burnin, sampfreq, and nsample unless explicitly replaced below.
template_lines <- read_lines(template_ctl)

# Helper: replace exactly one simple key = value directive in the template.
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

# Helper: replace the original three-line four-population species block with the
# two-lineage species model and the sample counts derived from the collapsed imap.
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

# Helper: replace the original multi-edge migration graph with the species-level
# bidirectional IM graph. These are the two migration rates required for the
# subsequent migration-aware gdi calculation.
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

# Create one empirical A00 control per chromosome. Each fit uses the same 20
# individuals and collapsed species map but a different chromosome-specific
# multilocus sequence file.
dir.create("fit/controls", recursive = TRUE, showWarnings = FALSE)
dir.create("fit/output", recursive = TRUE, showWarnings = FALSE)

for (chrom in chromosomes) {
  lines <- template_lines

  # Point BPP directly to the actual chromosome-specific empirical data and to
  # the collapsed ajan/alas sample map.
  lines <- replace_directive(lines, "seqfile", seqfiles[[chrom]])
  lines <- replace_directive(lines, "Imapfile", new_imap)

  # Give every chromosome a unique output prefix inside this workflow so fits do
  # not overwrite one another. The original template uses jobname rather than
  # separate outfile/mcmcfile directives.
  job_prefix <- normalizePath(file.path("fit/output", chrom), mustWork = FALSE)
  lines <- replace_directive(lines, "jobname", job_prefix)

  # Replace only the biological population model: four geographic populations
  # become two species, with bidirectional species-level migration.
  lines <- replace_species_block(lines)
  lines <- replace_migration_block(lines)

  # Match the final m3-prior2 demographic analysis used for our four-population
  # validation. A later prior-sensitivity run can repeat this fit with G(2,0.1).
  lines <- replace_directive(lines, "wprior", "2 0.01")

  # Add provenance and an explicit statement of the new fitted hypothesis.
  lines <- c(
    paste0("# Species-level ajan vs alas IM fit for ", chrom),
    paste0("# Empirical sequence file: ", seqfiles[[chrom]]),
    paste0("# Source sample map: ", source_imap),
    paste0("# Source BPP template: ", template_ctl),
    "# Interior and Seward are collapsed within species; migration is fitted in both directions.",
    "# Migration prior: W ~ Gamma(2, 0.01), matching m3-prior2.",
    lines
  )

  write_lines(lines, file.path("fit/controls", paste0(chrom, ".ctl")))
}

message("Built 9 empirical species-level ajan vs alas IM controls in fit/controls/")
message("Collapsed sample counts: ajan=", counts[["ajan"]], ", alas=", counts[["alas"]])
message("Migration prior: W ~ Gamma(2, 0.01) [m3-prior2]")
