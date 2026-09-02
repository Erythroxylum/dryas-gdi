# Species-level *ajanensis* vs *alaskensis* IM model

This workflow fits the two-lineage IM model that Yuttapong explicitly suggested for testing *Dryas ajanensis* versus *D. alaskensis* with gdi. His previous gdi calculation tested only Interior vs Seward within each species, not the species pair itself.

## Biological model

Interior and Seward samples are collapsed within each taxon:

```text
ajan_Interior + ajan_Seward -> ajan
alas_Interior + alas_Seward -> alas
```

The empirical BPP A00 model is:

```text
(ajan, alas)R;
```

with bidirectional migration:

```text
ajan -> alas
alas -> ajan
```

The scripts preserve the original chromosome sequence files, MCMC settings, mutation model, priors, and other estimation settings from the fitted four-population controls. Only the population assignment, species tree, and migration graph are changed.

Prefer using the four-population **prior3 / W ~ G(2,0.1)** empirical controls as templates because Yuttapong reported that this prior gave more stable root-age estimates than prior2 while recovering the same main demographic scenario.

## 1. Build the empirical controls

Provide the directory containing one final four-population empirical control per chromosome and the corresponding four-population imap:

```bash
Rscript ajan-alas/species-im/scripts/build_fit_controls.R \
  /path/to/four_pop_control_directory \
  /path/to/four_pop_imap.txt
```

The script collapses the original imap automatically, derives the number of ajan and alas samples, makes sequence-file paths absolute, and writes:

```text
ajan-alas/species-im/fit/controls/ch1.ctl ... ch9.ctl
ajan-alas/species-im/fit/imap/ajan_alas.imap.txt
```

Before running, inspect at least `ch1.ctl` to confirm that the source control, sequence file, priors, and MCMC settings are the intended ones.

## 2. Fit the chromosome-specific IM models

```bash
BPP_BIN=/full/path/to/bpp \
Rscript ajan-alas/species-im/scripts/run_fit.R
```

The fits are run sequentially and console output is written to:

```text
ajan-alas/species-im/fit/logs/
```

BPP posterior/MCMC outputs are written to `fit/output/` when the source controls contain `mcmcfile` and `outfile` directives.

## 3. Species-level gdi

After the fits complete, chromosome-specific posterior means will be extracted and used to simulate two three-sequence configurations:

```text
ajan, ajan, alas -> gdi_ajan
ajan, alas, alas -> gdi_alas
```

The divergence cutoff in both directions is `tau_R` from the fitted two-lineage model. The gdi simulation/parser scripts will be finalized after inspecting one completed BPP output so that posterior columns are mapped exactly rather than guessed.
