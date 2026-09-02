# Four-population ajanensis/alaskensis validation

This analysis reproduces the biological questions in Yuttapong's original ajan/alaskensis gdi calculation while using the same three-sequence Eq. 13 implementation as the other workflows in this repository.

The fitted model is:

```text
((ajan_Interior, ajan_Seward)J,
 (alas_Interior, alas_Seward)L)R;
```

with four migration edges:

```text
alas_Seward   -> ajan_Seward
ajan_Seward   -> alas_Seward
alas_Interior -> ajan_Interior
ajan_Interior -> alas_Interior
```

Yuttapong previously simulated 2+2+2+2 sequences and reported gdi values around 0.1 for Interior vs Seward within both species. Here, we test the same two merges separately with aab/abb three-sequence simulations:

```text
ajan_Interior, ajan_Interior, ajan_Seward -> gdi_ajan_Interior
ajan_Interior, ajan_Seward, ajan_Seward   -> gdi_ajan_Seward

alas_Interior, alas_Interior, alas_Seward -> gdi_alas_Interior
alas_Interior, alas_Seward, alas_Seward   -> gdi_alas_Seward
```

The cutoffs are `tau_J` for ajanensis and `tau_L` for alaskensis.

## Parameter source

Prefer the chromosome-specific posterior means from the four-population **W ~ G(2,0.1) / prior3** model because Yuttapong later considered its root-age estimates more stable than prior2. The scripts can also reproduce the prior2 result by supplying the prior2 posterior-mean table.

If you already have Yuttapong's `est-all_chr-mean.csv`, extract the required parameters with:

```bash
Rscript ajan-alas/four-pop-validation/scripts/extract_parameters.R \
  /path/to/est-all_chr-mean.csv
```

This writes:

```text
ajan-alas/four-pop-validation/parameters/four_pop_means.csv
```

## Generate and run simulations

```bash
Rscript ajan-alas/four-pop-validation/scripts/generate_controls.R
```

Then inspect a chromosome-1 ajan control and alaskensis control before running the batch. After validation:

```bash
BPP_BIN=/full/path/to/bpp \
Rscript ajan-alas/four-pop-validation/scripts/run_gdi.R
```

Results are written to:

```text
ajan-alas/four-pop-validation/output/gdi_four_pop_long.csv
ajan-alas/four-pop-validation/output/gdi_four_pop.csv
```

The expected validation target is the qualitative result Yuttapong reported: low Interior-vs-Seward gdi (approximately 0.1) within both species. Exact values may differ because the present implementation uses three focal sequences and scores the focal-pair terminal coalescence time directly.
