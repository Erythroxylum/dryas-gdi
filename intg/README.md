# Integrifolia gdi analysis

## Question

This analysis tests genealogical divergence between two geographically structured *Dryas integrifolia* populations:

- `intg_nGL_Nslope` — northern Greenland / North Slope population
- `intg_CAswGL` — Canadian Arctic / southwest Greenland population

Our working biological hypothesis is that these are strongly structured populations of the same species. The fitted population topology is paraphyletic with respect to `hook` (*D. hookeriana*), creating two possible biological interpretations: (1) *D. hookeriana* represents a derived/peripheral lineage that arose from within a structured *D. integrifolia*, or (2) gene flow involving `hook` and `intg_CAswGL` has affected the inferred population topology. The present gdi test focuses specifically on whether the two non-sister *integrifolia* populations show the level of genealogical divergence expected of distinct species.

## Fitted demographic model

Parameter source: `aih-prior3-s16-p4`.

We prefer `aih-prior3` over `aih-prior2` because the prior3 model gave a substantially more stable root-age posterior, whereas the prior2 root-age posterior was broad and close to its prior. The fitted guide-tree topology is:

```text
(ajan_alas,(intg_nGL_Nslope,(intg_CAswGL,hook)IH)IIH)R;
```

The A00 model included 18 directional migration parameters. In the prior3 analysis, the only migration edges supported across all chromosomes were the two directions between the *integrifolia* populations:

```text
intg_nGL_Nslope -> intg_CAswGL
intg_CAswGL     -> intg_nGL_Nslope
```

For gdi simulations we nevertheless retain the **full fitted MSC-M model**, including all modeled populations, internal nodes, and migration parameters, because all of these can affect the distribution of focal gene trees.

## Why three sequences?

Kornai et al. (2024), Eq. 13, recommends simulation of only three focal sequences when gene flow from other populations makes analytical calculation of gdi impractical. Yuttapong Thawornwattana confirmed that this is the appropriate design here: retain the full demographic model but simulate three sequences from the two lineages being tested, in two sets (`aab` and `abb`).

For the non-sister *integrifolia* populations, their most recent common ancestral population is node `IIH`, so `tau_IIH` is the divergence-time threshold in both directions.

### aab

```text
intg_nGL_Nslope, intg_nGL_Nslope, intg_CAswGL
```

This estimates:

```text
gdi_N = P(the two intg_nGL_Nslope sequences coalesce first and before tau_IIH)
```

### abb

```text
intg_nGL_Nslope, intg_CAswGL, intg_CAswGL
```

This estimates:

```text
gdi_CA = P(the two intg_CAswGL sequences coalesce first and before tau_IIH)
```

This treatment follows Kornai et al.'s explicit treatment of gdi between non-sister populations in a potentially paraphyletic species.

## Reproducible files

- `parameters/aih-prior3-s16-p4-means.csv` — chromosome-specific posterior means used to parameterize simulations.
- `imap/aab.imap.txt` and `imap/abb.imap.txt` — three-sequence sampling maps.
- `scripts/generate_controls.R` — generates the 18 chromosome/configuration-specific BPP controls.
- `scripts/run_gdi.R` — runs BPP simulations and scores Eq. 13 gdi.

Generate controls with:

```bash
Rscript intg/scripts/generate_controls.R
```

This creates:

```text
intg/controls/ch1_aab.ctl ... ch9_aab.ctl
intg/controls/ch1_abb.ctl ... ch9_abb.ctl
```

Each control retains the full four-population/18-migration-edge demographic model and simulates 1,000,000 gene trees of nominal locus length 50.

## Running the analysis

BPP 4.8.7 was used during development. Set the path to the executable, for example on Dawson's Mac:

```bash
export BPP_BIN=/Users/dawsonwhite/programs/bpp-4.8.7-macos-aarch64/bin/bpp
```

Then run:

```bash
Rscript intg/scripts/run_gdi.R
```

The script generates the controls if needed, runs all 18 simulations, calculates Eq. 13 gdi values, and writes:

```text
intg/output/gdi_intg.csv
intg/output/gdi_intg_long.csv
```

By default the large simulated gene-tree files are deleted after scoring. To retain them:

```bash
KEEP_TREES=true Rscript intg/scripts/run_gdi.R
```

## Additional test to consider

A separate sister-lineage comparison can evaluate `intg_CAswGL` versus `hook` using node `IH` as the divergence threshold. This addresses whether *D. hookeriana* is genealogically distinct from its immediate population-tree sister and is conceptually separate from the present non-sister *integrifolia* test.
