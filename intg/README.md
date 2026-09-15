# Integrifolia gdi analysis

## Questions

This analysis evaluates genealogical divergence among three lineages in the `aih-prior3-s16-p4` demographic model:

- `intg_nGL_Nslope` — northern Greenland / North Slope *Dryas integrifolia*
- `intg_CAswGL` — Canadian Arctic / southwest Greenland *D. integrifolia*
- `hook` — *D. hookeriana*

The fitted population topology is:

```text
(ajan_alas,(intg_nGL_Nslope,(intg_CAswGL,hook)IH)IIH)R;
```

This creates two complementary gdi questions:

1. Are the two geographically structured, non-sister *integrifolia* populations genealogically divergent enough to behave as distinct species?
2. Is *D. hookeriana* genealogically distinct from its immediate population-tree sister, `intg_CAswGL`?

Our working biological hypothesis is that the two *integrifolia* lineages are strongly structured populations of the same species. The paraphyletic topology permits at least two biological interpretations: (1) *D. hookeriana* represents a derived/peripheral lineage that arose from within a structured *D. integrifolia*, or (2) gene flow has affected the inferred population topology.

## Fitted demographic model

Parameter source: `aih-prior3-s16-p4`.

We prefer `aih-prior3` over `aih-prior2` because the prior3 model gave a substantially more stable root-age posterior, whereas the prior2 root-age posterior was broad and close to its prior.

The A00 model included 18 directional migration parameters. In the prior3 analysis, the only migration edges supported across all chromosomes were the two directions between the *integrifolia* populations:

```text
intg_nGL_Nslope -> intg_CAswGL
intg_CAswGL     -> intg_nGL_Nslope
```

For gdi simulations we retain the **full fitted MSC-M model**, including all modeled populations, internal nodes, and migration parameters, because all of these can affect the distribution of focal gene trees.

## Three-sequence simulations

Following Kornai et al. (2024), Eq. 13, gdi is estimated by retaining the full demographic model while simulating only three focal sequences at a time. Yuttapong Thawornwattana confirmed this design: use two complementary configurations, `aab` and `abb`, for each pair being tested.

### `intg_nGL_Nslope` vs `intg_CAswGL`

These are non-sister populations whose MRCA is node `IIH`, so `tau_IIH` is used as the divergence-time threshold in both directions.

```text
aab: intg_nGL_Nslope, intg_nGL_Nslope, intg_CAswGL
abb: intg_nGL_Nslope, intg_CAswGL, intg_CAswGL
```

These estimate:

```text
gdi_N  = P(two intg_nGL_Nslope sequences coalesce first and before tau_IIH)
gdi_CA = P(two intg_CAswGL sequences coalesce first and before tau_IIH)
```

This follows Kornai et al.'s treatment of gdi between non-sister populations in a potentially paraphyletic species.

### `intg_CAswGL` vs `hook`

These are sisters at node `IH`, so `tau_IH` is the divergence-time threshold.

```text
aab: intg_CAswGL, intg_CAswGL, hook
abb: intg_CAswGL, hook, hook
```

These estimate directional gdi for `intg_CAswGL` and `hook` relative to one another.

## Reproducible files

- `parameters/aih-prior3-s16-p4-means.csv` — chromosome-specific posterior means used to parameterize simulations.
- `imap/` — three-sequence sampling maps for both comparisons.
- `scripts/generate_controls.R` — generates 36 chromosome/configuration-specific BPP controls.
- `scripts/run_gdi.R` — runs BPP simulations and scores Eq. 13 gdi.
- `scripts/summarize_results.R` — reformats results to one row per chromosome, writes summary statistics, and produces a violin/point plot across chromosomes.

Generate controls with:

```bash
Rscript intg/scripts/generate_controls.R
```

Each control retains the full four-population/18-migration-edge demographic model and simulates 1,000,000 gene trees of nominal locus length 50.

## Running the analysis

BPP 4.8.7 was used during development. Set the executable path, for example:

```bash
export BPP_BIN=/Users/dawsonwhite/programs/bpp-4.8.7-macos-aarch64/bin/bpp
```

Then run:

```bash
Rscript intg/scripts/run_gdi.R
```

The script calculates all four directional gdi estimates for each chromosome and writes:

```text
intg/output/gdi_intg_long.csv
intg/output/gdi_intg.csv
```

By default large simulated gene-tree files are deleted after scoring. To retain them:

```bash
KEEP_TREES=true Rscript intg/scripts/run_gdi.R
```

## Summarizing existing results

The simulations do **not** need to be rerun to fix the wide output table or make figures. Starting from `gdi_intg_long.csv`, run:

```bash
Rscript intg/scripts/summarize_results.R
```

This writes:

```text
intg/output/gdi_intg.csv
intg/output/gdi_intg_summary.csv
intg/output/figures/gdi_intg_violin.pdf
intg/output/figures/gdi_intg_violin.png
```

`gdi_intg.csv` has one row per chromosome and four directional gdi columns. The violin plot summarizes the nine chromosome-specific estimates for each directional test, overlays the individual chromosome estimates and their mean, and marks the conventional 0.2 and 0.7 gdi thresholds.
