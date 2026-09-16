# Direct two-population D. integrifolia IM analysis

## Question

Seven of nine chromosome-level CASTER trees support a monophyletic *D. integrifolia* composed of two strongly structured geographic populations:

- `intg_nGL_Nslope` — northern Greenland / North Slope
- `intg_CAswGL` — Canadian Arctic / southwest Greenland

The broader `aih-prior3-s16-p4` BPP MSC-M model instead inferred a paraphyletic population history in which `intg_CAswGL` is sister to *D. hookeriana*.

This analysis directly fits the species hypothesis favored by the majority of chromosome-level CASTER trees:

```text
(intg_nGL_Nslope, intg_CAswGL)R;
```

with continuous migration in both directions between the two geographic populations.

The purpose is not to replace the broader AIH demographic model. It asks a complementary question: **given a monophyletic D. integrifolia species hypothesis, how much genealogical differentiation has accumulated between its two major geographic populations after accounting for migration?**

## Workflow

1. Start from the empirical chromosome-specific AIH multilocus sequence files and AIH sample map.
2. Retain only samples assigned to `intg_nGL_Nslope` or `intg_CAswGL`.
3. Retain loci containing the complete target sample set and write explicit two-population chromosome files under `fit/data/`.
4. Use the fitted `aih-prior3-s16-p4` BPP control as the template for priors and MCMC settings.
5. Fit a new two-population A00 isolation-with-migration model independently for each of the nine nuclear chromosomes.
6. Estimate chromosome-specific `theta`, `tau_R`, and reciprocal `W` values.
7. Use those newly fitted parameters in reciprocal `aab`/`abb` simulations for a direct within-*integrifolia* gdi analysis.

## Scripts

- `scripts/build_fit_controls.R` — subsets the empirical AIH chromosome files, writes the two-population imap, counts retained loci, and generates nine BPP controls.
- `scripts/run_fit.R` — runs chromosome-specific BPP fits in parallel; set `BPP_JOBS` to control concurrency.
- `scripts/extract_fit_parameters.R` — extracts posterior mean demographic parameters after the fits finish.

## Build controls

From the repository root:

```bash
Rscript intg/species-im/scripts/build_fit_controls.R \
  /path/to/aih-s16/multilocus \
  /path/to/aih-s16-p4.imap.txt \
  /path/to/bpp-a00-aih-prior3-s16-p4.ctl
```

Before starting BPP, inspect a generated control and confirm the population order, sample counts, priors, and chromosome-specific `nloci`:

```bash
cat intg/species-im/fit/controls/ch1.ctl
```

## Run fits

Example for a nine-way parallel run:

```bash
BPP_BIN=/full/path/to/bpp BPP_JOBS=9 \
  Rscript intg/species-im/scripts/run_fit.R
```

BPP is invoked with `--no-pin` so the operating system can schedule the independent chromosome processes.

## Extract posterior means

After all nine fits complete:

```bash
Rscript intg/species-im/scripts/extract_fit_parameters.R
```

This writes:

```text
intg/species-im/fit/species_im_means.csv
```

The subsequent gdi workflow should use these directly fitted chromosome-specific parameters rather than parameters from the paraphyletic four-population AIH model.
