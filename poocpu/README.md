# PoOcPu1 gdi analyses

This folder contains the two remaining Eurasian/global-model gdi comparisons that Dawson and Yuttapong agreed should use chromosome-specific posterior means from the global H4d model and the PoOcPu1 demographic submodel:

1. `octo_EU` vs `octo_Carp_MK`
2. `Pu` vs `octo_RU_SJ` (use the exact H4d population label for RU/SJ when parameterizing the model)

Yuttapong agreed that the H4d posterior estimates can be used with the PoOcPu1 submodel for these gdi calculations. Each directional gdi is calculated from a three-sequence simulation (`aab` and `abb`) under the full PoOcPu1 MSC-M model, following Kornai et al. (2024), Eq. 13.

## Important: model specification still required

The scripts are scaffolded but intentionally refuse to generate controls until the exact PoOcPu1 species tree, internal-node labels, migration edges, and chromosome-specific H4d posterior means are supplied. We should not reconstruct these from memory or simplify the fitted model.

Place the following files here before running:

```text
poocpu/parameters/h4d-s47-p9-means.csv
poocpu/model/poocpu_model.R
```

The parameter CSV should contain one row per chromosome (`ch1`-`ch9`) and one numeric column for every `tau`, `theta`, and `W` parameter used by PoOcPu1. `poocpu_model.R` must define the exact model template and the divergence node used for each focal pair.

Once those inputs are added, the workflow is:

```bash
Rscript poocpu/scripts/generate_controls.R
BPP_BIN=/full/path/to/bpp Rscript poocpu/scripts/run_gdi.R
Rscript poocpu/scripts/summarize_results.R
```

The intended output is a one-row-per-chromosome table with four directional gdi values and a violin/point figure analogous to the `intg` analysis.
