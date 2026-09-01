# Ajanensis vs alaskensis gdi

This folder is for the species-level gdi comparison between *Dryas ajanensis* and *D. alaskensis*.

The available fitted model `ajan-alas-m3-prior2-s20-p4` contains four sampled populations:

```text
((ajan_Interior, ajan_Seward), (alas_Interior, alas_Seward))
```

Yuttapong's existing gdi analysis tested only Interior vs Seward within each species and found low gdi (~0.1). He explicitly noted that this did **not** test ajanensis vs alaskensis. He later preferred retaining the original four-population model because the two-population/ghost analysis had little information for ajan-vs-alas gene flow.

## Remaining design decision

For a species-level three-sequence Eq. 13 simulation we must specify how `ajanensis` and `alaskensis` are represented while retaining the fitted four-population demographic model. We should not silently choose Interior or Seward as the species representative, nor merge fitted populations without defining how the demographic parameters transform.

The scripts here therefore use the same reusable gdi machinery but intentionally require a completed `model/ajan_alas_species_model.R` before generating controls.

Required inputs:

```text
ajan-alas/parameters/ajan-alas-m3-prior2-s20-p4-means.csv
ajan-alas/model/ajan_alas_species_model.R
```

Once the species-level sampling/model representation is agreed, run:

```bash
Rscript ajan-alas/scripts/generate_controls.R
BPP_BIN=/full/path/to/bpp Rscript ajan-alas/scripts/run_gdi.R
Rscript ajan-alas/scripts/summarize_results.R
```
