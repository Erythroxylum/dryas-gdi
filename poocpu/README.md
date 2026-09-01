# PoOcPu1 gdi analyses

This folder contains two gdi comparisons parameterized from chromosome-specific posterior means of the global H4d model and evaluated under the reduced PoOcPu1 MSC-M submodel:

1. `octo_EU` vs `octo_Carp_MK`
2. `Pu` vs `RU_SJ`

Each directional gdi is calculated from a three-sequence simulation (`aab` and `abb`) under the full PoOcPu1 model, following Kornai et al. (2024), Eq. 13.

## PoOcPu1 demographic model

The reduced population tree is:

```text
(Po,((octo_EU,octo_Carp_MK)Octo,(Pu,RU_SJ)Punc)AOcPu)R;
```

The three migration edges retained in this submodel are:

```text
octo_EU -> RU_SJ
RU_SJ   -> octo_EU
Po       -> Punc
```

The two gdi cutoffs are therefore:

```text
octo_EU vs octo_Carp_MK : tau_Octo
Pu vs RU_SJ             : tau_Punc
```

## Mapping PoOcPu1 onto H4d

The reduced-model labels are mapped onto H4d as follows. This mapping is important because the root of PoOcPu1 is not the root of the global H4d tree.

```text
Po       = H4d octo_Pobeda (tip 7)
Pu       = H4d punc (tip 9)
RU_SJ    = H4d octo_RU_SJ (tip 5)
Octo     = H4d node 17
Punc     = H4d node 19
AOcPu    = H4d node 16
R        = H4d node S (node 13)
```

Thus `tau_R` and `theta_R` in the PoOcPu1 controls come from H4d `tau:13:S` and `theta:13:S`, not from H4d's global node `R`.

The migration mappings are:

```text
W_EU_to_RUSJ = H4d W:4->5:octo_EU->octo_RU_SJ
W_RUSJ_to_EU = H4d W:5->4:octo_RU_SJ->octo_EU
W_Po_to_Punc = H4d W:7->19:octo_Pobeda->Punc
```

## Workflow

First pull the chromosome-specific posterior means out of Yuttapong's formatted H4d summary:

```bash
Rscript poocpu/scripts/extract_parameters.R /path/to/h4d-s47-p9-fmt.csv
```

This writes:

```text
poocpu/parameters/h4d-s47-p9-means.csv
```

Then generate the 36 controls (9 chromosomes x 2 comparisons x 2 directional configurations):

```bash
Rscript poocpu/scripts/generate_controls.R
```

Inspect at least one control from each comparison before launching the full run. Then:

```bash
BPP_BIN=/full/path/to/bpp Rscript poocpu/scripts/run_gdi.R
```

Finally:

```bash
Rscript poocpu/scripts/summarize_results.R
```

The analysis writes a long results table, a one-row-per-chromosome summary with four directional gdi values, and a violin/point figure analogous to the `intg` analysis. Large simulated gene-tree files are deleted after scoring unless `KEEP_TREES=true` is set.
