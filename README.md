# dryas-gdi

Reproducible chromosome-wise workflow for estimating migration-aware genealogical divergence index (gdi) from BPP multispecies-coalescent-with-migration (MSC-M) models in *Dryas*.

## Overview

The repository contains the demographic-model fits, parameter extraction, gene-tree simulations, and gdi summaries used to compare genealogical exclusivity among selected *Dryas* lineages. The workflow is downstream of broader population-genomic analyses and is intended to quantify lineage exclusivity under explicit migration models, not to treat gdi as a standalone species classifier.

The general workflow is:

1. Use admixture analyses, phylogenetic inference, sliding-window population trees, and geography to define population structure and formulate candidate demographic models.
2. Use pairwise BPP A00/IM analyses to screen candidate migration edges.
3. Fit larger chromosome-specific BPP A00 MSC-M models and estimate `theta`, `tau`, and migration rate `W`.
4. Extract chromosome-specific posterior means.
5. Simulate 1,000,000 three-sequence gene trees per test under the fitted model:
   - `aab`: two sequences from lineage A and one from B, estimating gdi for A.
   - `abb`: one sequence from A and two from B, estimating gdi for B.
6. Score gdi as the proportion of simulated trees in which the two focal sequences coalesce first and before the relevant divergence-time cutoff.
7. Compare chromosome-level gdi values across demographic models and with morphology, geography, admixture, and phylogenetic structure.

For sister pairs the cutoff is their divergence time. For non-sister pairs, the cutoff is the age of the ancestral population subtending the two focal lineages.

Reference: Kornai D., Jiao X., Ji J., Flouri T., Yang Z. 2024. Hierarchical Heuristic Species Delimitation Under the Multispecies Coalescent Model with Migration. *Systematic Biology* 73:1015–1037. https://doi.org/10.1093/sysbio/syae050

## Repository structure

```text
dryas-gdi/
├── intg/
│   ├── parameters/              # AIH chromosome-specific posterior means
│   ├── controls/                # full-AIH gdi simulation controls
│   ├── scripts/                 # control generation, simulation, summaries
│   ├── output/                  # full-AIH gdi results
│   └── species-im/              # direct two-population D. integrifolia IM analysis
│       ├── fit/                 # original W-prior-mean-2 empirical fits
│       ├── sensitivity/         # migration-prior sensitivity fits
│       ├── scripts/             # fit construction, fitting, parameter extraction
│       └── gdi/                 # prior2 vs prior20 gdi simulations
├── poocpu/
│   ├── model/                   # reduced H4d-derived demographic model
│   ├── parameters/              # chromosome-specific posterior means
│   ├── controls/                # gdi simulation controls
│   ├── scripts/                 # parameter extraction, simulation, summaries
│   └── output/
└── ajan-alas/
    ├── four-pop-validation/     # within-species Interior/Seward tests
    │   ├── model/
    │   ├── parameters/
    │   ├── controls/
    │   ├── scripts/
    │   └── output/
    └── species-im/              # direct D. ajanensis vs D. alaskensis IM analysis
        ├── fit/
        ├── scripts/
        └── gdi/
```

The analysis-specific directories preserve the fitted demographic model, extracted parameters, simulation controls, and resulting gdi tables so that each comparison can be reproduced independently.

## Analysis-specific details

### AIH / *D. integrifolia* and *D. hookeriana* — `intg/`

The full AIH analysis uses the fitted `aih-prior3-s16-p4` MSC-M model and retains all migration edges during gene-tree simulation.

Two comparisons are evaluated:

- northern Greenland/North Slope *D. integrifolia* vs. Canadian Arctic/southwestern Greenland *D. integrifolia*;
- Canadian Arctic/southwestern Greenland *D. integrifolia* vs. *D. hookeriana*.

In the fitted AIH topology, the two *D. integrifolia* populations are non-sister lineages. Their gdi test therefore uses the age of node `IIH`, the ancestral population subtending both focal populations, rather than a direct sister-pair divergence time.

Main workflow:

```text
intg/parameters/aih-prior3-s16-p4-means.csv
intg/scripts/generate_controls.R
intg/scripts/run_gdi.R
intg/scripts/summarize_results.R
```

### Direct two-population *D. integrifolia* IM model — `intg/species-im/`

Chromosome-level CASTER analyses frequently recover the two major *D. integrifolia* geographic populations together, whereas the broader AIH BPP model places the Canadian Arctic/southwestern Greenland population with *D. hookeriana*. To evaluate the monophyletic *D. integrifolia* hypothesis directly, a new two-population IM model was fit independently to each nuclear chromosome:

```text
(intg_nGL_Nslope, intg_CAswGL)R;
```

The empirical sequence files were reduced to the eight *D. integrifolia* samples in the original AIH imap. A locus was retained when both focal populations were represented; individual samples were allowed to be missing locus by locus.

The original direct fit used `W ~ Gamma(2,1)` (prior mean 2). Migration-rate posteriors were strongly prior-sensitive, so the same nine chromosome fits were repeated with:

```text
W ~ Gamma(2,0.1)   # prior mean 20
W ~ Gamma(2,0.01)  # prior mean 200
```

The mean-20 fits changed inferred `W` substantially but left `theta` and `tau` broadly stable. The mean-200 fits produced strong parameter coupling, very low ESS for `tau`, and unstable estimates on several chromosomes, indicating poor identifiability under the highly permissive prior. Final gdi sensitivity comparisons therefore use the two well-behaved parameterizations: prior mean 2 and prior mean 20.

Despite the tenfold shift in the migration prior, gdi was nearly unchanged:

```text
                         prior mean 2    prior mean 20
nGL/North Slope              0.286           0.294
CA/southwestern Greenland    0.366           0.374
```

This analysis shows that the direct *D. integrifolia* gdi result is robust even though the absolute migration-rate estimates are prior-sensitive.

Main workflow:

```text
intg/species-im/scripts/build_fit_controls.R
intg/species-im/scripts/run_fit.R
intg/species-im/scripts/extract_fit_parameters.R
intg/species-im/scripts/extract_sensitivity_parameters.R
intg/species-im/gdi/scripts/generate_controls.R
intg/species-im/gdi/scripts/run_gdi.R
```

### Eurasian *D. octopetala–D. punctata* comparisons — `poocpu/`

These analyses use a reduced H4d-derived MSC-M model:

```text
(Po,((octo_EU,octo_Carp_MK)Octo,(Pu,RU_SJ)Punc)AOcPu)R;
```

Focal tests are:

- European vs. Carpathian/Macedonian *D. octopetala*;
- *D. punctata* vs. the Russia/Svalbard-Japan lineage.

The reduced model retains three migration edges from H4d: `octo_EU -> RU_SJ`, `RU_SJ -> octo_EU`, and `Po -> Punc`. The focal pairs themselves therefore do not necessarily have direct migration edges, but gdi is simulated under the complete reduced migration graph.

Main workflow:

```text
poocpu/model/poocpu_model.R
poocpu/parameters/h4d-s47-p9-means.csv
poocpu/scripts/generate_controls.R
poocpu/scripts/run_gdi.R
poocpu/scripts/summarize_results.R
```

### *D. ajanensis–D. alaskensis* four-population validation — `ajan-alas/four-pop-validation/`

The four-population model is:

```text
((ajan_Interior,ajan_Seward)J,(alas_Interior,alas_Seward)L)R;
```

The fitted model retains four migration edges between geographically corresponding *ajanensis* and *alaskensis* populations. The gdi tests themselves compare Interior vs. Seward populations within each species:

- *D. ajanensis* Interior vs. Seward;
- *D. alaskensis* Interior vs. Seward.

These tests validate the behavior of strongly structured geographic populations before the species-level populations are collapsed.

Main workflow:

```text
ajan-alas/four-pop-validation/model/four_pop_model.R
ajan-alas/four-pop-validation/parameters/four_pop_means.csv
ajan-alas/four-pop-validation/scripts/generate_controls.R
ajan-alas/four-pop-validation/scripts/run_gdi.R
```

### Direct *D. ajanensis* vs. *D. alaskensis* IM model — `ajan-alas/species-im/`

Interior and Seward populations were collapsed within species and a direct two-species IM model was fit independently to each chromosome:

```text
(ajan,alas)R;
```

The fitted chromosome-specific `theta`, `tau_R`, and bidirectional `W` values were then used for reciprocal `aab`/`abb` simulations.

Main workflow:

```text
ajan-alas/species-im/scripts/build_fit_controls.R
ajan-alas/species-im/scripts/run_fit.R
ajan-alas/species-im/scripts/extract_fit_parameters.R
ajan-alas/species-im/gdi/scripts/generate_controls.R
ajan-alas/species-im/gdi/scripts/run_gdi.R
```

## Interpretation

Across these analyses, gdi is used as a quantitative description of genealogical exclusivity under explicitly fitted MSC-M histories. Low, intermediate, and high gdi values are interpreted together with morphology, geography, admixture, phylogenetic relationships, and ecological differentiation rather than as automatic taxonomic decisions.

Current results include low within-species Interior/Seward gdi in *D. ajanensis* and *D. alaskensis*, low direct species-level gdi between those two morphologically differentiated taxa, intermediate and prior-robust gdi between the two main *D. integrifolia* populations, and intermediate/asymmetric gdi among the Eurasian focal lineages.
