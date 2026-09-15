# dryas-gdi

Reproducible workflow for estimating the genealogical divergence index (gdi) across chromosomes from gene trees simulated under fitted BPP multispecies-coalescent-with-migration (MSC-M) demographic models for *Dryas*.

## General approach

The gdi analyses are downstream of a broader BPP model-building workflow:

1. **Population-tree exploration.** Use broad population sampling, admixture analyses, phylogenetic inference, and sliding-window population trees to identify stable lineage relationships, topological heterogeneity, and candidate demographic structure.

2. **Pairwise BPP A00/IM screening.** Fit pairwise isolation-with-migration models to identify candidate gene-flow events. These analyses are used to inform more complex demographic models; they are not simply concatenated into a global migration network because omitted populations can create ghost-gene-flow artifacts.

3. **Fit complex BPP A00 MSC-M models.** Estimate divergence times (`tau`), population sizes (`theta`), and migration rates (`W`) in broader demographic context, chromosome by chromosome.

4. **Summarize posterior behavior across chromosomes.** Compare chromosome-specific posterior estimates of `tau`, `theta`, and `W` to assess concordance, asymmetry, and outlying chromosomes.

5. **Simulate gene trees for gdi.** Following Kornai et al. (2024), retain the fitted MSC-M demographic model and chromosome-specific posterior mean parameter values, but sample only three sequences from the two focal populations at a time:
   - `aab`: two sequences from population A and one from population B
   - `abb`: one sequence from population A and two from population B

6. **Calculate gdi from simulations.** Under the Kornai et al. (2024) definition (`gdi_K`, Eq. 13), gdi is the proportion of simulated gene trees in which the two sequences from the focal population coalesce first and before the relevant population-divergence time.

7. **Compare gdi across chromosomes and with other evidence.** gdi is treated as a measure of genealogical exclusivity under the fitted MSC-M model rather than as an automatic species-delimitation rule. Results are interpreted alongside morphology, geography, admixture, phylogenetic structure, and demographic history.

Reference:

Kornai D., Jiao X., Ji J., Flouri T., Yang Z. 2024. Hierarchical Heuristic Species Delimitation Under the Multispecies Coalescent Model with Migration. *Systematic Biology* 73:1015–1037.  
https://doi.org/10.1093/sysbio/syae050

## Implemented analyses

### `intg/`

Migration-aware gdi analyses within the AIH clade using the fitted `aih-prior3-s16-p4` demographic model.

Focal comparisons include:

- *D. integrifolia*: northern Greenland & Brooks Range population vs. southwestern Greenland & central Canada population 
- *D. hookeriana* vs. *D. integrifolia* southwestern Greenland & central Canada population 

The non-sister *integrifolia* comparison uses the age of the relevant most recent common ancestral population as the gdi cutoff.

Main scripts:

- `intg/scripts/generate_controls.R`
- `intg/scripts/run_gdi.R`

### `poocpu/`

Migration-aware gdi analyses within the Eurasian *D. octopetala–D. punctata* complex using the fitted 'H4D' demographic model.

Focal comparisons include:

- European vs. Carpathian *D. octopetala* populations
- *D. punctata* vs. the Russia / Svalbard lineage of *D. octopetala*

Main scripts:

- `poocpu/scripts/generate_controls.R`
- `poocpu/scripts/run_gdi.R`

### `ajan-alas/four-pop-validation/`

Validation of geographic population structure within *D. ajanensis* and *D. alaskensis* using a fitted four-population model 'ajan-alas-m3-prior2-s20-p4'.

Focal comparisons include:

- Interior vs. Seward *D. ajanensis*
- Interior vs. Seward *D. alaskensis*

In response to the geographic structure of the phylogenetic inference ((ajan_seward, alas_seward),(ajan_interior, alas_interior)), these analyses test whether the geographic populations behave as independently exclusive lineages before collapsing them for the direct species-level comparison.

Main script:

- `ajan-alas/four-pop-validation/model/run_gdi.R`

### `ajan-alas/species-im/`

Direct species-level *D. ajanensis* vs. *D. alaskensis* analysis.

Interior and Seward populations were first collapsed within species. A new two-species isolation-with-migration model,

`(ajan, alas)R`

was then fit independently to each of the nine nuclear chromosomes to estimate chromosome-specific `theta`, `tau_R`, and bidirectional migration rates. These newly fitted species-level posterior means were used for reciprocal `aab` and `abb` gene-tree simulations and gdi estimation.

Main scripts:

- `ajan-alas/species-im/scripts/build_fit_controls.R`
- `ajan-alas/species-im/scripts/run_fit.R`
- `ajan-alas/species-im/scripts/extract_fit_parameters.R`
- `ajan-alas/species-im/gdi/scripts/generate_controls.R`
- `ajan-alas/species-im/gdi/scripts/run_gdi.R`

## Current results

The implemented analyses show strong heterogeneity in genealogical exclusivity among recent *Dryas* lineages:

- *D. ajanensis* Interior vs. Seward: low gdi
- *D. alaskensis* Interior vs. Seward: low gdi
- direct species-level *D. ajanensis* vs. *D. alaskensis*: uniformly low gdi across chromosomes
- European vs. Carpathian / Kola *D. octopetala*: intermediate gdi
- *D. punctata* vs. Russia / Svalbard–Japan: asymmetric low-to-intermediate gdi
- AIH comparisons: generally intermediate gdi

These results reinforce the use of gdi as one component of a broader population- and species-boundary framework rather than as a standalone taxonomic classifier.
