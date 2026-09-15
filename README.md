# dryas-gdi

Reproducible workflow for estimating the genealogical divergence index (gdi) across chromosomes from gene trees simulated under fitted BPP multispecies-coalescent-with-migration (MSC-M) demographic models for *Dryas*.

## General approach

The gdi analyses are downstream of a broader BPP model-building workflow:

1. **Population-tree exploration.** Use the larger population sampling to identify stable population relationships and candidate demographic structure.
2. **Pairwise BPP A00/IM screening.** Fit pairwise isolation-with-migration models to identify candidate gene-flow events. These pairwise analyses are used as evidence for building more complex demographic models; they are not simply concatenated into a global migration network because omitted populations can create ghost-gene-flow artifacts.
3. **Fit complex BPP A00 MSC-M models.** Estimate divergence times (`tau`), population sizes (`theta`), and migration rates (`W`) in the broader demographic context, chromosome by chromosome.
4. **Summarize posterior behavior across chromosomes.** Make figures for `tau`, `theta`, and `W` posterior estimates to assess chromosome-to-chromosome concordance and identify outliers.
5. **Simulate gene trees for gdi.** Following Kornai et al. (2024), use the full fitted MSC-M model and chromosome-specific posterior mean parameter values, but sample only three sequences from the two focal populations at a time: `aab` (two from A, one from B) and `abb` (one from A, two from B).
6. **Calculate gdi from simulations.** Under the Kornai et al. (2024) definition (`gdi_K`, Eq. 13), gdi is the proportion of simulated gene trees in which the two sequences from the focal population coalesce first and before the relevant population-divergence time.
7. **Compare gdi across chromosomes and other evidence.** gdi is treated as one line of evidence rather than an automatic species-delimitation rule and is interpreted alongside geography, morphology, admixture, and the fitted demographic history.

Reference: Kornai D., Jiao X., Ji J., Flouri T., Yang Z. 2024. Hierarchical Heuristic Species Delimitation Under the Multispecies Coalescent Model with Migration. *Systematic Biology* 73:1015–1037. https://doi.org/10.1093/sysbio/syae050

## Analyses

- [`intg/`](intg/) — non-sister *D. integrifolia* population comparison under the `aih-prior3-s16-p4` model.

Additional focal comparisons will be added as their simulation designs are finalized.

## Current analysis tasks

- [ ] Run and summarize the `intg` gdi simulations across chromosomes.
- [ ] Build chromosome-level posterior summary figures for `tau`, `theta`, and `W`.
- [ ] Add the remaining focal gdi comparisons after validating the `intg` workflow.
