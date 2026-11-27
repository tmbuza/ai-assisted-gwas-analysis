# Synthetic GWAS version 2: design and limits

This is a teaching simulation, not a fitted human population, genotype model, or empirical benchmark. Seed: 27112026. The bundled TSV and hashes define the reference release; exact regeneration across package versions is not guaranteed.

- 22,000 rows across 22 invented spans decreasing roughly from 215 to 42 Mb, with jitter. These are NOT reference-assembly chromosome lengths.
- Background marker counts depend on span and random density variation. Three regions receive 120 additional markers each.
- MAF is a mixture: 45% lower-frequency component and 55% common component. The component distributions and filtering consequences are intentional teaching choices.
- Most nominal per-variant N values are near 20,000; 18% are drawn from a lower-coverage component. A small frequency-dependent coverage reduction is applied.
- INFO depends on frequency with additional noise and a small poorer-quality component. This dependence is illustrative, not a universal relationship.
- SE uses 1/sqrt(2*N*MAF*(1-MAF)*INFO), except the fixed reference row.
- Latent standard-normal noise has distance-dependent adjacent correlation exp(-distance/50 kb) within chromosomes. Three smooth effect profiles use 120 kb width. This does not constitute estimated LD or causal-variant truth.
- Effects are generated relative to a simulated minor allele. When the named effect allele is major, beta is negated and EAF=1-MAF. This preserves z magnitude, p and MAF while respecting allele orientation.
- SYN000001 is overridden to A/G, beta .12, SE .02, EAF .25, INFO .98 and N 20000. Its preserved reference arithmetic is useful but it is an exception to the population simulation formula.
- Source p-values use the two-sided normal Wald formula. Source/derived agreement is therefore expected and is not an independent discovery.

Files synthetic-layout.tsv, synthetic-regions.tsv and synthetic-truth.tsv expose design details for maintainers. The analysis never reads these truth files to select peaks, filter variants, identify independent loci or assign genes.

The version is more varied than the original toy distribution while remaining explicitly fictional. Neither the exclusion rate nor signal strength is a claim about a typical human GWAS.
