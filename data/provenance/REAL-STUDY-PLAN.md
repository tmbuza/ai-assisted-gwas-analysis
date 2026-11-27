# Real-study acceptance plan

Status: not selected. This is the next content-development milestone.

Start with a published quantitative trait and full genome-wide summary statistics from a primary consortium or the GWAS Catalog. Do not use a Catalog significant-associations export as if it were full summary statistics.

Before selection, verify:

- Original publication and exact study/trait/model accession.
- Accessible source file and accompanying metadata; permitted redistribution and attribution.
- Complete background associations, not a p-value-selected hit list.
- Genome build, effect allele, other allele, beta and SE scale, source test, frequency and N definitions.
- Source imputation-quality field and its definition; if absent, design an explicit compatible extension rather than fabricate a score.
- Population descriptors, relatedness/ancestry adjustment, study exclusions and sample overlap.
- A practical import strategy for file size, compression and missing data.
- An importer that keeps source values, hashes and a row-count/exclusion audit.
- A named independently checked variant and verified figure interpretations.

Core full-file calculations can be retained for large sources while using a documented plotting strategy. Do not silently downsample, remove the null background, or use a selected subset for genome-wide calibration. This starter currently reads its table in memory; a multi-million-row source may need chunked import and plot aggregation before the real lesson is completed.

Regional LD colouring needs an appropriate population/build-matched LD reference and harmonized variants. Fine mapping, credible sets and polygenic scores need separately justified inputs and methods. Do not infer them from position alone.

Authoritative starting points:
- https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/
- https://www.ebi.ac.uk/gwas/downloads/summary-statistics
- https://github.com/gwas-catalog/summary-statistics-standard

The lesson contract is intentionally simpler than GWAS-SSF. Standard-compliant files may contain unsupported missing values, indels or effect types; a documented importer is required.
