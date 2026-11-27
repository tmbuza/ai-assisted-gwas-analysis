# AI-Assisted GWAS Analysis with Python

CDI scaffold for checking, plotting, verifying and explaining quantitative-trait GWAS summary statistics.

**Status:** runnable synthetic workflow and a rendered 24-page guide; real-study and own-study routes are scaffolded. No real GWAS dataset is bundled yet. Expanded lessons and verified real-data content follow in the next phase.

## Learning route

Synthetic → real study → your own study. Each route follows define, preview, meet inputs, check, analyse, verify/explain and report. Begin with `00-setup.qmd` or `docs/index.html`.

## Quick start

From the project root, with a Python 3.10+ environment active:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
python scripts/python/run_workflow.py --data-dir data/synthetic --output-dir results/synthetic --check-only
python scripts/python/run_workflow.py --data-dir data/synthetic --output-dir results/synthetic
python scripts/python/verify_variant.py data/synthetic results/synthetic SYN000001
quarto render
```

Quarto is installed separately. Rendering displays results and does not execute the analysis. To regenerate fictional inputs, run `python scripts/python/prepare_synthetic.py`.

## Bundled synthetic example

22,000 fictional variant rows across 22 autosomes; 15,760 retained under the included configuration. Per-variant N is simulated, not a count of variant rows. IDs and coordinates are not biological annotations.

Ten figure families: Manhattan, Q–Q, selected effect forest, frequency/INFO, effect-frequency, sample-size distribution, filtering audit and position-only regional zoom, SE-versus-N precision and source-versus-Wald p-value agreement. Plot data, QC flags and run records are saved alongside figures.

SYN000001: beta 0.12 per A allele relative to G, SE 0.02, z=6, p≈1.97317529e-9, 95% CI≈[0.08080072,0.15919928].

## Scope and limits

The starter accepts autosomal biallelic SNV quantitative-trait beta summaries with all documented fields. It is not a general importer or a genotype GWAS pipeline. It does not implement LD clumping, fine mapping, genetic risk prediction, binary-trait odds ratios or cross-study harmonization. Missing measurements must not be fabricated to satisfy the schema.

Full genome-wide files support global diagnostics; hits-only files do not. `summary_scope=selected_variants` suppresses Q–Q and lambda. Source p-values are retained even when they differ from derived normal Wald p-values.

The implementation loads the table in memory. Large real-study sources may need a verified streaming importer and a documented plotting strategy. See `data/provenance/REAL-STUDY-PLAN.md` before choosing a dataset.

## Folder guide

| Folder | Contents |
|---|---|
| `examples/` | Seven-step lessons for each route |
| `data/` | Bundled synthetic inputs and real/own-study instructions |
| `config/` | Visible filtering and display parameters |
| `scripts/python/` | Generation, validation, analysis and independent verification |
| `results/synthetic/` | Ten figures and numerical evidence |
| `templates/` | Reusable report and metadata starters |
| `library/` | Reference and provenance reading list |
| `assets/css/` | Shared CDI stylesheet |
| `docs/` | Rendered guide |

`reports/`, virtual environments and caches are ignored. No hosting domain is assumed. Configure repository/download URLs and domain settings when the new repository is created.

## Synthetic version 2

The generator now uses varied invented chromosome spans/densities, a mixture
of low-frequency and common variants, coverage concentrated near nominal
N=20000 with a lower tail, and INFO/SE patterns that depend on frequency and
coverage. These are documented teaching choices, not a fitted human population.
Effects and EAF are reversed consistently when the named effect allele is the
major allele. The fixed reference variant remains unchanged.

Seaborn is now used for style and distributions. All ten figure families are exported as PNG. Plotly and specialist-tool integrations remain pending.

See data/provenance/SIMULATION-DESIGN.md and VALIDATION.md. The full expanded
requirements stack has not been resolved together; the workflow's actual
NumPy/pandas/SciPy/Matplotlib/Seaborn environment is recorded in run-record.json.
