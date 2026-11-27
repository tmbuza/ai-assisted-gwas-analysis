# GWAS tools and installation

Copy these three files into the main GWAS project folder, replacing requirements.txt. The optional file refers to the main file in the same folder.

## Install the expanded working environment

With the project's virtual environment active:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
```

| Packages | Purpose | Scaffold status |
|---|---|---|
| NumPy, pandas, SciPy | Arrays, tables and statistical distributions | Used now |
| Matplotlib | Manhattan, Q–Q, intervals and other static plots | Used now |
| Seaborn | Statistical plot styling, distributions and heatmaps | Used for style and distribution plots |
| Plotly | Interactive HTML figures with hover/zoom | Dependency added; integration pending |
| adjustText | Place labels with fewer overlaps | Dependency added; integration pending |
| statsmodels | Statistical models, diagnostic and multiple-testing utilities | Dependency added; method-specific implementation pending |
| PyArrow | Columnar Parquet/Arrow input and output | Dependency added; streaming/import implementation pending |
| Requests | HTTP acquisition of source files | Dependency added; source-specific downloader pending |
| PyYAML | Parse source YAML study metadata | Dependency added; source-specific importer pending |
| tqdm | Progress bars for longer processing | Dependency added; integration pending |

Adding a dependency does not enable a method. The current workflow still reads its prepared TSV in memory, preserves source p-values, and performs the summary-statistics inspections described in the guide. The current static plotting code combines Matplotlib with Seaborn and exports PNG figures.

The version ranges are starting constraints, not an exact environment lock or a guarantee for every Python/OS combination. The existing numerical core was imported successfully in the preparation environment; the entire expanded and optional stacks have not been installed or resolved together. After successful project installation and smoke checks, record the environment actually used:

```bash
python -m pip freeze > requirements-lock.txt
```

Record the Python version, OS and environment purpose with the lockfile. A pip freeze from one platform is not automatically portable to another.

## Optional working tools

```bash
python -m pip install -r requirements-optional.txt
python -m pip check
```

This adds JupyterLab/ipykernel for exploration, openpyxl for spreadsheet supplements, pytest for development, pysam for indexed genomic data and Kaleido for Plotly static export. Install only when the corresponding work is planned; none is needed to reproduce the current synthetic analysis beyond the main numerical packages.

Kaleido v1 uses an external Chrome/Chromium installation. The official Plotly guide documents browser discovery and installation. Interactive Plotly HTML does not need Kaleido. Matplotlib already supports PNG, PDF and SVG; adding Seaborn does not require changing the image exporter.

## Specialist GWAS tools: separate scope and setup

| Tool | Potential role | What is required before adoption |
|---|---|---|
| GWASLab | Specialist summary-statistics QC, visualization and harmonization | Verify version compatibility and references; use a dedicated environment if dependencies conflict |
| PLINK 2 | Genotype QC and association/LD-related tasks | Separate platform-specific executable, appropriate genotype/reference inputs and a documented analysis plan |
| pysam / tabix access | Indexed genomic-file retrieval | Correctly compressed/indexed files and a verified coordinate/build convention |

GWASLab and PLINK are not installed by these requirements files. Their presence would not supply the matched LD, genotype data, phenotype definitions or allele harmonization needed for new analyses. LD clumping, fine mapping, polygenic scores and causal methods remain outside the current scaffold.

## Primary documentation

- Seaborn: https://seaborn.pydata.org/installing.html
- Plotly export and Chrome requirements: https://plotly.com/python/static-image-export/
- adjustText: https://adjusttext.readthedocs.io/en/stable/
- PyArrow/Parquet: https://arrow.apache.org/docs/python/parquet/
- pysam installation: https://pysam.readthedocs.io/en/stable/installation.html
- GWASLab: https://cloufield.github.io/gwaslab/
- PLINK 2: https://www.cog-genomics.org/plink/2.0/
