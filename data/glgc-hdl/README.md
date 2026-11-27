# Prepared GLGC 2013 HDL input

summary-statistics.tsv.gz contains 2,261,360 autosomal source SNVs with available 1000 Genomes EUR A1 reference frequency. It preserves source beta/SE, positive p-value text and reported N values. study-metadata.json declares GRCh37, unavailable INFO and the reference-frequency basis.

Original: https://csg.sph.umich.edu/willer/public/lipids2013/ ; Willer et al. (2013), doi:10.1038/ng.2797.

Reproduce with python scripts/python/prepare_glgc.py from the project root. Source checksum and exclusions are in data/provenance/glgc-hdl/. No significance-based preparation sampling is used. Do not treat reference frequency as study frequency or derived Wald p as the corrected source meta-analysis test.
