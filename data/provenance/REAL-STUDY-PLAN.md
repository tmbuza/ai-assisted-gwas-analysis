# Verified real-study source: GLGC HDL

The completed case study uses the GLGC 2013 joint GWAS/Metabochip HDL release. The source URL, checksum, original/prepared counts and exclusions are recorded in glgc-hdl/preparation-audit.json. See data/glgc-hdl/README.md and examples/real-study/ for the finished workflow and interpretation.

The original release has 2,447,441 rows. Preparation retains the 2,261,360 with valid autosomal coordinates and available reference frequency. Missing rsIDs are labelled by source coordinate and allele pair; the original source file remains the authority. No p-value threshold or random row subset is used to prepare the analysis input.

The guide preserves reference-frequency meaning, unavailable INFO, reported N decimals and source p-value text. Source methods are documented in the paper and authors' release page. Source data-use terms and citation requirements remain applicable.
