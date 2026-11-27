# Synthetic GWAS v2 validation

- All 24 guide pages render with Quarto 1.11.9.
- Reference dataset: 22,000 input variants, 16,344 after MAF, 15,760 after MAF and INFO.
- Ten figure families generated as PNG files.
- Every retained normal Wald p-value independently reproduced with standard-library erfc; derived association strengths and both interval endpoints checked.
- Simulated SE approximation, allele orientation/sign transformation and fixed reference exception checked against maintainer truth files.
- Shared histogram bins sum to the full input and retained counts; sequential filtering percentages verified.
- Lambda GC and threshold counts independently checked; input and script hashes match.
- SYN000001 independently verifies with unchanged beta .12, SE .02, z=6 and Wald p≈1.97317529008e-9.
- Invalid zero p, negative SE, duplicate ID, out-of-range INFO and identical alleles rejected.
- Selected-variant inputs suppress Q–Q/lambda and remove stale Q–Q exports including historical PDF/SVG copies.
- Synthetic report template renders with six figure paths; learner prompts intentionally remain in the reusable template.
- Figures visually reviewed, including chromosome labels, allele/position labels, shared bins, central Q–Q view and numerical-rounding scale in p-value agreement.

Detailed numerical record: data/provenance/validation.json. Local link/image audit: data/provenance/render-validation.json.

Only the actual numerical plotting environment was exercised; the entire expanded/optional requirements stack was not installed together. This remains a teaching simulation. No real GWAS dataset, estimated LD, independent-locus or causal claim is included.
