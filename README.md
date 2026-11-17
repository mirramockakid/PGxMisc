# PGxMisc
general helper scripts

## aldy_to_report.R

Summarises a directory of [Aldy](https://github.com/0xTCG/aldy) `.out` files into star-allele reports.
Requires R with the `readr`, `dplyr` and `tidyr` packages.

```bash
Rscript aldy_to_report.R --path ./results
```

`--path` is a directory; every `*.out` file in it is read. Four CSVs are written to the current working directory:

| File | Contents |
| --- | --- |
| `star_alleles.csv` | One row per sample, gene and Aldy solution (Sample, Gene, SolutionID, Major, Minor) |
| `ambiguous.csv` | Rows from samples/genes with more than one solution |
| `star_allele_freq.csv` | Count and frequency of each major star allele per gene |
| `suballele_freq.csv` | Count and frequency of each sub-allele per gene |

Frequencies count every solution of an ambiguous sample, so the denominator is the number of alleles counted per gene rather than 2 × samples.
