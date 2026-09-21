# data

Everything is read from inside the repository; paths are set in `src/swb/config.py`.

| Folder | Content | In git |
|---|---|---|
| `raw/snv/` | 480 per-run SNV VCFs, one folder per run: `TestCase <N>/<run>.vcf` | no |
| `raw/indel/` | 320 per-run indel VCFs | no |
| `truth/` | truth-set VCFs | yes |
| `reference/` | `TestCases.csv`, `TestCases_Indels.csv`, `sorted_exome_hc.bed.gz`, `stratification/*.bed.gz` | yes |
| `derived/` | small CSV caches and legacy tables used for parity; `*.parquet` caches are built by notebook 01 | CSVs only |

## Raw VCFs

`raw/` is gitignored. Copy the VCFs into `raw/snv/` (and `raw/indel/`) keeping one folder per run. Only notebooks 00 and 01 read `raw/snv/`.

- Names starting with `._` (macOS AppleDouble files) are ignored; `swb.io.real_files()` skips them.
- Only files ending in `.vcf` are read; `.vcf.gz` copies beside them are ignored.

## BED files

`sorted_exome_hc.bed.gz` and `stratification/*.bed.gz` are gzipped; `swb.io.read_bed` opens them with gzip.
