# data

Everything the notebooks read is inside this repository. Paths are set in `src/swb/config.py`; none is absolute. Sizes below are decimal (1 MB = 10^6 bytes). They come from `results/00_data_preparation/tables/inventory_summary.csv`, except where marked "measured" (file sizes read from disk).

## 1. What is committed and what is not

| Folder | Content | Size | In git |
|---|---|---:|---|
| `truth/` | truth-set VCFs: `hc_bed_filtered.recode.vcf` (1,161 SNVs), `hc_bed_filtered_INDELS.recode.vcf` (48 indels) | 1.8 MB | yes |
| `reference/` | `TestCases.csv` (864 rows), `TestCases_Indels.csv` (320 rows), `sequencing_yield.csv`, `sorted_exome_hc.bed.gz`, `stratification/*.bed.gz`, `validated_vcfs/`, `giab_regions/`, `Census_all.csv` | 857 MB (73 files) | yes |
| `derived/` | legacy tables and caches, used as parity targets (`sets_dict.csv`, `vcfcomparison_df_full.csv`, `union_metrics.csv`, ...), plus two files of values printed by the legacy notebooks: `legacy_indel_iou_recorded.csv` (pairwise IoU, cells 38, 39 and 43 of the indel notebook; parity target of `02_pairwise_iou_indels`) and `legacy_metrics_recorded.csv` (F1 summaries, factor tests and variance table, cells 126 to 134; parity target of `05_metrics_anova`) | 92.6 MB (35 files) | yes |
| `derived/sets_dict.parquet` | the variant-set cache built by notebook 01 (columns `key`, `variant`) | 1.4 MB | **no**, gitignored |
| `raw/snv/` | the 480 SNV VCFs, one folder per run (section 4) | 3.94 GB (960 files: 480 `.vcf` + 480 `.vcf.gz`) | **no**, gitignored |
| `raw/indel/` | the 320 indel VCFs, one folder per run | 179 MB (400 files: 320 `.vcf` + 80 `.vcf.gz`) | **no**, gitignored |
| `raw/vcf_snps/`, `raw/vcf_indels/` | older copies of the raw VCFs, see below | 3.41 GB and 172 MB (measured) | **yes**, tracked since commit `e89c3ab` |

Parts of `reference/` that no notebook reads (sizes measured): `giab_regions/` (uncompressed copies of the six `stratification/` BEDs, 257 MB) and the two `Census_all*.csv` files (identical to each other). `validated_vcfs/` (520 MB) is read by notebook 06, except its folder `NV` (9 files), which is not one of the study's five centers. The six `giab_regions/*.bed` files are stored with **Git LFS**; without `git-lfs` a clone holds small pointer files in their place (`git lfs pull` fetches the real files). Nothing else in the repository uses LFS.

**The raw VCFs are only partly in the repository.** `raw/vcf_snps/` holds 476 of the 480 SNV VCFs and `raw/vcf_indels/` holds all 320 indel VCFs (plus 10 `.vcf.gz`), under the old folder names. They are byte-identical (SHA-256) to the files in `raw/snv/` and `raw/indel/` on the author's disk, and they use the same `TestCase <N>/` folder layout. No notebook reads them. The four runs missing from `raw/vcf_snps/` are TestCase **26, 32, 38 and 44** (all EA, Altay + COSAP, BOWTIE, Mutect, duplicates MARK; the four `isTrimmed` x `baseRecalibration` combinations). Together the two old folders add about 3.6 GB to a checkout.

## 2. Notebooks 02 onward run from the committed cache alone

Only notebook 01 needs `raw/snv/`. It parses the 480 VCFs once and writes the cache `derived/sets_dict.parquet`. Every later notebook reads the cache and never re-parses a run VCF, unless it needs raw VCF fields and says so in its Purpose section. Notebook 00 lists and hashes `raw/` if it is present and runs without it.

Current status: the cache is gitignored, so a fresh clone does not contain `sets_dict.parquet`, and notebooks 02 and 07 (the notebooks after 01 that read the cache) stop with "run 01_variant_sets first" until notebook 01 has run. `02_pairwise_iou_indels` needs neither raw nor the parquet: it reads the committed `derived/sets_dict_indels.csv`, the legacy cache of the indel sets (no notebook parses the indel VCFs yet). The committed `derived/sets_dict.csv` holds the same sets in the legacy format (`Key`, `Identifier`); notebook 01's parity cell asserts they are equal for all 480 runs.

## 3. Obtaining the raw VCFs (to run notebook 01)

The raw VCFs are not on SRA. SRA holds the sequencing reads; the VCFs are the outputs of the 480 pipeline runs on those reads. To get them, use the copies in `raw/vcf_snps/` for 476 runs and regenerate the other four, or regenerate all 480.

### 3.1 Reads

SRA study **SRP162370** (SEQC2 somatic whole-exome sequencing, HCC1395 tumour and HCC1395BL normal cell lines). The run accessions are listed in Table 1 of the manuscript. The ten below are the accessions that appear in the run VCF file names; check them against Table 1. Each sample code (the `Sample` column of `TestCases.csv`) is one tumour/normal pair. Tumour and normal roles are from NCBI's SRA run information (checked 2026-09-21). In the run VCF file names the **normal comes first**, so a name never says which sample is the tumour: the analysis reads `##tumor_sample=` from the VCF header (`docs/legacy-deviations.md` D5).

| Sample code | Tumour (HCC1395) | Normal (HCC1395BL) |
|---|---|---|
| EA | SRR7890918 | SRR7890919 |
| FD | SRR7890879 | SRR7890880 |
| IL | SRR7890883 | SRR7890874 |
| LL | SRR7890850 | SRR7890851 |
| NC | SRR7890844 | SRR7890845 |

Download with the SRA Toolkit, for example `prefetch SRR7890919 && fasterq-dump --split-files SRR7890919`.

`reference/sequencing_yield.csv` holds the total bases (`total_gb`, tumour plus normal) and the tumour and normal coverage of each sample, as given in Table 1 of the manuscript. Notebook 03 prints them under the panels of its figure. They are not computed from anything in the repository. Checked against NCBI's run information: the summed bases of the two runs of LL, FD, NC and EA round to the tabulated values; for IL they sum to 65.88 Gb (Table 1: 65.8). Coverage cannot be checked there.

### 3.2 Pipeline configuration

`reference/TestCases.csv` has one row per run; `TestCaseNo` is the run number and the name of its folder. The 480 runs of notebook 01 are the full factorial 5 x 2 x 2 x 2 x 2 x 3 x 2, one run per combination:

| Column | Levels |
|---|---|
| `Sample` | EA, FD, IL, LL, NC |
| `HOW` (environment) | Altay + COSAP, Uhem + COSAP |
| `isTrimmed` | YES, NO |
| `baseRecalibration` | YES, NO |
| `Mapper` | BWA, BOWTIE |
| `Variant Caller` | Mutect, Strelka, SomaticSniper |
| `Duplicates` | MARK, DELETE |

All 480 runs use the hg38 reference. `TestCases.csv` has 864 rows; only the rows that have a folder in `raw/snv/` are runs of this study, and `TestCases_Indels.csv` (13 columns, run key `IDs` instead of `TestCaseNo`) describes the 320 indel runs.

**What is not recorded.** `M_Version` and `VC_Version` (tool versions) are empty for all 480 runs, and so are `Variant Filter`, `FilterPass`, `SOMATICfilter` and `remove_indels`. The repository therefore does not fix tool versions, and regenerated VCFs are not expected to be byte-identical to the originals. The acceptance test is notebook 01's parity cell: the number of variants of every run, after the legacy filters, must equal `Result VCF Count` in `derived/vcfcomparison_df_full.csv` (480/480, otherwise the notebook fails and writes nothing).

What notebook 01 uses from a regenerated VCF: the FILTER, REF and ALT columns. Mutect and Strelka records must carry `PASS` in FILTER, since 01 keeps only those; SomaticSniper is not PASS-filtered (its FILTER is `.` on every record). To reproduce the original counts the calls must also cover the same regions as the originals (their file names end in `.secondbed.vcf.recode.vcf`); the repository does not say which region file was used. Analyses that need the tumour sample read it from the header (`##tumor_sample=` in Mutect VCFs, the columns `NORMAL` and `TUMOR` in Strelka and SomaticSniper VCFs), see D5.

## 4. Directory layout expected under `raw/`

```
data/raw/
  snv/
    TestCase 296/
      snp_SRR7890919-SRR7890919bwa2-SRR7890918-SRR7890918bwa2_mutect.vcf.secondbed.vcf.recode.vcf
    TestCase 853/
      snp_SRR7890874_bowtie-SRR7890883_bowtie_strelka.vcf.secondbed.vcf.recode.vcf
    TestCase 606/
      snp_SRR7890880-SRR7890880bowtie-SRR7890879-SRR7890879bowtie_somaticsniper.vcf.filtered.recode.vcf
    ...                                   480 folders in all
  indel/
    135/
      indel_SRR7890919_bowtie-SRR7890918_bowtie_mutect.vcf.secondbed.vcf.recode.vcf
    ...                                   320 folders in all
```

Notebook 01 enforces:
- `raw/snv/` exists and holds exactly **480** run folders, each named `TestCase <TestCaseNo>` with a `TestCaseNo` from `TestCases.csv`; otherwise it stops with a message pointing here;
- each folder holds exactly one file ending in `.vcf`;
- the file name contains the run's caller (`mutect`, `strelka` or `somaticsniper`, in any case).

Other files are ignored: `.vcf.gz` copies beside the VCFs (the author's disk has them; you do not need them), `.DS_Store`, and every name starting with `._` (macOS AppleDouble files, which parse as empty).

Notebook 00 goes further when `raw/` is present: it compares the file lists with `derived/vcf_filenames.csv` (SNV) and `derived/vcf_filenames_indels.csv` (indel), stripping the legacy `vcf/TESTCASES_bedded/` and `vcf/TESTCASES_bedded_indels/` prefixes. Keep the original file names if you want that check to pass. No notebook parses `raw/indel/` yet; it is only inventoried.

To use the copies in the repository:

```bash
cp -R data/raw/vcf_snps data/raw/snv       # 476 folders; add TestCase 26, 32, 38, 44 before running notebook 01
cp -R data/raw/vcf_indels data/raw/indel   # 320 folders (the 10 extra .vcf.gz are ignored)
```

## 5. Verifying the files

Notebook 00 writes the size and SHA-256 of every file to `results/00_data_preparation/tables/`:
- `inventory_committed.csv`: `truth/`, `reference/` and the committed files of `derived/` (every committed file);
- `inventory_local.csv`: `raw/snv/`, `raw/indel/` and `derived/sets_dict.parquet` (1,361 files).

Commit both with the results so that they travel with the repository. Run this from the repository root (Python 3, standard library only). It checks your files against the author's checksums, skips `*.parquet` (rebuilt by notebook 01; its bytes can differ between pyarrow versions), and reports files that differ or are absent. Absent raw files are expected if you do not have `raw/`.

```python
import csv
import hashlib
import os


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


for table in ("inventory_committed", "inventory_local"):
    match, differ, absent = 0, [], []
    with open("results/00_data_preparation/tables/%s.csv" % table, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            path = row["path"]
            if path.endswith(".parquet"):
                continue
            if not os.path.isfile(path):
                absent.append(path)
            elif os.path.getsize(path) != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
                differ.append(path)
            else:
                match += 1
    print("{}: {} match, {} differ, {} absent".format(table, match, len(differ), len(absent)))
    for label, paths in (("DIFFERS", differ), ("ABSENT", absent)):
        for path in paths[:5]:
            print("  {} {}".format(label, path))
```

Expected on a correct checkout: `inventory_committed: N match, 0 differ, 0 absent`, N being the number of rows in that CSV, and for `inventory_local` either `1360 match, 0 differ, 0 absent` (all raw files present) or `0 differ` with the raw files reported absent. Without `git-lfs` the six `giab_regions/*.bed` files show as `DIFFERS`; run `git lfs pull` (no notebook reads them). On Windows set `git config core.autocrlf false` before cloning, or line-ending conversion changes the hashes of the text files. For a single file, compare `shasum -a 256 <file>` (macOS) or `sha256sum <file>` (Linux) with its row in the CSV.

You can also rerun notebook 00: `git diff --stat results/00_data_preparation/tables/inventory_committed.csv` prints nothing if the committed files are unchanged.
