# Changelog

All notable changes to PhageFlow are documented here.
This project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- **annotate / resistance**: outputs are now scoped by sample
  (`results/06_annotation/{sample}/{candidate}/`, same for `07_resistance`).
  Candidate IDs are numbered per sample, so in multi-sample runs
  `X_candidate_001` of the second sample overwrote — or, because of the
  skip-if-exists logic, silently reused — the first sample's annotation.
  Aggregate TSVs are keyed by (sample, candidate). `resistance` still reads
  the old layout.
- **viral-id**: geNomad is run with `--min-score = min(min_score, rescue_min_score)`.
  Previously geNomad dropped every contig below 0.7, so the 0.4–0.7 rescue tier
  could never contain anything.
- **viral-id**: genetic-code fallback from `virus_genes.tsv` read a non-existent
  column and split contig names on `|`, so it never matched.
- **coverage**: CoverM "Trimmed Mean" column was parsed as "Mean" (suffix match),
  so mean coverage was actually the trimmed mean and trimmed mean was always 0.
- **qc**: `est_coverage_50kb` in `qc_summary.tsv` was always empty (wrong key).
- **host-removal**: shell pipelines (bwa-mem2 | samtools …, read merging) now
  run with `bash -o pipefail` and quoted paths. A crash of bwa-mem2/seqtk was
  previously ignored and left truncated FASTQs that the next run skipped as done;
  paths with spaces broke the commands.
- **host-removal**: Kraken2 failures are no longer ignored (`check=False` removed).
- **host-removal**: Level B post-filter downloads reference genomes only for
  species-rank taxa (it used every taxon under Bacteria, e.g. whole phyla);
  Level A reports contaminants at genus rank using Kraken2 rank codes.
- **host-removal**: a missing Kraken2 DB raised `SystemExit`, which aborted the
  whole multi-sample `run` instead of failing one sample.
- **quality**: read-recruitment pipeline also uses `pipefail`; file handle leak fixed.
- **cli**: R2 discovery replaced *every* "R1" in the filename, so SRA
  accessions such as `SRR1234_R1` looked for `SRR2234_R2`. Only the mate token
  is replaced now.
- **cli**: removed `safety` and `report` commands, which imported modules that
  do not exist; `--from-module` choices now match the runner
  (`coverage` and `resistance` added, `safety`/`report` removed — they silently
  restarted the pipeline from `qc`).
- **cli**: `-o` is optional with `--project` (the documented
  `phageflow run --project DIR` failed with "missing option -o").
- **config**: `PHAGEFLOW_DB` was always overridden by the bundled config's
  relative database paths, so it never took effect; precedence now follows the
  README. Empty/`null` YAML sections no longer crash. Relative paths are resolved
  from the project root consistently (also for `kraken2`). CheckV version
  directory under `PHAGEFLOW_DB` is auto-detected.
- **download-databases**: Kraken2 tarball extracted with `filter="data"` (PEP 706).
- **assembly**: singletons are also given to MEGAHIT (`-r`), as documented;
  MEGAHIT honours `assembly.min_length`; progress labels no longer claim
  `--isolate`/`--no-mercy` modes that are not used.
- **release workflow**: the YAML was invalid (unindented code fences inside the
  release body), so the workflow could never run.

### Changed / workflow
- `phageflow run` writes `reports/pipeline_status.tsv` (running / completed /
  failed / skipped per sample and module); `phageflow status` accepts `-o`/`--reports-dir`.
- `phageflow run` exits non-zero when any sample fails; duplicate sample IDs are
  rejected up front; resuming with `--from-module` checks that required inputs exist.
- QC uses the fastp JSON as completion marker and deletes partial outputs on failure.
- New pytest suite (`tests/`) and ruff config; CI runs lint, YAML validation,
  tests on Python 3.10–3.12, all CLI `--help` pages and a wheel build.
- Release workflow verifies that the tag matches `__init__.py`, `pyproject.toml`
  and the conda recipe, and retries the tarball download instead of a fixed sleep.
- `environment.yml` minimum versions aligned with `conda-recipe/meta.yaml`.

## [1.0.0] - 2026-06-01

First production release. Complete 7-module pipeline validated end-to-end on
purified phage, virome, and metagenome datasets, including a multi-phage
recovery (4 distinct genomes from a single K. pneumoniae sample) on an HPC cluster.

### Pipeline
- **01 qc**: fastp + FastQC + MultiQC, read-length auto-detection (PE150/250/300)
- **02 host-removal**: bwa-mem2, Kraken2, and pass-through modes
- **03 assembly**: SPAdes + MEGAHIT + cd-hit-est, k-mer range by read length
- **03b coverage**: CoverM per-contig profiling
- **04 viral-id**: geNomad with 3 kb rescue tier, CrAss-like genetic-code detection
- **05 quality**: CheckV tiers, blastn/minimap2 dereplication, circular-rotation detection
- **06 annotate**: Pharokka + Phold two-tier cascade (--hyps strategy)
- **07 resistance**: dual-evidence AMR/VFDB + ACR + DefenseFinder + NetFlax
- **runner**: full pipeline with --from-module resume and per-sample error handling

### Infrastructure
- `phageflow download-databases`: one-command setup for CheckV, geNomad, Pharokka,
  Phold (incl. GPU structure DB), and Kraken2 standard-16GB; writes user config
- Config hierarchy: CLI > project > ~/.config/phageflow/ > bundled; `PHAGEFLOW_DB` env var
- Available on Bioconda: `conda install -c bioconda phageflow`

### Robustness fixes since initial bioconda submission
- CLI now loads user config from ~/.config/phageflow/config.yaml
- geNomad nested database path (genomad_db/genomad_db) resolved automatically
- Phold GPU database built during download-databases
- Phold auto-retries on CPU when Foldseek GPU search fails (no-GPU clusters)
- Phold falls back to Pharokka GBK for tiny/divergent genomes with no structural hits
