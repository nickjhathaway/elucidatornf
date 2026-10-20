# RNA-seq results

Written by elucidatornf `rnaseq_process_pf.nf` (version @PIPELINE_VERSION@) on @RUN_DATE@.
This directory holds @N_SAMPLES@ samples.

**Start with `summary/rnaseq_summary.tsv.gz`**, which has one row per sample covering read counts
through every step, QC, and Salmon mapping rates.

## Directories

| Directory | Contents |
| --- | --- |
| `summary/` | `rnaseq_summary.tsv.gz` for all samples; `<sample>_rnaseq_summary.tab.txt` per sample |
| `fastp_and_filter_info/` | fastp QC (`<sample>.fastp.json.gz`); host-filter read counts per pass (`_totalReadCounts`) and per contig (`_filteredByChrom`) |
| `human_removed_fastqs/` | reads left after host removal: `<sample>_kept_R1/_R2.fastq.gz` |
| `initial_quants/` | Salmon on the trimmed reads, **before** host removal; `all_quants.tsv.gz` combines all samples |
| `filtered_quants/` | Salmon on the host-removed reads; `all_quants.tsv.gz`; stage proportions from MuSiC (`pf_music_decomp_cell_proportions.tsv`) and InstaPrism (`pf_instaprism_fine_cell_proportions.tsv`, which also gives each fine stage's `coarse_stage`, and `pf_instaprism_course_cell_proportions.tsv`) |
| `kmer_extraction/` | k-mer read counts combined over all samples: `kmer_filter_*` (removing mycoplasma and Pf rRNA) and `kmer_sort_*` (sorting reads by gene) |
| `kmer_extraction/counts/` | the same count tables, per sample |
| `kmer_extraction/sorted_fastqs/<sample>/` | reads sorted by gene: `<gene>_R1/_R2.fastq.gz` pairs and `<gene>.fastq.gz` unpaired reads |
| `pipeline_info/` | Nextflow execution report, timeline and trace for each run |

## Reading the summary

- Counts are in **reads**, so each mate counts once. The Salmon columns count fragments (read pairs).
- `input_reads` = `fastp_removed_reads` + `host_reads` + `non_host_reads`.
- `host_reads` are the reads the two host filters removed: everything mapping to the listed host contigs (e.g. human and mycoplasma).
- `pf_reads` are the non-host reads the k-mer filter assigns to Pf (sets: @PF_SETS@).
- `sorted_gene_reads` and the per-gene columns (@GENE_SETS@) come from the k-mer sort. `kmer_extraction/kmer_sort_extractionCounts.tsv.gz` has every gene.

## Adding samples later

Rerunning into this directory skips every sample whose outputs are already here, and rebuilds the
combined tables, the deconvolution and the summary over all samples. After changing a reference or
a setting, use a new `--outdir` or pass `--rnaseq_skip_existing false`, so old and new results aren't mixed.

## References used in the latest run

| Input | Path |
| --- | --- |
| Reference directory (`--rnaseq_ref_files_dir`) | `@REF_FILES_DIR@` |
| Host 1 genome / contigs removed | `@HOST1_GENOME@` / `@HOST1_CONTIGS@` |
| Host 2 genome / contigs removed | `@HOST2_GENOME@` / `@HOST2_CONTIGS@` |
| Salmon index | `@SALMON_INDEX@` |
| k-mer filter table | `@KMER_FILTER_TABLE@` |
| k-mer sort table | `@KMER_SORT_TABLE@` |
| Non-unique k-mer tables used | @NON_UNIQUE_KMER_TABLES@ |
| Single-cell object | `@SINGLE_CELL@` |
