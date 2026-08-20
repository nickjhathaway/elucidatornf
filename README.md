# nickjhathaway/elucidatornf

## Introduction

**nickjhathaway/elucidatornf** is a collection of Nextflow workflows that stream
together tools from [SeekDeep](https://github.com/bailey-lab/SeekDeep),
[PathWeaver](https://github.com/bailey-lab/PathWeaver), and
[elucidator](https://github.com/nickjhathaway/elucidator) for amplicon
clustering, targeted assembly, and *P. falciparum* WGS/RNA-seq preprocessing.

This is a **multi-workflow project**: there is no single default pipeline.
Each analysis is its own entry script in the repository root — run the one you
need directly. `main.nf` is only a helper that lists the available workflows
(`nextflow run main.nf` prints the list below; it does not run an analysis).

| Entry script | What it does |
| --- | --- |
| `illumina_clustering.nf` | Illumina amplicon clustering (SeekDeep) |
| `nanopore_clustering.nf` | Nanopore amplicon clustering (SeekDeep) |
| `pathweaver_extract_regions.nf` | PathWeaver targeted region extraction + population clustering |
| `wgs_preprocess_pf.nf` | WGS read preprocessing / host filtering (*P. falciparum*) |
| `wgs_host_depletion_pf.nf` | Host read removal only, for submitting otherwise-raw reads |
| `rnaseq_process_pf.nf` | RNA-seq processing + deconvolution (*P. falciparum*) |
| `index_genomes.nf` | Build indexes for a directory of genomes |

## Usage

> [!NOTE]
> If you are new to Nextflow, see [this page](https://nf-co.re/docs/usage/installation)
> on how to set up Nextflow.

> [!WARNING]
> On Nextflow 26.04+ (strict syntax / parser v2), parameters given on the CLI are
> always strings. In practice the only thing to watch is **booleans**: to turn a
> flag *on*, pass it bare (e.g. `--do_variant_calling`); to turn it *off*, simply
> **omit it** — do not pass `--do_variant_calling false` (the string `"false"` is
> truthy and would still enable it). For fully typed inputs, use a
> [`-params-file`](https://www.nextflow.io/docs/latest/cli.html#run) (YAML/JSON).

### Common options

These apply to every workflow:

- `--outdir <dir>` — output directory (**required** by all but `index_genomes.nf`).
- `-profile <docker/singularity/apptainer/...>` — container engine / institutional profile.
- Resources (see [`conf/params/global.config`](conf/params/global.config)):
  `--max_task_cpus` / `--max_task_memory` / `--max_task_time` cap a single task;
  `--max_executor_cpus` / `--max_executor_memory` cap the whole-run pool;
  both are clamped to the detected machine. `--max_retry` sets automatic retries.

The amplicon workflows (`illumina`, `nanopore`) and `pathweaver` share an optional
**variant-calling** stage (see [`conf/params/variant_calling.config`](conf/params/variant_calling.config)):

- `--do_variant_calling` — enable variant calling (off by default).
- `--vc_primary_genome <name>`, `--meta_fnp <tsv>`, `--variant_calling_ncpus <n>`,
  `--vc_variant_frequency_cut_off`, `--vc_variant_occurrence_cut_off`,
  `--meta_fields_to_calc_pop_diffs`, `--vc_extra_args "..."`.

### Illumina amplicon clustering — `illumina_clustering.nf`

**Required:** `--input_fastq_dir`, `--primers_fnp`, `--genome_dir`, `--gff_dir`,
`--outdir`, `--illumina_clustering_paired_end_length`.

**Key options** (see [`conf/params/illumina.config`](conf/params/illumina.config) /
[`conf/params/amplicon.config`](conf/params/amplicon.config)):
`--illumina_clustering_min_sample_read_count`, `--illumina_clustering_trim_front`,
`--illumina_clustering_trim_back`, `--illumina_clustering_population_ncpus`,
`--amplicon_clustering_pop_clustering_extra_args`.

```bash
nextflow run elucidatornf/illumina_clustering.nf \
      --input_fastq_dir fastq/ \
      --primers_fnp primers.tsv \
      --genome_dir /tank/data/plasmodium/genomes/pf/genomes/ \
      --gff_dir /tank/data/plasmodium/genomes/pf/info/gff/ \
      --illumina_clustering_paired_end_length 250 \
      --outdir analysis_illumina \
      -profile <docker/singularity/apptainer/.../institute>
```

### Nanopore amplicon clustering — `nanopore_clustering.nf`

**Required:** `--input_fastq_dir`, `--primers_fnp`, `--genome_dir`, `--gff_dir`,
`--outdir`.

**Key options** (see [`conf/params/nanopore.config`](conf/params/nanopore.config) /
[`conf/params/amplicon.config`](conf/params/amplicon.config)):
`--nanopore_clustering_min_len`, `--nanopore_clustering_min_sample_read_count`,
`--nanopore_clustering_max_reads_use`, `--nanopore_extractor_use_inner_primers`,
`--nanopore_primers_errors_allowed`, `--nanopore_rename_key_fnp`,
`--nanopore_render_clustering_report`.

```bash
nextflow run elucidatornf/nanopore_clustering.nf \
      --input_fastq_dir fastq/ \
      --primers_fnp primers.tsv \
      --genome_dir /tank/data/plasmodium/genomes/pf/genomes/ \
      --gff_dir /tank/data/plasmodium/genomes/pf/info/gff/ \
      --outdir analysis_nanopore \
      --do_variant_calling \
      --vc_primary_genome Pf3D7 \
      -profile <docker/singularity/apptainer/.../institute>
```

### PathWeaver targeted region extraction — `pathweaver_extract_regions.nf`

**Required:** `--bams_dir`, `--genome_fnp`, `--outdir`, **and exactly one** way to
specify the target regions:

- `--pw_bed_fnp <bed>` — regions from a BED file, **or**
- `--primers_fnp <tsv>` — regions from a primers file, **or**
- `--pw_gene_ids <ids>` — regions from gene IDs, **or**
- `--pw_seqs_table <tsv>` (with `--pw_seqs_table_seqs_col`, `--pw_seqs_table_name_col`,
  `--pw_seqs_table_target_col`) — regions from a sequences table.

**Key options** (see [`conf/params/pathweaver.config`](conf/params/pathweaver.config)):
`--pw_samples_file`, `--bams_file_ending` (default `.sorted.bam`),
`--pw_overwrite_dir`, `--render_pw_report`, `--pw_ncpus`, `--pw_pop_clustering_ncpus`,
`--run_sub_segments_determination`, `--pw_pop_clustering_extra_args`.

```bash
nextflow run elucidatornf/pathweaver_extract_regions.nf \
      --bams_dir /data/pf/bams \
      --bams_file_ending .sorted.bam \
      --pw_bed_fnp regions.bed \
      --genome_fnp pf_genomes/genomes/Pf3D7.fasta \
      --outdir pathweaver_results \
      --do_variant_calling \
      --meta_fnp /data/pf/metadata/meta.tab.txt \
      --render_pw_report \
      -profile <docker/singularity/apptainer/.../institute>
```

### WGS preprocessing — `wgs_preprocess_pf.nf`

**Required:** `--input_fastq_dir`, `--outdir`, `--wgs_host1_filter_genome_fasta_fnp`,
`--wgs_host1_filter_contigs_fnp`, `--wgs_host2_filter_genome_fasta_fnp`,
`--wgs_host2_filter_contigs_fnp`, `--wgs_final_genome_fasta_fnp`.

**Key options** (see [`conf/params/wgs.config`](conf/params/wgs.config) /
[`conf/params/wgs_common.config`](conf/params/wgs_common.config)):
`--wgs_samples_keep_file`, `--wgs_rename_tsv`, `--wgs_skip_existing`,
`--wgs_host_filter_min_mapq`, `--wgs_host_filter_any_mate`,
`--wgs_host_filter_filter_with_unmapped_mate`.

```bash
nextflow run elucidatornf/wgs_preprocess_pf.nf \
      --input_fastq_dir fastq/ \
      --wgs_host1_filter_genome_fasta_fnp human.fasta \
      --wgs_host1_filter_contigs_fnp human_contigs.txt \
      --wgs_host2_filter_genome_fasta_fnp anopheles.fasta \
      --wgs_host2_filter_contigs_fnp anopheles_contigs.txt \
      --wgs_final_genome_fasta_fnp Pf3D7.fasta \
      --outdir wgs_preprocessed \
      -profile <docker/singularity/apptainer/.../institute>
```

### Host read removal only — `wgs_host_depletion_pf.nf`

The host-filter core of `wgs_preprocess_pf.nf` with **no fastp trimming and no final
mapping**, for submitting data that should be raw apart from having host reads stripped.
Reads that survive come back off the bams in their original orientation, so they are
byte identical (name, sequence, quality) to the input reads. Any read whose mate was
filtered off is discarded so the output stays strictly paired end.

**Required:** `--input_fastq_dir`, `--outdir`, `--wgs_host1_filter_genome_fasta_fnp`,
`--wgs_host1_filter_contigs_fnp`, `--wgs_host2_filter_genome_fasta_fnp`,
`--wgs_host2_filter_contigs_fnp`.

**Key options:** shares `--wgs_samples_keep_file`, `--wgs_rename_tsv`,
`--wgs_skip_existing` and the host-filter options in
[`conf/params/wgs_common.config`](conf/params/wgs_common.config), plus
`--wgs_depletion_save_unpaired` (default `false`) to keep the discarded singletons.

Outputs, under `--outdir`:

| Path | Contents |
| --- | --- |
| `host_depleted_fastqs/<sample>_R{1,2}.fastq.gz` | the reads to submit |
| `host_filter_info/host_depletion_summary.tsv.gz` | per sample: input reads, reads filtered off per host, singletons discarded, reads remaining |
| `host_filter_info/<sample>_filteredByChrom.tsv.gz` | which host contigs the filtered reads went to |
| `host_filter_info/<sample>_totalReadCounts.tsv.gz` | raw counts per filter pass |
| `host_filter_info/sample_name_key.tsv` | old → new sample name key from `--wgs_rename_tsv` |
| `discarded_singletons/<sample>_unpaired.fastq.gz` | only with `--wgs_depletion_save_unpaired` |

```bash
nextflow run elucidatornf/wgs_host_depletion_pf.nf \
      --input_fastq_dir fastq/ \
      --wgs_host1_filter_genome_fasta_fnp human.fasta \
      --wgs_host1_filter_contigs_fnp human_contigs.txt \
      --wgs_host2_filter_genome_fasta_fnp anopheles.fasta \
      --wgs_host2_filter_contigs_fnp anopheles_contigs.txt \
      --wgs_rename_tsv rename_key.tsv \
      --outdir wgs_host_depleted \
      -profile <docker/singularity/apptainer/.../institute>
```

### RNA-seq processing — `rnaseq_process_pf.nf`

**Required:** `--input_fastq_dir`, `--outdir`, `--rnaseq_host1_filter_genome_fasta_fnp`,
`--rnaseq_host1_filter_contigs_fnp`, `--rnaseq_host2_filter_genome_fasta_fnp`,
`--rnaseq_host2_filter_contigs_fnp`, `--rnaseq_salmon_index`,
`--single_cell_seurat_r_object_fnp`.

**Key options:** shares the host-filter options in
[`conf/params/wgs_common.config`](conf/params/wgs_common.config).

```bash
nextflow run elucidatornf/rnaseq_process_pf.nf \
      --input_fastq_dir fastq/ \
      --rnaseq_host1_filter_genome_fasta_fnp human.fasta \
      --rnaseq_host1_filter_contigs_fnp human_contigs.txt \
      --rnaseq_host2_filter_genome_fasta_fnp anopheles.fasta \
      --rnaseq_host2_filter_contigs_fnp anopheles_contigs.txt \
      --rnaseq_salmon_index salmon_index/ \
      --single_cell_seurat_r_object_fnp reference.rds \
      --outdir rnaseq_processed \
      -profile <docker/singularity/apptainer/.../institute>
```

### Index genomes — `index_genomes.nf`

**Required:** `--genome_dir` (a directory of genome FASTA files to index).

```bash
nextflow run elucidatornf/index_genomes.nf \
      --genome_dir /tank/data/plasmodium/genomes/pf/genomes/ \
      -profile <docker/singularity/apptainer/.../institute>
```

## Shell completion (bash)

A bash completion is provided under [`completions/`](completions/). Once enabled,
it detects the workflow by the `.nf` basename (so full or partial paths both
work) and tab-completes only that workflow's `--flags`, with file/dir completion
for path-valued flags and `-profile` values. For any command that is **not** one
of these workflows (other pipelines, `nextflow log`/`pull`/…), it delegates to
nextflow's own completion if one is installed — so source this **after** nextflow's
completion (e.g. after the `bash-completion` package) for that hand-off to work.

Enable it with one of:

```bash
# source from your shell startup
echo 'source /path/to/elucidatornf/completions/elucidatornf.bash' >> ~/.bashrc

# or install into the user completion dir (auto-loaded)
cp completions/elucidatornf.bash ~/.local/share/bash-completion/completions/nextflow
```

Then open a new shell (or `source ~/.bashrc`). Example:

```bash
nextflow run /path/to/elucidatornf/pathweaver_extract_regions.nf --<TAB>
# -> --bams_dir --pw_bed_fnp --genome_fnp --do_variant_calling --outdir ...
```

The flag lists mirror `conf/params/*.config`; if you add params, refresh them
with `./completions/generate_completion.sh`.

## Credits

nickjhathaway/elucidatornf was originally written by Nicholas Hathaway.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
