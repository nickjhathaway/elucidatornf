# nickjhathaway/elucidatornf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nickjhathaway/elucidatornf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--amplicon_samples_keep_file` for `illumina_clustering.nf`, `nanopore_clustering.nf` and `illumina_read_qc.nf`: only process the samples listed (input names or sample key sample/replicate names), like `--wgs_samples_keep_file` and `--pw_samples_file`
- `--download_plasmodium_genomes` for `illumina_clustering.nf`, `nanopore_clustering.nf` and `pathweaver_extract_regions.nf`: download, gunzip and index [plasmodium_genomes](https://github.com/nickjhathaway/plasmodium_genomes) into `--plasmodium_genomes_dir` (default `<outdir>/plasmodium_genomes`) and use it for `--genome_dir`/`--gff_dir` (pathweaver: `--genome_fnp` names the genome, default Pf3D7); an already extracted and indexed directory is reused
- `download_plasmodium_genomes.nf`: just the download/gunzip/index of plasmodium_genomes into `--plasmodium_genomes_dir`, to set up a shared directory once
- `--illumina_single_end` for `illumina_clustering.nf`: one fastq per sample rather than pairs (e.g. SRA runs merged before upload), extracted with `SeekDeep extractor`, whose reports now match `extractorPairedEnd`'s (needs a SeekDeep with that change); `--illumina_clustering_paired_end_length` isn't needed with it
- `illumina_clustering.nf` and `nanopore_clustering.nf` stop before extraction when a sample isn't in `--meta_fnp` (or the file has no `sample` column or repeats a sample), rather than failing in population clustering

- `--nanopore_trim_primers_fnp` for `nanopore_clustering.nf`: trim the reads of the amplified targets (`--primers_fnp`) to one or more regions per target, each marked off by its own primers
- `nanopore_clustering.nf` no longer guesses `--primerWithinStart` from the input (`AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR`, still used for Illumina), reads are cut to their target before the primers are looked for so it isn't needed; `--nanopore_extractor_primer_within_start` still sets it. There's no `extractionReports/outSeekDeepExtractorFlags.txt` for nanopore runs anymore

### `Fixed`

- The extractor flag guessing (`AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR`) only looks at the input files the run processes (kept samples, properly paired/named files), so other files in `--input_fastq_dir`, e.g. single end fastqs, no longer make it fail
- Resetting `--outdir` no longer fails with `DirectoryNotEmptyException` when it's on NFS, and leaves the execution trace of the run in progress in `pipeline_info` instead of deleting it

### `Dependencies`

### `Deprecated`
