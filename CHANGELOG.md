# nickjhathaway/elucidatornf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nickjhathaway/elucidatornf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--amplicon_samples_keep_file` for `illumina_clustering.nf`, `nanopore_clustering.nf` and `illumina_read_qc.nf`: only process the samples listed (input names or sample key sample/replicate names), like `--wgs_samples_keep_file` and `--pw_samples_file`

### `Fixed`

- The extractor flag guessing (`AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR`) only looks at the input files the run processes (kept samples, properly paired/named files), so other files in `--input_fastq_dir`, e.g. single end fastqs, no longer make it fail

### `Dependencies`

### `Deprecated`
