# Sample key: renaming, replicates and MIDs

`illumina_clustering.nf` and `nanopore_clustering.nf` take an optional sample key with
`--amplicon_sample_key_fnp`. It maps the **input names** the pipeline finds to the **sample** and
**replicate** names used for clustering, and can split pooled inputs by MID (barcode).

Without a key, every input is its own sample with one replicate, both named after the input.

## Format

Tab (or comma) delimited, with a header. Blank lines and lines starting with `#` are skipped.

| column | required | meaning |
|---|---|---|
| `input` | yes | the input name (see [Input names](#input-names)) |
| `sample` | yes | the final sample name |
| `replicate` | no | the replicate name; filled in automatically when blank |
| `fbarcode` | no | forward MID of this sample within the input |
| `rbarcode` | no | reverse MID, only for dual barcoding (needs `fbarcode`) |

```tsv
input	sample	replicate	fbarcode	rbarcode
Run1_LibA_S1	P001		ACGTACGTAC	TTGACCATGC
Run1_LibA_S1	P002		TGCATGCATG	CCAGTTGACA
Run2_LibA_S3	P001		ACGTACGTAC	TTGACCATGC
Sample9_S4	P009
```

Here `P001` has two replicates (one in each run), `P002` and `P009` have one, and `Run1_LibA_S1`
holds two pooled samples separated by their MIDs.

### Old rename format

A headerless 2 column file (`old_name  new_name`) still works and is a plain rename: each input
becomes one sample with one replicate, both named `new_name`. `--nanopore_rename_key_fnp` is the
older name of the parameter; give only one of it and `--amplicon_sample_key_fnp`.

## Input names

| | input name |
|---|---|
| Illumina | the fastq pair prefix: `Run1_LibA_S1` for `Run1_LibA_S1_R1_001.fastq.gz` |
| Nanopore | the name in `meta/rawInputSampleNames.tsv` (the bam `SM`, e.g. `barcode01`, or the file name; prefixed with the top directory when `--nanopore_sample_name_prefix_top_dir` is on) |

Key entries that don't match any input are warned about and ignored. For nanopore the warning
points out entries that only match before the top directory prefix. To check names ahead of time:

```bash
SeekDeep getPossibleSampleNamesFromRawInput --inputDir demux_dir --technology nanopore --recursive --prefixWithTopDir
```

## How replicates are used

1. **Extraction** runs once per input. Inputs with barcodes in the key are demultiplexed with only
   their own MIDs (an id file per input in `meta/midIdFiles/`), and the extractor writes one file
   per MID named by replicate, `<target><replicate>.fastq.gz`.
2. **Per-sample clustering** (`SeekDeep clusterDown`, `AMPLICON_CLUSTER_BY_EXPECTED_ERRORS` for
   Illumina, `AMPLICON_CLUSTER_BY_KMER_SIMILARITY` for nanopore) runs once per **replicate** and
   target, tagged with `--sample` and `--replicate`.
3. **Population clustering** (`SeekDeep processClusters`) groups the replicates by sample:
   - a replicate needs `--illumina_clustering_min_sample_read_count` /
     `--nanopore_clustering_min_sample_read_count` reads to be used;
   - by default a haplotype must be found in **every** passing replicate of a sample to be kept.
     Add `--runsRequired <n>` to `--amplicon_clustering_pop_clustering_extra_args` to relax this.

So replicates of one sample are clustered separately and combined at the population step. Inputs
that are **not** replicates should get different `sample` names; giving two inputs the same
`sample` means their haplotypes are intersected, not pooled.

## Replicate names

When `replicate` is blank:

- the sample has one row → the replicate is the sample name;
- the sample has several rows → `<input>` plus the barcodes (`<input>-<fbarcode>-<rbarcode>`),
  which can't collide.

Name them yourself if you want readable per-replicate outputs (e.g. `P001-rep1`, `P001-rep2`).

## Validation rules

The run stops before extraction if any of these fail:

- an `input` + barcodes combination appears only once;
- replicate names are unique, and no replicate is named after a different sample;
- within an input, either every row has barcodes or none do, and an input without barcodes has
  only one row (more rows would need MIDs to separate them);
- an input doesn't mix single and dual barcoded rows;
- names have no `/`, `\`, `?` or whitespace; barcodes are only `A`, `C`, `G`, `T`, `N`;
- only the columns above are allowed.

An input left out of the key can't share a name with a sample or replicate in the key (it would be
read as that sample), so add it to the key.

## Outputs

- `meta/sample_name_key.tsv`: input, sample, replicate, barcodes and `extraction_name` for every
  sample/replicate that had extracted reads.
- `meta/midIdFiles/`: the per-input MID id files, and `allMids_ids.tsv` (every barcode in the key,
  used to guess the extractor flags).
- The extraction reports are per input; in the HTML reports they're labelled by replicate when the
  input holds a single sample, and per-MID counts replace per-input counts for barcoded inputs.

## Illumina vs nanopore

The key is read and applied the same way in both. Differences:

- **Input names**: fastq pair prefix vs scanned name (see above).
- **Barcode length (Illumina)**: MIDs take up read length, so the longest `fbarcode` + `rbarcode`
  in the key is used when deciding which targets overlap enough to stitch. Targets that can't be
  stitched are clustered as `<target>-R1` and `<target>-R2`, each with its own replicates.
