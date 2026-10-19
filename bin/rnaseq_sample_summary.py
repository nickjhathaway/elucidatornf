#!/usr/bin/env python3
"""Per-sample read accounting and QC summary for rnaseq_process_pf.nf.

Every count is in reads (each mate counted), which is how fastp, BamFilterByChroms
and extractByCountingUniqKmersFromSets all count, so the stages add up:

    input_reads = fastp_removed_reads + host_reads + non_host_reads
    non_host_reads = sum(filter_<set>_reads) + filter_undetermined_reads + filter_other_reads

host_reads is whatever the two host filters removed: the contigs listed for each
host, e.g. human plus mycoplasma.

salmon counts fragments (pairs), so its counts are labelled as fragments.
Writes one header line and one row, tab separated, to --out.
Standard library only, since the pipeline image has no conda python.
"""

import argparse
import csv
import gzip
import json
import sys

# rows of an extractionCounts table that are bookkeeping rather than k-mer sets
KMER_SPECIAL_ROWS = {"undetermined", "smallLenReads", "poorQualityCount", "containsNs", "inverseCounts", "multiHitCounts"}
FASTP_REMOVED_REASONS = ["low_quality_reads", "too_many_N_reads", "too_short_reads", "too_long_reads",
                         "low_complexity_reads", "adapter_dimer_reads"]


def open_text(fnp):
    return gzip.open(fnp, "rt") if str(fnp).endswith(".gz") else open(fnp)


def read_json(fnp):
    with open_text(fnp) as fh:
        return json.load(fh)


def read_tsv(fnp):
    with open_text(fnp) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def comma_list(value):
    return [v.strip() for v in value.split(",") if v.strip()]


def pct(num, den):
    return f"{100.0 * num / den:.4f}" if den > 0 else "NA"


def host_counts(rows):
    """BamFilterByChroms totalReadCounts, combined over the four filter passes."""
    filtered = {"host1": 0, "host2": 0}
    kept = 0
    for row in rows:
        count = int(float(row["count"]))
        if row["condition"].startswith("filtered"):
            host = row["filter"].split("_", 1)[0]
            filtered[host] = filtered.get(host, 0) + count
        elif row["condition"].startswith("kept"):
            # the last pass does not write unmapped reads separately, so its kept
            # reads include them; together the kept reads are what goes on
            kept += count
    return filtered, kept


def kmer_counts(rows):
    """Counts per target from the last iteration of an extractionCounts table."""
    iterations = [row["iteration"] for row in rows]
    last = "final" if "final" in iterations else "initial"
    counts = {}
    total = 0
    for row in rows:
        if row["iteration"] != last:
            continue
        counts[row["target"]] = int(float(row["count"]))
        total = int(float(row["totalReads"]))
    return counts, total


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fastp_json", required=True)
    parser.add_argument("--host_counts", required=True, help="combined BamFilterByChroms totalReadCounts")
    parser.add_argument("--kmer_filter_counts", required=True, help="extractionCounts of the k-mer filter pass")
    parser.add_argument("--kmer_sort_counts", required=True, help="extractionCounts of the k-mer sort pass")
    parser.add_argument("--initial_salmon_meta", required=True)
    parser.add_argument("--initial_salmon_lib_format", required=True)
    parser.add_argument("--filtered_salmon_meta", required=True)
    parser.add_argument("--pf_sets", required=True, help="filter-pass sets that are Pf reads, comma separated")
    parser.add_argument("--sort_exclude_sets", required=True, help="sort-pass sets that are not genes, comma separated")
    parser.add_argument("--gene_sets", required=True, help="sort-pass sets to report individually, comma separated")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out = {"sample": args.sample}

    # ---- fastp ----
    fastp = read_json(args.fastp_json)
    before = fastp["summary"]["before_filtering"]
    after = fastp["summary"]["after_filtering"]
    filtering = fastp.get("filtering_result", {})
    input_reads = int(before["total_reads"])
    passed = int(after["total_reads"])
    out["input_reads"] = input_reads
    out["fastp_passed_reads"] = passed
    out["fastp_removed_reads"] = input_reads - passed
    out["fastp_removed_pct"] = pct(input_reads - passed, input_reads)
    for reason in FASTP_REMOVED_REASONS:
        out[f"fastp_{reason}"] = int(filtering.get(reason, 0))
    out["q30_rate_before"] = before["q30_rate"]
    out["q30_rate_after"] = after["q30_rate"]
    out["gc_content_after"] = after["gc_content"]
    out["mean_read_length_after"] = (after["read1_mean_length"] + after.get("read2_mean_length", after["read1_mean_length"])) / 2
    out["adapter_trimmed_reads"] = int(fastp.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0))
    out["insert_size_peak"] = fastp.get("insert_size", {}).get("peak", "NA")

    # ---- host depletion ----
    filtered, non_host = host_counts(read_tsv(args.host_counts))
    host = sum(filtered.values())
    out["host1_filtered_reads"] = filtered.get("host1", 0)
    out["host2_filtered_reads"] = filtered.get("host2", 0)
    out["host_reads"] = host
    out["host_pct"] = pct(host, input_reads)
    out["non_host_reads"] = non_host
    out["non_host_pct"] = pct(non_host, input_reads)
    if host + non_host != passed:
        print(f"warning: {args.sample}: host + non-host reads ({host + non_host}) != fastp passed reads ({passed})",
              file=sys.stderr)

    # ---- k-mer filter pass: what the non-host reads are ----
    filter_counts, filter_total = kmer_counts(read_tsv(args.kmer_filter_counts))
    filter_sets = sorted(t for t in filter_counts if t not in KMER_SPECIAL_ROWS)
    accounted = 0
    for set_name in filter_sets:
        out[f"filter_{set_name}_reads"] = filter_counts[set_name]
        accounted += filter_counts[set_name]
    out["filter_undetermined_reads"] = filter_counts.get("undetermined", 0)
    accounted += filter_counts.get("undetermined", 0)
    # multi-hit and quality-filtered reads, which are dropped
    out["filter_other_reads"] = filter_total - accounted
    pf_reads = 0
    for set_name in comma_list(args.pf_sets):
        if set_name not in filter_counts:
            print(f"warning: {args.sample}: pf set {set_name} not in the k-mer filter counts", file=sys.stderr)
        pf_reads += filter_counts.get(set_name, 0)
    out["pf_reads"] = pf_reads
    out["pf_pct"] = pct(pf_reads, input_reads)
    if filter_total != non_host:
        print(f"warning: {args.sample}: k-mer filter input ({filter_total}) != non-host reads ({non_host})",
              file=sys.stderr)

    # ---- k-mer sort pass: reads per gene ----
    sort_counts, sort_total = kmer_counts(read_tsv(args.kmer_sort_counts))
    excluded = set(comma_list(args.sort_exclude_sets)) | KMER_SPECIAL_ROWS
    out["sort_input_reads"] = sort_total
    out["sorted_gene_reads"] = sum(c for t, c in sort_counts.items() if t not in excluded)
    out["sorted_gene_pct_of_pf"] = pct(out["sorted_gene_reads"], pf_reads)
    for gene in comma_list(args.gene_sets):
        if gene not in sort_counts:
            print(f"warning: {args.sample}: gene set {gene} not in the k-mer sort counts", file=sys.stderr)
        out[f"{gene}_reads"] = sort_counts.get(gene, 0)
    out["sort_undetermined_reads"] = sort_counts.get("undetermined", 0)

    # ---- salmon ----
    initial_meta = read_json(args.initial_salmon_meta)
    lib_format = read_json(args.initial_salmon_lib_format)
    filtered_meta = read_json(args.filtered_salmon_meta)
    out["salmon_initial_fragments"] = initial_meta["num_processed"]
    out["salmon_initial_percent_mapped"] = f"{initial_meta['percent_mapped']:.4f}"
    out["salmon_filtered_fragments"] = filtered_meta["num_processed"]
    out["salmon_filtered_percent_mapped"] = f"{filtered_meta['percent_mapped']:.4f}"
    out["salmon_library_type"] = lib_format.get("expected_format", "NA")
    out["salmon_compatible_fragment_ratio"] = lib_format.get("compatible_fragment_ratio", "NA")
    out["salmon_strand_mapping_bias"] = lib_format.get("strand_mapping_bias", "NA")

    with open(args.out, "w") as fh:
        fh.write("\t".join(out.keys()) + "\n")
        fh.write("\t".join(str(v) for v in out.values()) + "\n")


if __name__ == "__main__":
    main()
