process WGS_BUILD_FINAL_SAMPLE_SUMMARY {

  tag "${sample_id}_final_summary"
  label 'process_single'

  // publish wherever you want; you already have fastp_trim_info_dir
  publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "*.tsv.gz"

  input:
  tuple val(sample_id),
        val(pub_dir),
        path(fastp_json),
        path(totalReadCounts_tsv_gz),
        path(flagstat_txt)

  output:
  path("${sample_id}_final_summary.tsv.gz")

  script:
  """
  set -euo pipefail

  cat_or_zcat () {
    local f="\$1"
    if [[ "\$f" == *.gz ]]; then
      zcat < "\$f"
    else
      cat "\$f"
    fi
  }

  # ---- fastp counts ----
  before_reads=\$(cat_or_zcat "${fastp_json}" | jq -r '.summary.before_filtering.total_reads')
  after_reads=\$(cat_or_zcat "${fastp_json}"  | jq -r '.summary.after_filtering.total_reads')
  fastp_lost=\$(( before_reads - after_reads ))

  # ---- host filter counts from combined totalReadCounts ----
  # columns: bam  condition  count  frac  total
  host_filtered_pairs=\$(cat_or_zcat "${totalReadCounts_tsv_gz}" | awk 'NR==1{next} \$2=="filteredPairs"{s+=\$3} END{print s+0}')
  host_kept_pairs=\$(cat_or_zcat "${totalReadCounts_tsv_gz}"     | awk 'NR==1{next} \$2=="keptPairs"{s+=\$3} END{print s+0}')

  host_filtered_pct="NA"
  host_kept_pct="NA"
  if [ "\$after_reads" -gt 0 ]; then
    host_filtered_pct=\$(awk -v a="\$host_filtered_pairs" -v b="\$after_reads" 'BEGIN{printf "%.4f", (100.0*a)/b}')
    host_kept_pct=\$(awk -v a="\$host_kept_pairs" -v b="\$after_reads" 'BEGIN{printf "%.4f", (100.0*a)/b}')
  fi

  # ---- final BAM primary mapped ----
  bam_primary_mapped=\$(awk '/ primary mapped /{print \$1+0; exit}' "${flagstat_txt}")
  bam_properly_paired=\$(awk '/ properly paired /{print \$1+0; exit}' "${flagstat_txt}")

  bam_primary_mapped_pct="NA"
  bam_properly_paired_pct="NA"
  if [ "\$host_kept_pairs" -gt 0 ]; then
    bam_primary_mapped_pct=\$(awk -v a="\$bam_primary_mapped" -v b="\$host_kept_pairs" 'BEGIN{printf "%.4f", (100.0*a)/b}')
    bam_properly_paired_pct=\$(awk -v a="\$bam_properly_paired" -v b="\$host_kept_pairs" 'BEGIN{printf "%.4f", (100.0*a)/b}')
  fi

  # ---- write summary ----
  {
    printf "sample_id\\tfastp_before_reads\\tfastp_after_reads\\tfastp_lost_reads\\thost_filtered_pairs\\thost_filtered_pct_of_fastp_after\\thost_kept_pairs\\thost_kept_pct_of_fastp_after\\tfinal_bam_primary_mapped_reads\\tfinal_bam_primary_mapped_pct\\tfinal_bam_properly_paired_mapped_reads\\tfinal_bam_properly_paired_mapped_pct\\n"
    printf "%s\\t%d\\t%d\\t%d\\t%d\\t%s\\t%d\\t%s\\t%d\\t%s\\t%d\\t%s\\n" \\
      "${sample_id}" "\$before_reads" "\$after_reads" "\$fastp_lost" \\
      "\$host_filtered_pairs" "\$host_filtered_pct" "\$host_kept_pairs" "\$host_kept_pct" \\
      "\$bam_primary_mapped" "\$bam_primary_mapped_pct" "\$bam_properly_paired" "\$bam_properly_paired_pct"
  } | gzip -c > ${sample_id}_final_summary.tsv.gz
  """
}
