process WGS_HOST_DEPLETION_SAMPLE_SUMMARY {

    tag "${sample_id}_depletion_summary"
    label 'process_single'

    input:
    tuple val(sample_id),
          path(total_read_counts_tsv_gz),
          path(pairing_stats)

    output:
    path("${sample_id}_host_depletion_summary.tab.txt")

    script:
    """
    # totalReadCounts columns: bam  condition  count  frac  total  filter
    # counts are in reads, and the 'total' column is only populated on the
    # non-zero rows, so take the max per filter stage
    input_reads=\$(gzip -dc "${total_read_counts_tsv_gz}" \\
        | awk -F'\\t' 'NR>1 && \$6=="host1_minimap2_filt" && \$5+0>m {m=\$5+0} END{print m+0}')

    host1_filtered=\$(gzip -dc "${total_read_counts_tsv_gz}" \\
        | awk -F'\\t' 'NR>1 && \$6 ~ /^host1_/ && \$2 ~ /^filtered/ {s+=\$3} END{print s+0}')

    host2_filtered=\$(gzip -dc "${total_read_counts_tsv_gz}" \\
        | awk -F'\\t' 'NR>1 && \$6 ~ /^host2_/ && \$2 ~ /^filtered/ {s+=\$3} END{print s+0}')

    host_filtered=\$(( host1_filtered + host2_filtered ))

    # pairing stats columns: sample_id kept_r1_reads kept_r2_reads final_pairs final_reads singletons_discarded
    kept_reads=\$(awk -F'\\t' 'NR==2{print \$2+\$3}' "${pairing_stats}")
    final_pairs=\$(awk -F'\\t' 'NR==2{print \$4+0}' "${pairing_stats}")
    final_reads=\$(awk -F'\\t' 'NR==2{print \$5+0}' "${pairing_stats}")
    singletons=\$(awk -F'\\t' 'NR==2{print \$6+0}' "${pairing_stats}")

    pct () {
        awk -v a="\$1" -v b="\$2" 'BEGIN{ if(b+0>0) printf "%.4f", (100.0*a)/b; else printf "NA" }'
    }

    host1_pct=\$(pct "\$host1_filtered" "\$input_reads")
    host2_pct=\$(pct "\$host2_filtered" "\$input_reads")
    host_pct=\$(pct "\$host_filtered" "\$input_reads")
    final_pct=\$(pct "\$final_reads" "\$input_reads")

    {
        printf "sample_id\\tinput_reads\\thost1_filtered_reads\\thost1_filtered_pct\\thost2_filtered_reads\\thost2_filtered_pct\\thost_filtered_reads\\thost_filtered_pct\\tkept_reads\\tsingletons_discarded\\tfinal_pairs\\tfinal_reads\\tfinal_pct_of_input\\n"
        printf "%s\\t%d\\t%d\\t%s\\t%d\\t%s\\t%d\\t%s\\t%d\\t%d\\t%d\\t%d\\t%s\\n" \\
            "${sample_id}" "\$input_reads" \\
            "\$host1_filtered" "\$host1_pct" \\
            "\$host2_filtered" "\$host2_pct" \\
            "\$host_filtered" "\$host_pct" \\
            "\$kept_reads" "\$singletons" "\$final_pairs" "\$final_reads" "\$final_pct"
    } > ${sample_id}_host_depletion_summary.tab.txt
    """
}


process WGS_HOST_DEPLETION_RUN_SUMMARY {

    tag "host_depletion_run_summary"
    label 'process_single'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*.tsv.gz"

    input:
    path(sample_summaries)
    val(pub_dir)

    output:
    path("host_depletion_summary.tsv.gz")

    script:
    """
    elucidator rBind \\
        --contains _host_depletion_summary.tab.txt \\
        --delim tab \\
        --header \\
        --out host_depletion_summary.tsv.gz
    """
}
