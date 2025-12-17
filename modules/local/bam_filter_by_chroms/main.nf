process BAM_FILTER_BY_CHROMS {

    tag "${sample_id} filter bam ${host_chroms}"

    label 'process_single'

    input:
    tuple val(sample_id), path(bam), path(host_chroms), val(write_unmapped_sep)

    output:
    tuple val(sample_id),
          path("${sample_id}_kept_R1.fastq.gz"),
          path("${sample_id}_kept_R2.fastq.gz"),
          path("${sample_id}_unmapped_R1.fastq.gz"),
          path("${sample_id}_unmapped_R2.fastq.gz"),
          path("${sample_id}_filteredByChrom.tab.txt"),
          path("${sample_id}_totalReadCounts.tab.txt")

    script:

    def any_mate_flag                  = params.host_filter_any_mate                  ? '--any_mate' : ''
    def filter_with_unmapped_mate_flag = params.host_filter_filter_with_unmapped_mate ? '--filterWithUnmappedMate' : ''
    def unmapped_flag                  = write_unmapped_sep                           ? '--writeOutUnmappedSeparately' : ''


    """
    elucidator BamFilterByChroms \\
        --chroms "${host_chroms}" \\
        ${any_mate_flag} \\
        --minMappingQuality ${params.host_filter_min_mapq} \\
        ${filter_with_unmapped_mate_flag} \\
        --doNotWriteFilterOff \\
        ${unmapped_flag} \\
        --bam "${bam}" \\
        --overWrite \\
        --out "${sample_id}"

    # BamFilterByChroms generates:
    #   \${sample_id}_kept_R1.fastq.gz
    #   \${sample_id}_kept_R2.fastq.gz
    #   \${sample_id}_filteredOff_R1.fastq.gz (if not --doNotWriteFilterOff)
    #   \${sample_id}_unmapped_R1.fastq.gz (if --writeOutUnmappedSeparately)
    #   \${sample_id}_filteredByChrom.tab.txt
    #   \${sample_id}_totalReadCounts.tab.txt

    # Ensure all expected fastq outputs exist, even if empty, this will make workflow flow
    touch_empty_gz () {
        [ -f "\$1" ] || gzip -c </dev/null > "\$1"
    }

    touch_empty_gz "${sample_id}_kept_R1.fastq.gz"
    touch_empty_gz "${sample_id}_kept_R2.fastq.gz"
    touch_empty_gz "${sample_id}_unmapped_R1.fastq.gz"
    touch_empty_gz "${sample_id}_unmapped_R2.fastq.gz"
    """
}
