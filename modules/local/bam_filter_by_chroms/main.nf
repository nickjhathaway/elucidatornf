process BAM_FILTER_BY_CHROMS {

    tag "${sample_id} ${filter_name} filter bam"

    label 'process_single'

    input:
    tuple val(sample_id), val(filter_name), path(bam), path(bam_bai), path(host_chroms), val(write_unmapped_sep)

    output:
    tuple val(sample_id),
          path("${sample_id}_${filter_name}_kept_R1.fastq.gz"),
          path("${sample_id}_${filter_name}_kept_R2.fastq.gz"),
          path("${sample_id}_${filter_name}_unmapped_R1.fastq.gz"),
          path("${sample_id}_${filter_name}_unmapped_R2.fastq.gz"),
          path("${sample_id}_${filter_name}_filteredByChrom.tab.txt"),
          path("${sample_id}_${filter_name}_totalReadCounts.tab.txt")

    script:

    def any_mate_flag                  = params.wgs_host_filter_any_mate                  ? '--any_mate' : ''
    def filter_with_unmapped_mate_flag = params.wgs_host_filter_filter_with_unmapped_mate ? '--filterWithUnmappedMate' : ''
    def unmapped_flag                  = write_unmapped_sep                           ? '--writeOutUnmappedSeparately' : ''


    """
    elucidator BamFilterByChroms \\
        --chroms "${host_chroms}" \\
        ${any_mate_flag} \\
        --minMappingQuality ${params.wgs_host_filter_min_mapq} \\
        ${filter_with_unmapped_mate_flag} \\
        --doNotWriteFilterOff \\
        ${unmapped_flag} \\
        --bam "${bam}" \\
        --overWrite \\
        --out "${sample_id}_${filter_name}"

    # BamFilterByChroms generates:
    #   \${sample_id}_${filter_name}_kept_R1.fastq.gz
    #   \${sample_id}_${filter_name}_kept_R2.fastq.gz
    #   \${sample_id}_${filter_name}_filteredOff_R1.fastq.gz (if not --doNotWriteFilterOff)
    #   \${sample_id}_${filter_name}_unmapped_R1.fastq.gz (if --writeOutUnmappedSeparately)
    #   \${sample_id}_${filter_name}_filteredByChrom.tab.txt
    #   \${sample_id}_${filter_name}_totalReadCounts.tab.txt

    # Ensure all expected fastq outputs exist, even if empty, this will make workflow flow
    touch_empty_gz () {
        [ -f "\$1" ] || gzip -c </dev/null > "\$1"
    }

    mv ${sample_id}_${filter_name}_filteredByChrom.tab.txt original_${sample_id}_${filter_name}_filteredByChrom.tab.txt
    elucidator addColumn --file original_${sample_id}_${filter_name}_filteredByChrom.tab.txt --newColumnName filter --header --delim tab --element ${filter_name} --out ${sample_id}_${filter_name}_filteredByChrom.tab.txt
    mv ${sample_id}_${filter_name}_totalReadCounts.tab.txt original_${sample_id}_${filter_name}_totalReadCounts.tab.txt
    elucidator addColumn --file original_${sample_id}_${filter_name}_totalReadCounts.tab.txt --newColumnName filter --header --delim tab --element ${filter_name} --out ${sample_id}_${filter_name}_totalReadCounts.tab.txt

    touch_empty_gz "${sample_id}_${filter_name}_kept_R1.fastq.gz"
    touch_empty_gz "${sample_id}_${filter_name}_kept_R2.fastq.gz"
    touch_empty_gz "${sample_id}_${filter_name}_unmapped_R1.fastq.gz"
    touch_empty_gz "${sample_id}_${filter_name}_unmapped_R2.fastq.gz"
    """
}
