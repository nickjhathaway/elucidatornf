process COMBINE_KEPT_FILTERED_FASTQS{
    tag "${sample_id} combine kept fastqs"

    label 'process_single'

    input:
    tuple val(sample_id), path(kept_host1_minimap2_r1),
                          path(kept_host1_minimap2_r2),
                          path(kept_host1_bwa_r1),
                          path(kept_host1_bwa_r2),
                          path(kept_host2_minimap2_r1),
                          path(kept_host2_minimap2_r2),
                          path(kept_host2_bwa_r1),
                          path(kept_host2_bwa_r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_kept_R1.fastq.gz"),
          path("${sample_id}_kept_R2.fastq.gz")

    script:

    """
    cat ${kept_host1_minimap2_r1} ${kept_host1_bwa_r1} ${kept_host2_minimap2_r1} ${kept_host2_bwa_r1} > ${sample_id}_kept_R1.fastq.gz
    cat ${kept_host1_minimap2_r2} ${kept_host1_bwa_r2} ${kept_host2_minimap2_r2} ${kept_host2_bwa_r2} > ${sample_id}_kept_R2.fastq.gz
    """
}
