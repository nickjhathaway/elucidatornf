process MAP_MINIMAP2 {

    tag "${sample_id}_align_minimap2"
    label 'process_medium_low_memory'

    input:
    tuple val(sample_id), path(genome_mmi_fnp), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """

    minimap2 -x sr -a --secondary=no -t ${task.cpus} \\
      "${genome_mmi_fnp}" "${r1}" "${r2}" \\
      2> ${sample_id}.minimap2.log.txt \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam \\
      2> ${sample_id}.samtools.log.txt

    samtools index ${sample_id}.sorted.bam
    """
}
