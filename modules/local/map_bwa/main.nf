process MAP_BWA {

    tag "${sample_id}_align_bwa"
    label 'process_medium_low_memory'

    input:
    tuple val(sample_id), path(genome_fasta_fnp), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """

    bwa mem -M  -t ${task.cpus} \\
    -R "@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:illumina\\tSM:${sample_id}\\tPU:${sample_id}" \\
      "${genome_fasta_fnp}" "${r1}" "${r2}" \\
      2> ${sample_id}.bwa.log.txt \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam \\
      2> ${sample_id}.samtools.log.txt

    samtools index ${sample_id}.sorted.bam
    """
}
