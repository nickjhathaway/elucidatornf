process MAP_BWA {

    tag "${sample_id}_align_bwa"
    label 'process_medium_low_memory'

    input:
    tuple val(sample_id), path(genome_fasta_fnp), path(genome_index_files), path(r1), path(r2)

    output:
    tuple val(sample_id),
            path("${sample_id}.sorted.bam"),
            path("${sample_id}.sorted.bam.bai"),
            path("${sample_id}.flagstat.txt")
    script:
    """

    bwa mem -M  -t ${task.cpus} \\
    -R "@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:illumina\\tSM:${sample_id}\\tPU:${sample_id}" \\
      "${genome_fasta_fnp}" "${r1}" "${r2}" \\
      2> ${sample_id}.bwa.log.txt \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam \\
      2> ${sample_id}.samtools.log.txt

    samtools index ${sample_id}.sorted.bam
    samtools flagstat -@ ${task.cpus} ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt

    """
}

process MAP_BWA_MARKDUP {

    tag "${sample_id}_align_bwa_markdup"
    label 'process_medium_low_memory'

    input:
    tuple val(sample_id), path(genome_fasta_fnp), path(genome_index_files), path(r1), path(r2)

    output:
    tuple val(sample_id),
            path("${sample_id}.sorted.bam"),
            path("${sample_id}.sorted.bam.bai"),
            path("${sample_id}.flagstat.txt"),
            path("${sample_id}.markdup_stats.txt")
    script:
    // the samtools stages overlap with bwa in the pipe, so they share the allocation
    def aux_cpus = Math.max(1, (task.cpus as int).intdiv(2))
    """

    bwa mem -M  -t ${task.cpus} \\
    -R "@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:illumina\\tSM:${sample_id}\\tPU:${sample_id}" \\
      "${genome_fasta_fnp}" "${r1}" "${r2}" \\
      2> ${sample_id}.bwa.log.txt \\
    | samtools fixmate -@ ${aux_cpus} -m -u - - \\
      2> ${sample_id}.samtools.fixmate.log.txt \\
    | samtools sort -@ ${aux_cpus} -u -T ${sample_id}.srt_tmp - \\
      2> ${sample_id}.samtools.sort.log.txt \\
    | samtools markdup -@ ${aux_cpus} --write-index \\
      -f ${sample_id}.markdup_stats.txt \\
      - "${sample_id}.sorted.bam##idx##${sample_id}.sorted.bam.bai" \\
      2> ${sample_id}.samtools.markdup.log.txt

    samtools flagstat -@ ${task.cpus} ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt

    """
}
