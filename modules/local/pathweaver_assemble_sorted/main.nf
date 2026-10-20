// PathWeaver assembly of one gene/gene family's k-mer-sorted reads for one sample.
// In the sorted fastqs R2 is on the opposite strand to R1, hence --revCompMate.
// Only output_aboveCutOff.fasta is kept; it is written empty when nothing assembles,
// since PathWeaver then writes no fasta at all.

process PATHWEAVER_ASSEMBLE_SORTED {

    tag "${sample_id} ${gene} PathWeaver"
    label 'process_medium_low_memory'
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*_output_aboveCutOff.fasta"

    input:
    // singles: the set's unpaired reads, [] when there are none
    tuple val(sample_id), val(gene), path(r1), path(r2), path(singles), val(pub_dir)

    output:
    tuple val(sample_id), val(gene), path("${gene}_output_aboveCutOff.fasta"), emit: assemblies

    script:
    def singles_arg = singles ? "--fastq ${singles}" : ''
    """
    PathWeaver SeqsExtractPathways \\
        ${singles_arg} \\
        --fastq1 ${r1} \\
        --fastq2 ${r2} \\
        --revCompMate \\
        --dout ${gene}_PathWeaver \\
        --overWriteDir \\
        --numThreads ${task.cpus} \\
        --sampName ${sample_id} \\
        ${params.rnaseq_pathweaver_args}

    if [ -f ${gene}_PathWeaver/${sample_id}/output_aboveCutOff.fasta ]; then
        cp ${gene}_PathWeaver/${sample_id}/output_aboveCutOff.fasta ${gene}_output_aboveCutOff.fasta
    else
        touch ${gene}_output_aboveCutOff.fasta
    fi
    """

    stub:
    """
    touch ${gene}_output_aboveCutOff.fasta
    """
}
