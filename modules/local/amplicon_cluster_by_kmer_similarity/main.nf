process AMPLICON_CLUSTER_BY_KMER_SIMILARITY {
    label 'process_low'

    cpus   { ncpus }

    input:
    tuple path(fastq_fnp), val(sample_name), val(target_name), val(ncpus)


    output:
    tuple val(sample_name), val(target_name), path("${sample_name}_${target_name}_output.fastq.gz"), emit: output_results

    script:
    """
    SeekDeep kluster \
            --fastqgz ${fastq_fnp} \
            --sample ${sample_name} \
            --target ${target_name} \
            --dout klusterOut \
            --numThreads ${ncpus}
    ln -s klusterOut/output.fastq.gz ${sample_name}_${target_name}_output.fastq.gz
    """
}
