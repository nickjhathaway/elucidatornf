process AMPLICON_CLUSTER_BY_EXPECTED_ERRORS {
    tag "${sample_name}_${target_name}_cluster_by_expected_errors"

    label 'process_low'

    input:
    tuple path(fastq_fnp), val(sample_name), val(target_name), val(error_profile)


    output:
    tuple val(sample_name), val(target_name), path("${sample_name}_${target_name}_output.fastq.gz"), emit: output_results

    script:
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    SeekDeep clusterDown \
            --${error_profile} \
            --fastqgz ${fastq_fnp} \
            --sample ${sample_name} \
            --target ${target_name} \
            --dout custer_out \
            ${extra_args}
    ln -s custer_out/output.fastq.gz ${sample_name}_${target_name}_output.fastq.gz
    """
}
