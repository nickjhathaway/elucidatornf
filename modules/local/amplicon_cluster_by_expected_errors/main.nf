process AMPLICON_CLUSTER_BY_EXPECTED_ERRORS {
    tag "${replicate_name}_${target_name}_cluster_by_expected_errors"

    label 'process_low'

    input:
    tuple path(fastq_fnp), val(sample_name), val(replicate_name), val(target_name), val(error_profile)


    output:
    tuple val(sample_name), val(target_name), path("${replicate_name}_${target_name}_output.fastq.gz"), emit: output_results

    script:
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    SeekDeep clusterDown \
            --no_run_log \
            --${error_profile} \
            --fastqgz ${fastq_fnp} \
            --sample ${sample_name} \
            --replicate ${replicate_name} \
            --target ${target_name} \
            --out ${replicate_name}_${target_name}_output.fastq.gz \
            ${extra_args}
    """
}
