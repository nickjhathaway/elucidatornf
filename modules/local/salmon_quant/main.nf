process SALMON_QUANT {

    publishDir "${pubdir}", mode: 'copy', overwrite: true

    tag "${sample_id}_salmon_quant"
    label 'process_medium_low_memory'

    input:
    tuple val(sample_id), path(salmon_index), path(r1), path(r2), path(pubdir)

    output:
    tuple val(sample_id),
            path("${sample_id}_salmon_quants_quant.sf.gz"),
            path("${sample_id}_salmon_quants_lib_format_counts.json"),
            path("${sample_id}_salmon_quants_meta_info.json")
    script:
    """

    salmon quant \
        -i salmon_index  \
        -l A \
        -1 ${r1} -2 ${r2} \
        -p ${task.cpus} \
        --gcBias \
        -o quants_salmon
    pigz quants_salmon/quant.sf -c > ${sample_id}_quant.sf.gz
    ln -s quants_salmon/lib_format_counts.json ${sample_id}_salmon_lib_format_counts.json
    ln -s aux_info/meta_info.json ${sample_id}_salmon_meta_info.json
    """
}
