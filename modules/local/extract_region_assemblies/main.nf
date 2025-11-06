process EXTRACT_REGION_ASSEMBLIES {
    // label 'process_medium'
    cpus   { params.pw_ncpus }

    input:
    tuple val(sample), path(bam_fnp), path(bam_bai_fnp), path(bedfile_fnp), path(genome_dir_fnp), val(primary_genome), val(dirstub)

    output:
    path "${sample}_${dirstub}", emit: sample_results


    script:
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    PathWeaver BamExtractPathwaysFromRegion \
        --bamExtractTrimToRegion \
        --bed ${bedfile_fnp} \
        --genomeDir ${genome_dir_fnp} \
        --primaryGenome ${primary_genome} \
        --bam ${bam_fnp} \
        --dout ${sample}_${dirstub} \
        --overWriteDir \
        --numThreads ${params.pw_ncpus} \
        $extra_args
    """
}
