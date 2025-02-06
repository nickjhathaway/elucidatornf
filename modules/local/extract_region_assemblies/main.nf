process EXTRACT_REGION_ASSEMBLIES {
    // label 'process_medium'
    cpus   { params.pw_ncpus }

    
    // publishDir "${params.pw_results_dir}", mode: 'copy', overwrite: true, pattern: "${sample}_${dirstub}"

    input:
    tuple val(sample), path(bam_fnp), path(bam_bai_fnp), path(bedfile_fnp), path(genome_dir_fnp), val(primary_genome), val(dirstub)
    
    output:
    path "${sample}_${dirstub}", emit: sample_results


    script:
    """
    PathWeaver BamExtractPathwaysFromRegion \
        --bamExtractTrimToRegion \
        --bed ${bedfile_fnp} \
        --genomeDir ${genome_dir_fnp} \
        --primaryGenome ${primary_genome} \
        --bam ${bam_fnp} \
        --dout ${sample}_${dirstub} \
        --overWriteDir \
        --numThreads ${params.pw_ncpus};
    """
}
