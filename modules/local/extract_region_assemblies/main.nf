process EXTRACT_REGION_ASSEMBLIES {
    tag "${sample} extract region assemble"

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true

    cpus { params.pw_ncpus }

    input:
    tuple val(sample),
          path(bam_fnp),
          path(bam_bai_fnp),
          path(bedfile_fnp),
          path(genome_dir_fnp),
          val(primary_genome),
          val(dirstub),
          val(pub_dir)

    output:
    path "${sample}_${dirstub}", emit: sample_results

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def outdir     = "${sample}_${dirstub}"

    // Note: samples whose results already exist are filtered out in the calling
    // workflow (the existing output is channelled directly), so this process only
    // runs for samples that actually need computing -- no copy-in / copy-out.
    """
    set -euo pipefail

    PathWeaver BamExtractPathwaysFromRegion \\
        --bamExtractTrimToRegion \\
        --bed ${bedfile_fnp} \\
        --genomeDir ${genome_dir_fnp} \\
        --primaryGenome ${primary_genome} \\
        --bam ${bam_fnp} \\
        --dout "${outdir}" \\
        --overWriteDir \\
        --numThreads ${params.pw_ncpus} \\
        ${extra_args}
    rm -f ${outdir}/extractionLog.json
    rm -f ${outdir}/final/coiCounts.tab.txt
    """
}
