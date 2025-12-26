process EXTRACT_REGION_ASSEMBLIES {
    tag "${sample} extract region assemble"

    publishDir "${pub_dir}", mode: 'copy', overwrite: true

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
    def pub_outdir = "${pub_dir}/${outdir}"

    """
    set -euo pipefail

    outdir="${outdir}"
    pub_outdir="${pub_outdir}"

    if [[ -d "\$pub_outdir" ]]; then
        echo "[EXTRACT_REGION_ASSEMBLIES] Found existing results: \$pub_outdir"
        rm -rf "\$outdir"
        cp -a "\$pub_outdir" "\$outdir"
        exit 0
    fi

    PathWeaver BamExtractPathwaysFromRegion \\
        --bamExtractTrimToRegion \\
        --bed ${bedfile_fnp} \\
        --genomeDir ${genome_dir_fnp} \\
        --primaryGenome ${primary_genome} \\
        --bam ${bam_fnp} \\
        --dout "\$outdir" \\
        --overWriteDir \\
        --numThreads ${params.pw_ncpus} \\
        ${extra_args}
    rm ${outdir}/extractionLog.json
    rm ${outdir}/final/coiCounts.tab.txt
    """
}
