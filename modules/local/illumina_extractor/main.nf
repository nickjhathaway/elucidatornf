process ILLUMINA_EXTRACTOR {
    tag "${sample_name}_illumina_extractor"
    label 'process_single'


    input:
    tuple path(fastq1_fnp), path(fastq2_fnp), val(sample_name), path(primers_fnp), path(overlap_status_fnp), path(length_cut_offs_per_target), path(auto_flags_fnp)

    output:
    tuple val(sample_name), path("extraction/*.fastq.gz"), emit: fastqs_per_target, optional : true
    tuple val(sample_name), path("extraction/"), emit: sample_dir

    path "${sample_name}_extractionProfile.tab.txt", emit: extraction_profile_per_target
    path "${sample_name}_extractionStats.tab.txt", emit: extraction_stats
    path "${sample_name}_allFailedPrimerCounts.tab.txt", emit: failed_primer_counts
    path "${sample_name}_processPairsCounts.tab.txt", emit: process_pair_counts

    script:
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    SeekDeep extractorPairedEnd \
            --dout extraction \
            --overlapStatusFnp ${overlap_status_fnp} \
            --rename \
            --fastq1gz ${fastq1_fnp} \
            --fastq2gz ${fastq2_fnp} \
            --id ${primers_fnp} \
            --sampleName ${sample_name} \
            --lenCutOffs ${length_cut_offs_per_target} \
            ${extra_args} \$(cat ${auto_flags_fnp})
    ln -s extraction/extractionProfile.tab.txt ${sample_name}_extractionProfile.tab.txt
    ln -s extraction/extractionStats.tab.txt ${sample_name}_extractionStats.tab.txt
    ln -s extraction/allFailedPrimerCounts.tab.txt ${sample_name}_allFailedPrimerCounts.tab.txt
    ln -s extraction/processPairsCounts.tab.txt ${sample_name}_processPairsCounts.tab.txt

    """
}
