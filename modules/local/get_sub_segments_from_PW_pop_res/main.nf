process GET_SUB_SEGMENTS_FROM_FASTA {
    label 'process_low'

    def output_dir = "subSegments"
    // publishDir "${params.pw_results_dir}", mode: 'copy', overwrite: true, pattern: "${sample}_${dirstub}"

    input:
    tuple path(per_sample_seqs), path (bedfile_fnp), path(genome_dir_fnp), val(primary_genome), val (region_id), val (correction_occurence_cut_off), val (low_freq_cut_off)

    output:
    path "${output_dir}", emit: sub_segments_res_dir
    path "${region_id}_0_ref_variable_expanded_genomic.bed", emit: ref_variable_expanded_genomic_0_bed
    path "${region_id}_0_ref_sharedLocs_genomic.bed", emit: ref_sharedLocs_genomic_0_bed


    script:
    def genome_fnp = "${genome_dir_fnp}/${primary_genome}.fasta"
    """
    elucidator createSharedSubSegmentsFromRefSeqs \
        --genome ${genome_fnp} \
        --refBedFnp ${bedfile_fnp} \
        --refBedID ${region_id} \
        --fastagz ${per_sample_seqs} \
        --dout ${output_dir} \
        --overWriteDir \
        --correctionOccurenceCutOff ${correction_occurence_cut_off} \
        --lowFreqCutOff ${low_freq_cut_off};
        
    # rename files so they have the target_id in them 
    cat "${output_dir}/subRegionInfo/0_ref_variable_expanded_genomic.bed" > "${region_id}_0_ref_variable_expanded_genomic.bed"
    cat "${output_dir}/subRegionInfo/0_ref_sharedLocs_genomic.bed" > "${region_id}_0_ref_sharedLocs_genomic.bed"
    """
}
