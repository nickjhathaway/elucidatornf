// makes the same reports as VARIANT_CALL_ON_HAP_TABLE from the variant calls population clustering made per target
process COMBINE_VARIANT_CALLS {
    label 'process_medium'

    publishDir { "${pub_results_dir}" }, mode: 'copy', overwrite: true


    input:
    path bedfile_fnp
    path input_results
    path genome_dir_fnp
    val primary_genome
    path gff_fnp
    path known_amino_acid_changes_fnp, name: "known_amino_acid_changes.tsv"
    val variant_frequency_cut_off
    val variant_occurrence_cut_off
    path meta_fnp, name: "meta.tsv"
    // the per target variant calls, staged under perTargetVariantCalls/, a target's calls are in perTargetVariantCalls/<per_target_root>/<target>/<per_target_sub_dir>
    path per_target_variant_calls, stageAs: "perTargetVariantCalls/*"
    val per_target_root
    val per_target_sub_dir
    val pub_results_dir
    val extra_args

    output:
    path "reports", emit: reports


    script:
    def extra_extra_args = task.ext.args ? task.ext.args : ''
    def per_target_dir = "" == per_target_root ? "perTargetVariantCalls" : "perTargetVariantCalls/${per_target_root}"
    def per_target_sub_dir_arg = "" == per_target_sub_dir ? "" : "--perTargetVariantCallsSubDir ${per_target_sub_dir}"
    """
    meta_arg="--metaFnp ${meta_fnp}"
    known_amino_acid_changes_arg="--knownAminoAcidChangesFnp ${known_amino_acid_changes_fnp}"
    if [ ! -s "\$(readlink -f ${meta_fnp})" ]; then meta_arg=""; fi
    if [ ! -s "\$(readlink -f ${known_amino_acid_changes_fnp})" ]; then known_amino_acid_changes_arg=""; fi

    # only the regions supplied by the bed file were variant called
    cut -f4 ${bedfile_fnp} > regions.txt
    elucidator tableExtractElementsWithLevels --levels regions.txt \
        --column p_name --file ${input_results} --delim tab \
        --header --out extracted_${input_results} --overWrite

    SeekDeep variantCallOnSeqAndProtein \
        --genomicLocations ${bedfile_fnp} \
        --resultsFnp extracted_${input_results} \
        --genome ${genome_dir_fnp}/${primary_genome}.fasta \
        --gff ${gff_fnp} \
        \${known_amino_acid_changes_arg} \
        --variantFrequencyCutOff ${variant_frequency_cut_off} \
        --variantOccurrenceCutOff ${variant_occurrence_cut_off} \
        --ignoreSubFields site:LabCross,site:LabControl,site:LabContaminated \
        --dout variantCalls \
        --numThreads ${params.variant_calling_ncpus} \
        \${meta_arg} \
        --exportLabIsolateSeqs \
        --perTargetVariantCallsDir ${per_target_dir} ${per_target_sub_dir_arg} \
        ${extra_args}\
        ${extra_extra_args}

    ln -s variantCalls/reports
    """
}
