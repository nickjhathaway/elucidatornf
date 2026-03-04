process VARIANT_CALL_ON_HAP_TABLE {
    label 'process_medium'

    publishDir "${pub_results_dir}", mode: 'copy', overwrite: true


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
    val getting_pairwise_comps
    val meta_fields_to_calc_pop_diffs
    val pub_results_dir
    val extra_args

    output:
    path "reports", emit: reports
    // path "variantCalls/runLog*", emit: run_log


    script:
    // def meta_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--metaFnp ${meta_fnp}"
    // def known_amino_acid_changes_arg = "EMPTY_FILE.txt" == "${known_amino_acid_changes_fnp}" ? "" : "--knownAminoAcidChangesFnp ${known_amino_acid_changes_fnp}"
    def getting_pairwise_comps_arg = getting_pairwise_comps ? "--getPairwiseComps" : ""
    def meta_fields_to_calc_pop_diffs_arg = "" == meta_fields_to_calc_pop_diffs ? "" : "--metaFieldsToCalcPopDiffs ${meta_fields_to_calc_pop_diffs}"
    """
    meta_arg="--metaFnp ${meta_fnp}"
    known_amino_acid_changes_arg="--knownAminoAcidChangesFnp ${known_amino_acid_changes_fnp}"
    if [ ! -s "\$(readlink -f ${meta_fnp})" ]; then meta_arg=""; fi
    if [ ! -s "\$(readlink -f ${known_amino_acid_changes_fnp})" ]; then known_amino_acid_changes_arg=""; fi

    # only variant call on the regions supplied by the bed file
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
        ${getting_pairwise_comps_arg} \
        ${meta_fields_to_calc_pop_diffs_arg}\
        ${extra_args}

    ln -s variantCalls/reports
    """
}


