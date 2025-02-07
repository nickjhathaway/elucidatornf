process VARIANT_CALL_ON_HAP_TABLE {
    label 'process_medium'

    publishDir "${pub_results_dir}", mode: 'copy', overwrite: true


    input:
    path bedfile_fnp
    path input_results
    path genome_dir_fnp
    val primary_genome
    path gff_fnp
    path known_amino_acid_changes_fnp
    val variant_frequency_cut_off
    val variant_occurrence_cut_off
    path meta_fnp
    val getting_pairwise_comps
    val meta_fields_to_calc_pop_diffs
    val pub_results_dir
    val extra_args

    output:
    path "reports", emit: reports


    script:
    def meta_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--metaFnp ${meta_fnp}"
    def known_amino_acid_changes_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--knownAminoAcidChangesFnp ${known_amino_acid_changes_fnp}"
    def getting_pairwise_comps_arg = getting_pairwise_comps ? "" : "--getPairwiseComps"
    def meta_fields_to_calc_pop_diffs_arg = "" == meta_fields_to_calc_pop_diffs ? "" : "--metaFieldsToCalcPopDiffs ${meta_fields_to_calc_pop_diffs}"
    """

    SeekDeep variantCallOnSeqAndProtein \
        --genomicLocations ${bedfile_fnp} \
        --resultsFnp ${input_results} \
        --genome ${genome_dir_fnp}/${primary_genome}.fasta \
        --gff ${gff_fnp} \
        ${known_amino_acid_changes_arg} \
        --variantFrequencyCutOff ${variant_frequency_cut_off} \
        --variantOccurrenceCutOff ${variant_occurrence_cut_off} \
        --ignoreSubFields site:LabCross,site:LabControl,site:LabContaminated \
        --dout variantCalls \
        --numThreads ${params.variant_calling_ncpus} \
        ${meta_arg} \
        --exportLabIsolateSeqs \
        ${getting_pairwise_comps_arg} \
        ${meta_fields_to_calc_pop_diffs_arg}\
        ${extra_args}

    ln -s variantCalls/reports
    """
}
//SeekDeep variantCallOnSeqAndProtein --overWriteDir --genomicLocations /tank/data/plasmodium/falciparum/beds/targeted_panels/MADHITO/MADHITO_Pf3D7_withPrimers.bed --resultsFnp allSelectedClustersInfo.tab.txt.gz --genome /tank/data/genomes/plasmodium/genomes/pf/genomes/Pf3D7.fasta --gff /tank/data/genomes/plasmodium/genomes/pf/info/gff/Pf3D7.gff --knownAminoAcidChangesFnp /tank/data/genomes/plasmodium/genomes/pf/info/pf_drug_resistant_aaPositions.tsv --variantFrequencyCutOff 0.0099 --variantOccurrenceCutOff 10 --ignoreSubFields site:LabCross,site:LabControl,site:LabContaminated --dout variantCalls --numThreads 12 --metaFnp /tank/data/plasmodium/falciparum/pfdata/metadata/metaByBioSample.tab.txt --exportLabIsolateSeqs --selectTargets Pf3D7_03_v3_regions.txt --getPairwiseComps --metaFieldsToCalcPopDiffs region,country,subRegion,secondaryRegion


