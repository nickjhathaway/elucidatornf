process GEN_TARGET_INFO_FROM_GENOMES_NANOPORE {
    // label 'process_medium'
    cpus   { ncpus }


    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "genome_extraction/forSeekDeep"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "genome_extraction/locationsByGenome"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "genome_extraction/allExtractionCounts.tab.txt"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "targets_with_extractins.txt"


    input:
    path primers_fnp
    val pub_dir
    path genome_dir
    path gff_dir
    val errors_allowed
    val ncpus

    output:
    path "genome_extraction/forSeekDeep", emit: for_seek_deep_info
    path "genome_extraction/locationsByGenome", emit: locations_by_genome
    path "genome_extraction/allExtractionCounts.tab.txt", emit: all_extraction_counts
    path "targets_with_extractins.txt", emit: targets_with_extractins

    script:
    if (file("${pub_dir}/genome_extraction/forSeekDeep").exists() && file("${primers_fnp}").lastModified() < file("${pub_dir}/genome_extraction/forSeekDeep").lastModified()){
        """
        mkdir genome_extraction
        cp -r "${pub_dir}/genome_extraction/forSeekDeep" genome_extraction/
        cp -r "${pub_dir}/genome_extraction/locationsByGenome" genome_extraction/
        cp "${pub_dir}/genome_extraction/allExtractionCounts.tab.txt" genome_extraction/
        elucidator tableExtractCriteria --file genome_extraction/allExtractionCounts.tab.txt --delim tab --columnName extractionCounts --header --cutOff 0 | elucidator printCol --file STDIN --delim tab --header --columnName target --sort --unique  > targets_with_extractins.txt
        """
    } else {
        """
        SeekDeep genTargetInfoFromGenomes \
                --primers ${primers_fnp} \
                --longRangeAmplicon \
                --genomeDir ${genome_dir} \
                --gffDir ${gff_dir} \
                --dout genome_extraction \
                --errors ${errors_allowed} \
                --numThreads ${ncpus} \
                --useBlast
        elucidator tableExtractCriteria --file genome_extraction/allExtractionCounts.tab.txt --delim tab --columnName extractionCounts --header --cutOff 0 | elucidator printCol --file STDIN --delim tab --header --columnName target --sort --unique  > targets_with_extractins.txt
        """
    }
}

process GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA {
    // label 'process_medium'
    cpus   { ncpus }


    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "forSeekDeep"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "locationsByGenome"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "allExtractionCounts.tab.txt"


    input:
    path primers_fnp
    val pub_dir
    path genome_dir
    path gff_dir
    val paired_end_length
    val errors_allowed
    val ncpus

    output:
    path "forSeekDeep", emit: for_seek_deep_info
    path "locationsByGenome", emit: locations_by_genome
    path "extraction/allExtractionCounts.tab.txt", emit: all_extraction_counts

    script:
    """
    SeekDeep genTargetInfoFromGenomes \
            --primers ${primers_fnp} \
            --pairedEndLength ${paired_end_length} \
            --genomeDir ${genome_dir} \
            --gffDir ${gff_dir} \
            --dout extraction \
            --errors ${errors_allowed} \
            --numThreads ${ncpus} \
            --useBlast
    ln -s extraction/forSeekDeep
    ln -s extraction/locationsByGenome
    """
}
