process GEN_TARGET_INFO_FROM_GENOMES_NANOPORE {
    // label 'process_medium'
    cpus   { ncpus }


    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "forSeekDeep"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "locationsByGenome"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "allExtractionCounts.tab.txt"


    input:
    path primers_fnp
    path pub_dir
    path genome_dir
    path gff_dir
    val errors_allowed
    val ncpus

    output:
    path "forSeekDeep", emit: for_seek_deep_info
    path "locationsByGenome", emit: locations_by_genome
    path "allExtractionCounts.tab.txt", emit: all_extraction_counts

    script:
    """
    SeekDeep genTargetInfoFromGenomes \
            --primers ${primers_fnp} \
            --longRangeAmplicon \
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

process GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA {
    // label 'process_medium'
    cpus   { ncpus }


    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "forSeekDeep"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "locationsByGenome"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "allExtractionCounts.tab.txt"


    input:
    path primers_fnp
    path pub_dir
    path genome_dir
    path gff_dir
    val paired_end_length
    val errors_allowed
    val ncpus

    output:
    path "forSeekDeep", emit: for_seek_deep_info
    path "locationsByGenome", emit: locations_by_genome
    path "allExtractionCounts.tab.txt", emit: all_extraction_counts

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
