process CONCATENATE_ILLUMINA_EXTRACTOR {
    label 'process_single'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "all*.tsv.gz"

    input:
    path extraction_profile_per_target_fnps
    path extraction_stats_fnps
    path failed_primer_counts_fnps
    path process_pair_counts_fnps
    val pub_dir

    output:
    path "allExtractionProfile.tsv.gz", emit: all_extraction_profile_per_target
    path "allExtractionStats.tsv.gz", emit: all_extraction_stats
    path "allFailedPrimerCounts.tsv.gz", emit: all_failed_primer_counts
    path "allProcessPairsCounts.tsv.gz", emit: all_process_pair_counts
    path "targets_with_passing_reads.txt", emit: targets_with_passing_reads

    script:
    """
    elucidator rBind --recursive --contains extractionProfile.tab.txt --delim tab --header --out STDOUT | elucidator trimContent --overWrite --file STDIN --delim tab --header --trimAt "(" --out allExtractionProfile.tsv.gz
    elucidator rBind --recursive --contains extractionStats.tab.txt   --delim tab --header --out STDOUT | elucidator trimContent --overWrite --file STDIN --delim tab --header --trimAt "(" --out allExtractionStats.tsv.gz
    elucidator rBind --recursive --contains allFailedPrimerCounts.tab.txt --delim tab --header --out STDOUT | elucidator trimContent --overWrite --file STDIN --delim tab --header --trimAt "(" --out allFailedPrimerCounts.tsv.gz
    elucidator rBind --recursive --contains processPairsCounts.tab.txt   --delim tab --header --out STDOUT | elucidator trimContent --overWrite --file STDIN --delim tab --header --trimAt "(" --out allProcessPairsCounts.tsv.gz

    elucidator tableExtractCriteria --file allExtractionProfile.tsv.gz --delim tab --header --columnName good | elucidator printCol --file STDIN --delim tab --header --columnName target --sort --unique --out targets_with_passing_reads.txt
    """
}
