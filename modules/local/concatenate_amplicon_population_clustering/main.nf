process CONCATENATE_AMPLICON_POPULATION_CLUSTERING {
    label 'process_single'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "allSelectedClustersInfo.tsv.gz"
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "allele_data.tsv.gz"
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "allSampleClustersInfo.tsv.gz"

    input:
    path clustered_fnps
    path sample_clusters_fnps
    val pub_dir


    output:
    path "allSelectedClustersInfo.tsv.gz", emit: all_selected_clusters_info
    path "allele_data.tsv.gz", emit: allele_data
    path "allSampleClustersInfo.tsv.gz", emit: all_sample_clusters_info

    script:
    // there are per replicate columns (R1_*, R2_*, ...) up to the most replicates a sample has for a target,
    // so targets can differ in their columns when samples have replicates, --fill adds the missing ones as NA
    """
    elucidator rBind --recursive --contains selectedClustersInfo.tab.txt.gz --delim tab --header --fill --out allSelectedClustersInfo.tsv.gz
    elucidator tableExtractColumns --file allSelectedClustersInfo.tsv.gz --delim tab --header --columns s_Sample,p_name,h_Consensus,c_ReadCnt | \
            elucidator tableRenameColumns --rename_key "s_Sample:library_sample_name;p_name:target_name;h_Consensus:seq;c_ReadCnt:reads" --file STDIN --delim tab --header --out allele_data.tsv.gz
    # every sample cluster, final and excluded, with its status, exclusion reasons, sequence and quality
    elucidator rBind --recursive --contains sampleClustersInfo.tab.txt.gz --delim tab --header --fill --out allSampleClustersInfo.tsv.gz
    """
}
