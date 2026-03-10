process CONCATENATE_AMPLICON_POPULATION_CLUSTERING {
    label 'process_single'

    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "allSelectedClustersInfo.tsv.gz"
    publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "allele_data.tsv.gz"

    input:
    path clustered_fnps
    val pub_dir


    output:
    path "allSelectedClustersInfo.tsv.gz", emit: all_selected_clusters_info
    path "allele_data.tsv.gz", emit: allele_data

    script:
    """
    elucidator rBind --recursive --contains selectedClustersInfo.tab.txt.gz --delim tab --header --out allSelectedClustersInfo.tsv.gz
    elucidator tableExtractColumns --file allSelectedClustersInfo.tsv.gz --delim tab --header --columns s_Sample,p_name,h_Consensus,c_ReadCnt | \
            elucidator tableRenameColumns --rename_key "s_Sample:library_sample_name;p_name:target_name;h_Consensus:seq;c_ReadCnt:reads" --file STDIN --delim tab --header --out allele_data.tsv.gz
    """
}
