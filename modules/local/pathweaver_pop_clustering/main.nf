process PATHWEAVER_POP_CLUSTERING {
    label 'process_high'
    cpus   { params.pw_pop_clustering_ncpus }

    publishDir "${pub_results_dir}", mode: 'copy', overwrite: true, pattern : "reports"
    publishDir "${pub_results_dir}", mode: 'copy', overwrite: true, pattern : "info"

    input:
    val dirstub
    path meta_fnp
    path res_folders
    val pub_results_dir

    output:
    path "popClustering/", emit: pop_clustering_dir
    path "popClustering/reports/allSelectedClustersInfo.tab.txt.gz", emit: all_selected_clusters_info
    path "popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz", emit: all_slim_selected_clusters_info
    path "reports/", emit: pop_clustering_reports_dir
    path "info/", emit: pop_clustering_info_dir
    path "popClustering/reports/targets_with_results.txt", emit: targets_with_results

    script:
    def meta_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--groupingsFile ${meta_fnp}"
    """
    PathWeaver runProcessClustersOnRecon --overWriteDir --dout popClustering --pat _${dirstub} --numThreads ${params.pw_pop_clustering_ncpus} ${meta_arg}
    elucidator tableExtractCriteria --file popClustering/reports/seqsPerTargetGatheredDetailed.tab.txt --delim tab --columnName nSamples --header  --cutOff 1 | elucidator printCol --file STDIN --delim tab --header --columnName target > popClustering/reports/targets_with_results.txt
    elucidator tableExtractColumns --file popClustering/reports/allSelectedClustersInfo.tab.txt.gz --columns s_Sample::sample,p_name::target,h_popUID::allele_name,h_Consensus::seq,c_ReadCnt::read_count --delim tab --header --out popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz
    ln -s popClustering/reports
    ln -s popClustering/info
    """
}
