process PATHWEAVER_POP_CLUSTERING {
    label 'process_high'
    cpus   { params.pw_pop_clustering_ncpus }

    publishDir { "${pub_results_dir}" }, mode: 'copy', overwrite: true, pattern : "reports"
    publishDir { "${pub_results_dir}" }, mode: 'copy', overwrite: true, pattern : "info"

    input:
    val dirstub
    path meta_fnp
    path assemblies_root
    val pub_results_dir

    output:
    path "popClustering/", emit: pop_clustering_dir
    path "popClustering/reports/allSelectedClustersInfo.tab.txt.gz", emit: all_selected_clusters_info
    path "popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz", emit: all_slim_selected_clusters_info
    path "reports/", emit: pop_clustering_reports_dir
    path "info/", emit: pop_clustering_info_dir
    path "popClustering/reports/targets_with_results.txt", emit: targets_with_results
    path "output_sampleMetaData.tab.txt", optional: true, emit: output_sample_meta_fnp

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def meta_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--groupingsFile ${meta_fnp}"
    """
    PathWeaver runProcessClustersOnRecon --inputDirectory ${assemblies_root} --overWriteDir --dout popClustering --pat _${dirstub} --numThreads ${params.pw_pop_clustering_ncpus} ${meta_arg} ${extra_args}
    elucidator tableExtractCriteria --file popClustering/reports/seqsPerTargetGatheredDetailed.tab.txt --delim tab --columnName nSamples --header  --cutOff 1 | elucidator printCol --file STDIN --delim tab --header --columnName target > popClustering/reports/targets_with_results.txt
    elucidator tableExtractColumns --file popClustering/reports/allSelectedClustersInfo.tab.txt.gz --columns s_Sample::library_sample_name,p_name::target_name,h_popUID::microhaplotype_name,h_Consensus::seq,c_ReadCnt::reads --delim tab --header --out popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz
    ln -s popClustering/reports
    ln -s popClustering/info
    if [ -f popClustering/info/sampleMetaData.tab.txt ]; then ln -s popClustering/info/sampleMetaData.tab.txt output_sampleMetaData.tab.txt; fi;

    """
}

process PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED {
    label 'process_high'
    cpus   { params.pw_pop_clustering_ncpus }

    publishDir { "${pub_results_dir}" }, mode: 'copy', overwrite: true, pattern : "reports"
    publishDir { "${pub_results_dir}" }, mode: 'copy', overwrite: true, pattern : "info"

    input:
    val dirstub
    path meta_fnp
    path assemblies_root
    val pub_results_dir
    path trim_bed
    path genome_2bit_fnp

    output:
    path "popClustering/", emit: pop_clustering_dir
    path "popClustering/reports/allSelectedClustersInfo.tab.txt.gz", emit: all_selected_clusters_info
    path "popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz", emit: all_slim_selected_clusters_info
    path "reports/", emit: pop_clustering_reports_dir
    path "info/", emit: pop_clustering_info_dir
    path "popClustering/reports/targets_with_results.txt", emit: targets_with_results
    path "output_sampleMetaData.tab.txt", optional: true, emit: output_sample_meta_fnp

    script:
    def meta_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--groupingsFile ${meta_fnp}"
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    PathWeaver runProcessClustersOnRecon --inputDirectory ${assemblies_root} --genome2bit ${genome_2bit_fnp} --trimBedFnp ${trim_bed} --overWriteDir --dout popClustering --pat _${dirstub} --numThreads ${params.pw_pop_clustering_ncpus} ${meta_arg} ${extra_args}
    elucidator tableExtractCriteria --file popClustering/reports/seqsPerTargetGatheredDetailed.tab.txt --delim tab --columnName nSamples --header  --cutOff 1 | elucidator printCol --file STDIN --delim tab --header --columnName target > popClustering/reports/targets_with_results.txt
    elucidator tableExtractColumns --file popClustering/reports/allSelectedClustersInfo.tab.txt.gz --columns s_Sample::library_sample_name,p_name::target_name,h_popUID::microhaplotype_name,h_Consensus::seq,c_ReadCnt::reads --delim tab --header --out popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz
    ln -s popClustering/reports
    ln -s popClustering/info
    if [ -f popClustering/info/sampleMetaData.tab.txt ]; then ln -s popClustering/info/sampleMetaData.tab.txt output_sampleMetaData.tab.txt; fi;
    """
}
