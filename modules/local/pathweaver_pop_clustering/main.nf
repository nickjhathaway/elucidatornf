process PATHWEAVER_POP_CLUSTERING {
    tag 'process_high'

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    publishDir "${params.pw_results_dir}", mode: 'copy', overwrite: true, pattern: "popClustering"

    input:
    val dirstub
    path res_folders
    
    output:
    path "popClustering", emit: pop_clustering_dir


    script:
    def meta_arg = params.meta_fnp == null ? "" : "--groupingsFile ${params.meta_fnp}"
    """
    PathWeaver runProcessClustersOnRecon --overWriteDir --dout popClustering --pat _${dirstub} --numThreads ${params.pw_pop_clustering_ncpus} ${meta_arg}

    """
}
