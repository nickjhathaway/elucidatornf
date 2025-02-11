process AMPLICON_POPULATION_CLUSTERING {
    label 'process_medium'

    cpus   { ncpus }

    input:
    tuple path(fastqs_fnp), val(target_name), val(ncpus), path(meta_fnp, name: 'meta.tsv'), path(previous_pop_fnp, name: 'previous_pop.fasta'), val(sample_min_read_count)


    output:
    path "${target_name}_selectedClustersInfo.tab.txt.gz", emit: clustering_results

    script:

    """
    meta_arg="--groupingsFile ${meta_fnp}"
    previous_pop_arg="--previousPop ${previous_pop_fnp}"
    if [ ! -s "\$(readlink -f ${meta_fnp})" ]; then meta_arg=""; fi
    if [ ! -s "\$(readlink -f ${previous_pop_fnp})" ]; then previous_pop_arg=""; fi
    SeekDeep processClusters \
            --flatMasterDir \
            --strictErrors \
            --dout analysis \
            --fastqgz output.fastq.gz \
            --allowHomopolymerCollapse \
            --removeOneSampOnlyOneOffHaps \
            --excludeCommonlyLowFreqHaplotypes \
            --excludeLowFreqOneOffs \
            --rescueExcludedOneOffLowFreqHaplotypes \
            --replicateMinTotalReadCutOff ${sample_min_read_count} \
            \${meta_arg} \
            --experimentName ${target_name} \
            --numThreads ${ncpus}\
            \${previous_pop_arg}
    ln -s analysis/selectedClustersInfo.tab.txt.gz ${target_name}_selectedClustersInfo.tab.txt.gz
    """
}
