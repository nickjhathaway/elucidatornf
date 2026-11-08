process AMPLICON_POPULATION_CLUSTERING {
    tag "${target_name}_pop_clustering"
    label 'process_medium'

    cpus   { ncpus }

    input:
    tuple path(fastqs_fnp), val(target_name), val(ncpus), path(meta_fnp, name: 'meta.tsv'), path(previous_pop_dir), val(sample_min_read_count), path(primers_fnp, name: 'primers.tsv'), val(panel_name)


    output:
    path "${target_name}_selectedClustersInfo.tab.txt.gz", emit: clustering_results
    path "${target_name}_sampleClustersInfo.tab.txt.gz", emit: sample_clusters_info
    path "${target_name}_PMO.json.gz", emit: pmo

    script:
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    meta_arg="--groupingsFile ${meta_fnp}"
    previous_pop_arg="--previousPop ${previous_pop_dir}/${target_name}.fasta"
    if [ ! -s "\$(readlink -f ${meta_fnp})" ]; then meta_arg=""; fi
    if [ ! -f ${previous_pop_dir}/${target_name}.fasta ]; then previous_pop_arg=""; fi
    # targets clustered with separate mates are named <target>-R1 and <target>-R2 but have the one set of primers
    primer_target="${target_name}"
    if ! awk '{print \$1}' primers.tsv | grep -qxF "${target_name}"; then primer_target=\$(echo "${target_name}" | sed -E 's/-R[12]\$//'); fi
    SeekDeep processClusters \
            --flatMasterDir \
            --strictErrors \
            --dout analysis \
            --fastqgz output.fastq.gz \
            --allowHomopolymerCollapse \
            --replicateMinTotalReadCutOff ${sample_min_read_count} \
            \${meta_arg} \
            --experimentName ${target_name} \
            --numThreads ${ncpus}\
            --previousPop-largeBaseIndel 0.99 --previousPop-oneBaseIndel 0.99 --previousPop-twoBaseIndel 0.99 \
            \${previous_pop_arg}\
            --writePMO --pmoPrimers primers.tsv --pmoPrimerTarget \${primer_target} --pmoPanelName ${panel_name} \
            ${extra_args}
    ln -s analysis/selectedClustersInfo.tab.txt.gz ${target_name}_selectedClustersInfo.tab.txt.gz
    ln -s analysis/sampleClustersInfo.tab.txt.gz ${target_name}_sampleClustersInfo.tab.txt.gz
    ln -s analysis/PMO.json.gz ${target_name}_PMO.json.gz
    """
}
