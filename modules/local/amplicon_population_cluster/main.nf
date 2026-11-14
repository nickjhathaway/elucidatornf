process AMPLICON_POPULATION_CLUSTERING {
    tag "${target_name}_pop_clustering"
    label 'process_medium'

    cpus   { ncpus }

    input:
    // vc_genome is the genome the target is variant called against, blank to not call it
    tuple path(fastqs_fnp), val(target_name), val(ncpus), path(meta_fnp, name: 'meta.tsv'), path(previous_pop_dir), val(sample_min_read_count), path(primers_fnp, name: 'primers.tsv'), val(panel_name), val(vc_genome)
    // variant calling on the final haplotypes, the same as VARIANT_CALL_ON_HAP_TABLE, an empty file for the genome dir turns it off
    // the known amino acid changes are for vc_primary_genome, the beds dir has <genome>_variantCall.bed from ASSIGN_TARGET_GENOMES
    tuple path(vc_genome_dir), val(vc_primary_genome), path(vc_gff_dir, name: 'vc_gffs'), path(vc_known_amino_acid_changes_fnp, name: 'vc_known_amino_acid_changes.tsv'), path(vc_beds_dir, name: 'vc_beds'), val(vc_args)


    output:
    path "${target_name}_selectedClustersInfo.tab.txt.gz", emit: clustering_results
    path "${target_name}_sampleClustersInfo.tab.txt.gz", emit: sample_clusters_info
    path "${target_name}_PMO.json.gz", emit: pmo
    tuple val(vc_genome), path("variantCalls/${target_name}"), optional: true, emit: variant_calls

    script:
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    meta_arg="--groupingsFile ${meta_fnp}"
    previous_pop_arg="--previousPop ${previous_pop_dir}/${target_name}.fasta"
    if [ ! -s "\$(readlink -f ${meta_fnp})" ]; then meta_arg=""; fi
    if [ ! -f ${previous_pop_dir}/${target_name}.fasta ]; then previous_pop_arg=""; fi
    variant_call_args=""
    # the genome dir keeps its name so the vcf reference lines match VARIANT_CALL_ON_HAP_TABLE's
    if [ -d "${vc_genome_dir}" ] && [ -n "${vc_genome}" ]; then
        variant_call_args="--genome ${vc_genome_dir}/${vc_genome}.fasta --gff vc_gffs/${vc_genome}.gff --genomicLocations vc_beds/${vc_genome}_variantCall.bed ${vc_args}"
        if [ "${vc_genome}" == "${vc_primary_genome}" ] && [ -s "\$(readlink -f vc_known_amino_acid_changes.tsv)" ]; then variant_call_args="\${variant_call_args} --knownAminoAcidChangesFnp vc_known_amino_acid_changes.tsv"; fi
    fi
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
            \${variant_call_args} \
            ${extra_args}
    ln -s analysis/selectedClustersInfo.tab.txt.gz ${target_name}_selectedClustersInfo.tab.txt.gz
    ln -s analysis/sampleClustersInfo.tab.txt.gz ${target_name}_sampleClustersInfo.tab.txt.gz
    ln -s analysis/PMO.json.gz ${target_name}_PMO.json.gz
    # targets not in the genomic locations bed aren't variant called
    if [ -d analysis/variantCalling ]; then mkdir -p variantCalls && ln -s ../analysis/variantCalling variantCalls/${target_name}; fi
    """
}
