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
    path samples_fnp, name: "pw_samples.txt"
    path primers_fnp, name: "pmo_primers.tsv"
    path genome_2bit_fnp
    // variant calling on the final haplotypes, the same as VARIANT_CALL_ON_HAP_TABLE, an empty file for the genome dir turns it off
    tuple path(vc_genome_dir), val(vc_primary_genome), path(vc_gff_fnp, name: 'vc_genome.gff'), path(vc_known_amino_acid_changes_fnp, name: 'vc_known_amino_acid_changes.tsv'), path(vc_bed_fnp, name: 'vc_genomic_locations.bed'), val(vc_args)

    output:
    path "popClustering/", emit: pop_clustering_dir
    path "popClustering/reports/allSelectedClustersInfo.tab.txt.gz", emit: all_selected_clusters_info
    path "popClustering/reports/allSampleClustersInfo.tab.txt.gz", emit: all_sample_clusters_info
    path "popClustering/reports/allPMO.json.gz", emit: all_pmo
    path "popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz", emit: all_slim_selected_clusters_info
    path "reports/", emit: pop_clustering_reports_dir
    path "info/", emit: pop_clustering_info_dir
    path "popClustering/reports/targets_with_results.txt", emit: targets_with_results
    path "output_sampleMetaData.tab.txt", optional: true, emit: output_sample_meta_fnp

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def meta_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--groupingsFile ${meta_fnp}"
    """
    samples_arg="--samples ${samples_fnp}"
    if [ ! -s "\$(readlink -f ${samples_fnp})" ]; then samples_arg=""; fi
    primers_arg="--pmoPrimers pmo_primers.tsv"
    if [ ! -s "\$(readlink -f pmo_primers.tsv)" ]; then primers_arg=""; fi
    variant_call_args=""
    # the genome dir keeps its name so the vcf reference lines match VARIANT_CALL_ON_HAP_TABLE's
    if [ -d "${vc_genome_dir}" ]; then
        variant_call_args="--genome ${vc_genome_dir}/${vc_primary_genome}.fasta --gff vc_genome.gff --genomicLocations vc_genomic_locations.bed ${vc_args}"
        if [ -s "\$(readlink -f vc_known_amino_acid_changes.tsv)" ]; then variant_call_args="\${variant_call_args} --knownAminoAcidChangesFnp vc_known_amino_acid_changes.tsv"; fi
    fi
    PathWeaver runProcessClustersOnRecon --writePMO \${primers_arg} --pmoPanelName ${params.pmo_panel_name} --genome2bit ${genome_2bit_fnp} --inputDirectory ${assemblies_root} --overWriteDir --dout popClustering --pat _${dirstub} --numThreads ${params.pw_pop_clustering_ncpus} ${meta_arg} \${samples_arg} \${variant_call_args} ${extra_args}
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
    path samples_fnp, name: "pw_samples.txt"
    path primers_fnp, name: "pmo_primers.tsv"
    // variant calling on the final haplotypes, the same as VARIANT_CALL_ON_HAP_TABLE, an empty file for the genome dir turns it off
    tuple path(vc_genome_dir), val(vc_primary_genome), path(vc_gff_fnp, name: 'vc_genome.gff'), path(vc_known_amino_acid_changes_fnp, name: 'vc_known_amino_acid_changes.tsv'), path(vc_bed_fnp, name: 'vc_genomic_locations.bed'), val(vc_args)

    output:
    path "popClustering/", emit: pop_clustering_dir
    path "popClustering/reports/allSelectedClustersInfo.tab.txt.gz", emit: all_selected_clusters_info
    path "popClustering/reports/allSampleClustersInfo.tab.txt.gz", emit: all_sample_clusters_info
    path "popClustering/reports/allPMO.json.gz", emit: all_pmo
    path "popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz", emit: all_slim_selected_clusters_info
    path "reports/", emit: pop_clustering_reports_dir
    path "info/", emit: pop_clustering_info_dir
    path "popClustering/reports/targets_with_results.txt", emit: targets_with_results
    path "output_sampleMetaData.tab.txt", optional: true, emit: output_sample_meta_fnp

    script:
    def meta_arg = "EMPTY_FILE.txt" == "${meta_fnp}" ? "" : "--groupingsFile ${meta_fnp}"
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    samples_arg="--samples ${samples_fnp}"
    if [ ! -s "\$(readlink -f ${samples_fnp})" ]; then samples_arg=""; fi
    primers_arg="--pmoPrimers pmo_primers.tsv"
    if [ ! -s "\$(readlink -f pmo_primers.tsv)" ]; then primers_arg=""; fi
    variant_call_args=""
    # the genome dir keeps its name so the vcf reference lines match VARIANT_CALL_ON_HAP_TABLE's
    if [ -d "${vc_genome_dir}" ]; then
        variant_call_args="--genome ${vc_genome_dir}/${vc_primary_genome}.fasta --gff vc_genome.gff --genomicLocations vc_genomic_locations.bed ${vc_args}"
        if [ -s "\$(readlink -f vc_known_amino_acid_changes.tsv)" ]; then variant_call_args="\${variant_call_args} --knownAminoAcidChangesFnp vc_known_amino_acid_changes.tsv"; fi
    fi
    PathWeaver runProcessClustersOnRecon --writePMO \${primers_arg} --pmoPanelName ${params.pmo_panel_name} --inputDirectory ${assemblies_root} --genome2bit ${genome_2bit_fnp} --trimBedFnp ${trim_bed} --overWriteDir --dout popClustering --pat _${dirstub} --numThreads ${params.pw_pop_clustering_ncpus} ${meta_arg} \${samples_arg} \${variant_call_args} ${extra_args}
    elucidator tableExtractCriteria --file popClustering/reports/seqsPerTargetGatheredDetailed.tab.txt --delim tab --columnName nSamples --header  --cutOff 1 | elucidator printCol --file STDIN --delim tab --header --columnName target > popClustering/reports/targets_with_results.txt
    elucidator tableExtractColumns --file popClustering/reports/allSelectedClustersInfo.tab.txt.gz --columns s_Sample::library_sample_name,p_name::target_name,h_popUID::microhaplotype_name,h_Consensus::seq,c_ReadCnt::reads --delim tab --header --out popClustering/reports/slim_allSelectedClustersInfo.tab.txt.gz
    ln -s popClustering/reports
    ln -s popClustering/info
    if [ -f popClustering/info/sampleMetaData.tab.txt ]; then ln -s popClustering/info/sampleMetaData.tab.txt output_sampleMetaData.tab.txt; fi;
    """
}
