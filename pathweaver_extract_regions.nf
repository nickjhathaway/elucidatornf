#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include { SETUP_PLASMODIUM_GENOMES } from './subworkflows/local/plasmodium_genomes.nf'
include { PATHWEAVER_EXTRACT_REGIONS_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'
include { PATHWEAVER_EXTRACT_REGIONS_WITH_PRIMERS_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'
include { PATHWEAVER_EXTRACT_REGIONS_WITH_SEQS_TABLE_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'
include { PATHWEAVER_EXTRACT_REGIONS_WITH_GENE_IDS_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'



workflow {

main:
    // Validate inputs
    def download_genomes = PlasmodiumGenomes.enabled(params)
    if (params.bams_dir == null  || (params.genome_fnp == null && !download_genomes) || params.outdir == null) {
        error "flags '--bams_dir', '--genome_fnp' (or '--download_plasmodium_genomes'), and '--outdir' must be specified!"
    }

    sample_file = params.pw_samples_file
    if (params.pw_samples_file == null){
        sample_file = file(params.empty_file_fnp)
    }

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (params.pw_overwrite_dir && results_dir_obj.exists()){
        // keep genomes downloaded into the output directory
        OutputDir.clear(results_dir_obj, download_genomes ? PlasmodiumGenomes.rootDir(params) : null, params.trace_timestamp, [workflow.launchDir, workflow.workDir])
    } else if (file("${results_dir_obj}/pipeline_info").exists()){
        //if the output does exist and are not overwirte, do remove the pipeline run info so the new run time info is logged
        // only the earlier runs' files, nextflow already has this run's trace open in it (on NFS deleting it leaves a
        // .nfs placeholder that stops pipeline_info from being removed)
        OutputDir.clear(file("${results_dir_obj}/pipeline_info"), null, params.trace_timestamp)
        //since the sub regions can change based on the input, have to delete the sub region 
        sub_region_info_dir = file("${results_dir_obj}/PathWeaverResults/subVarRegions")
        if(sub_region_info_dir.exists()){
            sub_region_info_dir.deleteDir()
        }
    }

    results_dir_obj.mkdirs()

    if(params.primers_fnp == null && params.pw_bed_fnp == null && params.pw_seqs_table == null && params.pw_gene_ids == null){
        error "at least one of the following flags '--primers_fnp' , '--pw_seqs_table', '--pw_gene_ids' or '--pw_bed_fnp' must be specified!"
    }

    // --download_plasmodium_genomes, --genome_fnp is then the name of the downloaded genome to use (default Pf3D7)
    def genome_fnp = params.genome_fnp
    def genomes_ready = channel.value(true)
    if (download_genomes) {
        def genomes_root = PlasmodiumGenomes.rootDir(params)
        SETUP_PLASMODIUM_GENOMES(genomes_root, params.plasmodium_genomes_version)
        def genome_name = params.genome_fnp ? file(params.genome_fnp.toString()).name.replaceAll(/\.fasta$/, "") : PlasmodiumGenomes.DEFAULT_GENOME
        if (!PlasmodiumGenomes.genomeNames(genomes_root).contains(genome_name)) {
            error "--genome_fnp '${genome_name}' isn't one of the downloaded genomes in ${PlasmodiumGenomes.genomeDir(genomes_root)}: ${PlasmodiumGenomes.genomeNames(genomes_root).join(", ")}"
        }
        genome_fnp = PlasmodiumGenomes.genomeDir(genomes_root).resolve("${genome_name}.fasta").toString()
        genomes_ready = SETUP_PLASMODIUM_GENOMES.out.ready
    }
    if( params.primers_fnp != null){
        PATHWEAVER_EXTRACT_REGIONS_WITH_PRIMERS_FULL(sample_file, params.bams_dir, params.primers_fnp, genome_fnp, params.outdir, params.meta_fnp, genomes_ready)
    } else if (params.pw_bed_fnp != null){
        PATHWEAVER_EXTRACT_REGIONS_FULL(sample_file, params.bams_dir, params.pw_bed_fnp, genome_fnp, params.outdir, params.meta_fnp, genomes_ready)
    } else if (params.pw_gene_ids != null){
        PATHWEAVER_EXTRACT_REGIONS_WITH_GENE_IDS_FULL(sample_file, params.bams_dir, params.pw_gene_ids, genome_fnp, params.outdir, params.meta_fnp, genomes_ready)
    } else if (params.pw_seqs_table != null){
        if (params.pw_seqs_table_seqs_col == null  || params.pw_seqs_table_name_col == null || params.pw_seqs_table_target_col == null) {
            error "flags '--pw_seqs_table_seqs_col', '--pw_seqs_table_name_col', and '--pw_seqs_table_target_col' must be specified if supplying '--pw_seqs_table'!"
        }
        PATHWEAVER_EXTRACT_REGIONS_WITH_SEQS_TABLE_FULL(
            sample_file,
            params.bams_dir,
            params.pw_seqs_table,
            params.pw_seqs_table_seqs_col,
            params.pw_seqs_table_name_col,
            params.pw_seqs_table_target_col,
            genome_fnp,
            params.outdir,
            params.meta_fnp,
            genomes_ready)
    }
}

