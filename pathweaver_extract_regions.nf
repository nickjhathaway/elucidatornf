#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include { PATHWEAVER_EXTRACT_REGIONS_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'
include { PATHWEAVER_EXTRACT_REGIONS_WITH_PRIMERS_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'
include { PATHWEAVER_EXTRACT_REGIONS_WITH_SEQS_TABLE_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'
include { PATHWEAVER_EXTRACT_REGIONS_WITH_GENE_IDS_FULL } from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'



workflow {

main:
    // Validate inputs
    if (params.bams_dir == null  || params.genome_fnp == null || params.outdir == null) {
        error "flags '--bams_dir', '--genome_fnp', and '--outdir' must be specified!"
    }

    sample_file = params.pw_samples_file
    if (params.pw_samples_file == null){
        sample_file = file(params.empty_file_fnp)
    }

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (params.pw_overwrite_dir && results_dir_obj.exists()){
        results_dir_obj.deleteDir()
    } else if (file("${results_dir_obj}/pipeline_info").exists()){
        //if the output does exist and are not overwirte, do remove the pipeline run info so the new run time info is logged
        file("${results_dir_obj}/pipeline_info").deleteDir()
    }

    results_dir_obj.mkdirs()
    if(params.primers_fnp == null && params.pw_bed_fnp == null && params.pw_seqs_table == null && params.pw_gene_ids == null){
        error "at least one of the following flags '--primers_fnp' , '--pw_seqs_table', '--pw_gene_ids' or '--pw_bed_fnp' must be specified!"
    }
    if( params.primers_fnp != null){
        PATHWEAVER_EXTRACT_REGIONS_WITH_PRIMERS_FULL(sample_file, params.bams_dir, params.primers_fnp, params.genome_fnp, params.outdir, params.meta_fnp)
    } else if (params.pw_bed_fnp != null){
        PATHWEAVER_EXTRACT_REGIONS_FULL(sample_file, params.bams_dir, params.pw_bed_fnp, params.genome_fnp, params.outdir, params.meta_fnp)
    } else if (params.pw_gene_ids != null){
        PATHWEAVER_EXTRACT_REGIONS_WITH_GENE_IDS_FULL(sample_file, params.bams_dir, params.pw_gene_ids, params.genome_fnp, params.outdir, params.meta_fnp)
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
            params.genome_fnp,
            params.outdir,
            params.meta_fnp)
    }
}

