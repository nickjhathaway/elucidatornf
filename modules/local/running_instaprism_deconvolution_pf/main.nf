process RUN_INSTAPRISM_DECONVOLUTION_PF {

    publishDir "${pubdir}", mode: 'copy', overwrite: true

    tag "run_instaprism_deconvolution_pf"

    input:
    tuple path(all_quant_fnp), path(single_cell_r_object), val(pubdir)

    output:
    path("pf_instaprism_fine_cell_proportions.tsv"), emit: pf_instaprism_fine_cell_proportions_fnp
    path("pf_instaprism_course_cell_proportions.tsv"), emit: pf_instaprism_fine_cell_proportions_fnp

    script:
    //right now code is very specific to pf malaria cell atlas data, @todo make more flexible
    """
    running_instaprism_deconvolution_pf.R \
            --all_quants ${all_quant_fnp} \
            --coarse_cell_label STAGE_LR \
            --fine_cell_label STAGE_HR2 \
            --output_course pf_instaprism_fine_cell_proportions.tsv \
            --output_fine pf_instaprism_fine_cell_proportions.tsv \
            --single_cell_r_object ${single_cell_r_object}\
            --use_top_diff \
            --top_diff_gene_amount 200
    """
}
