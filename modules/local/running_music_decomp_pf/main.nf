process RUN_MUSIC_DECOMP_PF {

    publishDir { "${pubdir}" }, mode: 'copy', overwrite: true

    tag "run_music_decomp_pf"

    input:
    tuple path(all_quant_fnp), path(single_cell_r_object), val(pubdir)

    output:
    path("pf_music_decomp_cell_proportions.tsv"), emit: pf_music_decomp_cell_proportions_fnp

    script:
    """
    running_music_decomp_pf.R --single_cell_r_object  ${single_cell_r_object}\
        --all_quants ${all_quant_fnp} \
        --output pf_music_decomp_cell_proportions.tsv
    """
}
