process CONCATENATE_SUB_SEGMENT_LOCS {
    label 'process_single'

    publishDir "${pub_results_dir}", mode: 'copy', overwrite: true


    input:
    path variable_genomic_loc_fnps
    path conserved_genomic_loc_fnps
    path gff_fnp
    val pub_results_dir

    output:
    path "allVariableRegions.bed", emit: all_variable_regions
    path "allConservedRegions.bed", emit: all_conserved_regions
    path "combinedSubRegions.bed", emit: combined_all_sub_regions

    script:
    """
    elucidator rBind --contains _0_ref_variable_expanded_genomic.bed --delim tab | elucidator bedCoordSort --bed STDIN --out raw_allVariableRegions.bed
    elucidator rBind --contains _0_ref_sharedLocs_genomic.bed --delim tab | elucidator bedCoordSort --bed STDIN --out raw_allConservedRegions.bed

    elucidator bedGetIntersectingGenesInGff --gff ${gff_fnp} --extraAttributes description --overWrite --bed raw_allVariableRegions.bed --out allVariableRegions.bed
    elucidator bedGetIntersectingGenesInGff --gff ${gff_fnp} --extraAttributes description --overWrite --bed raw_allConservedRegions.bed --out allConservedRegions.bed

    cat allVariableRegions.bed allConservedRegions.bed  | elucidator bedCoordSort --bed STDIN --out  combinedSubRegions.bed
    """
}


