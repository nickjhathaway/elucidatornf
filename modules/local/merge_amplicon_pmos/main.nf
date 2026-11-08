process MERGE_AMPLICON_PMOS {
    label 'process_single'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "allPMO.json.gz"

    input:
    path pmo_fnps
    path(meta_fnp, name: 'meta.tsv')
    val pub_dir


    output:
    path "allPMO.json.gz", emit: pmo

    script:
    // the per target PMOs are over the same samples, the sample meta is added to the specimens of the merged PMO
    """
    SeekDeep mergeMinimumPmos --contains _PMO.json.gz --out merged.json.gz
    if [ -s "\$(readlink -f ${meta_fnp})" ]; then
        SeekDeep addSpecimenMetaToPmo --pmo merged.json.gz --meta ${meta_fnp} --out allPMO.json.gz
    else
        mv merged.json.gz allPMO.json.gz
    fi
    """
}
