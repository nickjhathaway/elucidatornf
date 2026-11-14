process MERGE_AMPLICON_PMOS {
    label 'process_single'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "allPMO.json.gz"

    input:
    path pmo_fnps
    path(meta_fnp, name: 'meta.tsv')
    path locations_by_genome // locationsByGenome of genTargetInfoFromGenomes, which genome each target is in
    path genome_dir
    path gff_dir
    val preferred_genomes // comma separated genomes to place targets in first when they're in several, can be blank
    val pub_dir


    output:
    path "allPMO.json.gz", emit: pmo

    script:
    // the per target PMOs are over the same samples, each target gets its genome, location and genes, and the sample meta is added to the specimens
    def preferred_arg = preferred_genomes ? "--preferredGenomes ${preferred_genomes}" : ""
    """
    SeekDeep mergeMinimumPmos --contains _PMO.json.gz --out merged.json.gz
    SeekDeep addTargetGenomeInfoToPmo --pmo merged.json.gz --locations ${locations_by_genome} --genomeDir ${genome_dir} --gffDir ${gff_dir} ${preferred_arg} --out withGenomes.json.gz
    if [ -s "\$(readlink -f ${meta_fnp})" ]; then
        SeekDeep addSpecimenMetaToPmo --pmo withGenomes.json.gz --meta ${meta_fnp} --out allPMO.json.gz
    else
        mv withGenomes.json.gz allPMO.json.gz
    fi
    """
}
