// which genome each target is variant called against, and a bed of each genome's targets for --genomicLocations
process ASSIGN_TARGET_GENOMES {
    label 'process_single'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true

    input:
    path locations_by_genome // locationsByGenome of genTargetInfoFromGenomes
    val genomes // comma separated, a target goes to the first it's in
    val pub_dir

    output:
    path "assignedTargetGenomes", emit: assigned

    script:
    """
    SeekDeep assignTargetGenomes --locations ${locations_by_genome} --genomes ${genomes} --dout assignedTargetGenomes --verbose
    """
}
