#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {INDEX_GENOMES_DIR} from './modules/local/index_genomes_dir'


workflow {

main:
    // Validate inputs
    if (params.genome_dir == null) {
        error "flag '--genome_dir' must be specified!"
    }
    INDEX_GENOMES_DIR(file(params.genome_dir), params.max_cpus)
}

