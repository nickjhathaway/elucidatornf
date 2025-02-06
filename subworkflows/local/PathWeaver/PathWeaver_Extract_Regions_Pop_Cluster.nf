#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {EXTRACT_REGION_ASSEMBLIES} from '../../../modules/local/extract_region_assemblies'
include {PATHWEAVER_POP_CLUSTERING} from '../../../modules/local/pathweaver_pop_clustering'


workflow PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp  // Path to the BED file
    genome_fnp // Path to the genome file 
    results_dir // Directory to store results

    main:


    // Create output directory if not exists
    new File(results_dir).mkdirs()

    // Load samples from file and create a channel
    samples = Channel
        .fromPath("${samples_file}")
        .splitText()

    // Construct a channel of tuples (BAM, BAI, sample)
    input_ch = samples.map { samp -> 
        samp = samp.trim()  // Remove any whitespace
        tuple(
            samp,
            file("${bams_dir}/${samp}.sorted.bam"), 
            file("${bams_dir}/${samp}.sorted.bam.bai"), 
            file("${bed_fnp}"),
            file("${genome_fnp}").getParent(),
            file("${genome_fnp}").getBaseName(), 
            file("${results_dir}").baseName
        )
    }

    // run PathWeaver on each 
    PATHWEAVER_POP_CLUSTERING(file("${results_dir}").baseName,EXTRACT_REGION_ASSEMBLIES(input_ch) | collect)
}
