#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nickjhathaway/elucidatornf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This is a MULTI-WORKFLOW project: there is no single default pipeline.
    This main.nf does NOT run any analysis -- it just lists the available entry
    workflows so you can pick one. It also exists so nf-core tooling
    (e.g. `nf-core pipelines schema build`) recognizes this repo as a pipeline.
----------------------------------------------------------------------------------------
*/

workflow {
    main:
    log.info(
"""
==================================================================================
  nickjhathaway/elucidatornf -- available entry workflows
==================================================================================
  This main.nf does not run anything. Run one of the workflows below directly
  with `nextflow run <file> [params]` (pass no params to see its required flags):

    illumina_clustering.nf          Illumina amplicon clustering (SeekDeep)
    nanopore_clustering.nf          Nanopore amplicon clustering (SeekDeep)
    pathweaver_extract_regions.nf   PathWeaver region extraction + population clustering
    wgs_preprocess_pf.nf            WGS preprocessing (P. falciparum)
    rnaseq_process_pf.nf            RNA-seq processing (P. falciparum)
    index_genomes.nf                Index a directory of genomes

  Example:
    nextflow run pathweaver_extract_regions.nf --bams_dir bams/ --genome_fnp ref.fasta --pw_bed_fnp targets.bed --outdir results
==================================================================================
"""
    )
}
