#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {PATHWEAVER_EXTRACT_REGIONS_FULL} from './subworkflows/local/PathWeaver/PathWeaver_Extract_Regions_Pop_Cluster.nf'


workflow {

  main:
  // Validate inputs
  if (params.pw_samples_file == null || params.bams_dir == null  || params.pw_bed_fnp == null || params.genome_fnp == null || params.pw_results_dir == null) {
      error "flags '--pw_samples_file', '--bams_dir', '--pw_bed_fnp', '--genome_fnp', and '--pw_results_dir' must be specified!"
  }

  PATHWEAVER_EXTRACT_REGIONS_FULL(params.pw_samples_file, params.bams_dir, params.pw_bed_fnp, params.genome_fnp, params.pw_results_dir, params.meta_fnp)
  
}

