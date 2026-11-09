#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EXTRACT_REGION_ASSEMBLIES } from '../../../modules/local/extract_region_assemblies'
include { CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS } from '../../../modules/local/concatenate_extract_region_assemblies_results/main.nf'
include { PATHWEAVER_POP_CLUSTERING } from '../../../modules/local/pathweaver_pop_clustering'
include { PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED } from '../../../modules/local/pathweaver_pop_clustering'
include { COMBINE_VARIANT_CALLS } from '../../../modules/local/combine_variant_calls/main.nf'
include { popClusteringVariantCallArgs; noPopClusteringVariantCalls } from '../variant_call_setup.nf'
include { COMBINE_VARIANT_CALLS as COMBINE_VARIANT_CALLS_ON_TRIMMED } from '../../../modules/local/combine_variant_calls/main.nf'
include { WAIT_FOR_ASSEMBLIES } from '../../../modules/local/wait_for_assemblies'


// A sample's assemblies are considered complete (and safe to skip re-running)
// only when these key result files are all present in its output dir.
def pwAssemblyComplete(dir) {
    file("${dir}/final/allFinal.fasta").exists() &&
    file("${dir}/final/basicInfoPerRegion.tab.txt").exists() &&
    file("${dir}/partial/allPartial.fasta").exists()
}


// what population clustering needs to call variants on its final haplotypes, an empty file for the genome dir when variant calling is off
def variantCallSetup(bed_ch, genome_fnp, gff_fnp, known_amino_acid_changes_fnp, pairwise_comps) {
    if (!params.do_variant_calling) {
        return noPopClusteringVariantCalls()
    }
    // Variant calling needs the genome's GFF (derived from the genome file name).
    // Fail early with a clear message instead of letting the variant caller error
    // out downstream on a non-existent annotation file.
    if (!file(gff_fnp).exists()) {
        error """
Variant calling requested (--do_variant_calling) but the genome GFF was not found:
  ${gff_fnp}
Expected at <genome_dir>/info/gff/${file("${genome_fnp}").getBaseName()}.gff. Add the GFF for '${file("${genome_fnp}").getBaseName()}' and re-run.
"""
    }
    def vc_args = popClusteringVariantCallArgs(pairwise_comps)
    return bed_ch.map { bed -> tuple(file("${genome_fnp}").getParent(), file("${genome_fnp}").getBaseName(), file(gff_fnp), file(known_amino_acid_changes_fnp), bed, vc_args) }
}


workflow PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp_ch  // Path to the BED file

    genome_fnp // Path to the genome file
    results_dir // Directory to store results
    meta_fnp //meta file

    main:

    // get possible genome
    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    genome_twobit_fnp = file("${file("${genome_fnp}").getParent()}/${genome_base_name}.2bit")
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }

    // Create output directories
    pop_clus_dir = file("${results_dir}/reports/populationClustering/")
    reports_dir = file("${results_dir}/reports/rawResults/")
    pop_clus_dir.mkdirs()
    reports_dir.mkdirs()
    // Load samples from file and create a channel
    if ("EMPTY_FILE.txt" == file("${samples_file}").name){
        samples = channel.fromPath("${bams_dir}/*${params.bams_file_ending}")
            .map{ samp_file ->
                file(samp_file).name.replaceAll("${params.bams_file_ending}", "")
        }
    } else {
        samples = channel
            .fromPath("${samples_file}")
            .splitText()
            .map{samp ->
                samp.trim() // Remove any whitespace
                }
    }


    // Construct a channel of tuples (BAM, BAI, bed, genome_dir, primary_genome, results_dir)
    input_ch = bed_fnp_ch.combine(samples).map {bed_fnp, samp ->
        tuple(
            samp,
            file("${bams_dir}/${samp}${params.bams_file_ending}"),
            file("${bams_dir}/${samp}${params.bams_file_ending}.bai"),
            file("${bed_fnp}"),
            file("${genome_fnp}").getParent(),
            file("${genome_fnp}").getBaseName(),
            file("${results_dir}").baseName,
            reports_dir.toString()
        )
    }

    // Skip samples whose assemblies already exist: channel the existing output
    // directly instead of re-running the (expensive) extraction + copy-in/out.
    assembly_inputs = input_ch.branch { sample, bam_fnp, bam_bai_fnp, bed_fnp, gdir, primary_genome, dirstub, pub_dir ->
        done: pwAssemblyComplete("${pub_dir}/${sample}_${dirstub}")
            return file("${pub_dir}/${sample}_${dirstub}")
        todo: true
            return tuple(sample, bam_fnp, bam_bai_fnp, bed_fnp, gdir, primary_genome, dirstub, pub_dir)
    }

    // run PathWeaver only on samples that still need it
    EXTRACT_REGION_ASSEMBLIES(assembly_inputs.todo)

    // all assemblies = freshly computed + already-existing
    assemblies_ch = EXTRACT_REGION_ASSEMBLIES.out.sample_results.mix(assembly_inputs.done)

    // Pass the SINGLE parent dir (rawResults) instead of staging tens of thousands of
    // child dirs -- a huge nxf_stage makes bash segfault, and one dir is bind-mounted
    // into containers natively. Because that dir is the publishDir target (populated
    // asynchronously), gate downstream on a completion guard that waits until every
    // expected assembly is fully written before the directory is read.
    expected_names = assemblies_ch.map { d -> d.name }.collectFile(name: 'expected_assembly_dirs.txt', newLine: true)
    WAIT_FOR_ASSEMBLIES(expected_names, reports_dir)
    assemblies_root = WAIT_FOR_ASSEMBLIES.out.ready.first().map { _ready -> reports_dir }

    // concatenate the results files; publish outside rawResults since it's a `path` input
    // here and below, so writing into it changes its mtime and busts -resume
    CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(assemblies_root, "${results_dir}/reports/concatenatedAssemblies".toString(), file("${results_dir}").baseName, file("${samples_file}"))

    // run population clustering, calling variants on the final haplotypes when variant calling is on
    PATHWEAVER_POP_CLUSTERING(file("${results_dir}").baseName, file("${meta_fnp}"), assemblies_root, pop_clus_dir.toString(), file("${samples_file}"), file(params.primers_fnp ?: params.empty_file_fnp), genome_twobit_fnp,
        variantCallSetup(bed_fnp_ch, genome_fnp, gff_fnp, known_amino_acid_changes_fnp, !params.vc_no_pairwise_comps))

    channel.empty().set { variant_calling_reports_dir }

    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/variantCalls")
        meta_fnp_for_variant_calling_ch = channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).name ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{ dir -> file("${dir}/info/sampleMetaData.tab.txt") }
        }
        //the reports across targets from the calls population clustering made
        COMBINE_VARIANT_CALLS (
            bed_fnp_ch,
            PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir,
            "popClustering",
            "variantCalling",
            variant_call_dir.toString(),
            params.vc_extra_args
        )
        // only set this when the module actually runs
        variant_calling_reports_dir = COMBINE_VARIANT_CALLS.out.reports
    }
    emit:
    samples = samples
    assemblies = assemblies_ch
    pop_clustering_res_pop_clustering_dir = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir
    pop_clustering_res_all_selected_clusters_info = PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info
    pop_clustering_res_targets_with_results = PATHWEAVER_POP_CLUSTERING.out.targets_with_results
    variant_calling_reports_dir = variant_calling_reports_dir
}

workflow PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_WITH_TRIM_BED {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp_ch  // Path to the BED file
    trim_bed_fnp_ch  // Path to the BED file
    genome_fnp // Path to the genome file
    results_dir // Directory to store results
    meta_fnp //meta file

    main:

    // get possible genome
    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    genome_dir_info = file("${genome_fnp}").getParent()
    genome_twobit_fnp = file("${genome_dir_info}/${genome_base_name}.2bit")
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }

    // Create output directories
    full_pop_clus_dir = file("${results_dir}/reports/full/populationClustering/")
    trimmed_pop_clus_dir = file("${results_dir}/reports/populationClustering/")
    reports_dir = file("${results_dir}/reports/rawResults/")
    full_pop_clus_dir.mkdirs()
    trimmed_pop_clus_dir.mkdirs()
    reports_dir.mkdirs()
    // Load samples from file and create a channel
    if ("EMPTY_FILE.txt" == file("${samples_file}").name){
        samples = channel.fromPath("${bams_dir}/*${params.bams_file_ending}")
            .map{ samp_file ->
                file(samp_file).name.replaceAll("${params.bams_file_ending}", "")
        }
    } else {
        samples = channel
            .fromPath("${samples_file}")
            .splitText()
            .map{samp ->
                samp.trim() // Remove any whitespace
                }
    }

    // Construct a channel of tuples (BAM, BAI, bed, genome_dir, primary_genome, results_dir)
    input_ch = bed_fnp_ch.combine(samples).map {bed_fnp, samp ->
        tuple(
            samp,
            file("${bams_dir}/${samp}${params.bams_file_ending}"),
            file("${bams_dir}/${samp}${params.bams_file_ending}.bai"),
            file("${bed_fnp}"),
            file("${genome_fnp}").getParent(),
            file("${genome_fnp}").getBaseName(),
            file("${results_dir}").baseName,
            reports_dir.toString()
        )
    }
    // Skip samples whose assemblies already exist: channel the existing output
    // directly instead of re-running the (expensive) extraction + copy-in/out.
    assembly_inputs = input_ch.branch { sample, bam_fnp, bam_bai_fnp, bed_fnp, gdir, primary_genome, dirstub, pub_dir ->
        done: pwAssemblyComplete("${pub_dir}/${sample}_${dirstub}")
            return file("${pub_dir}/${sample}_${dirstub}")
        todo: true
            return tuple(sample, bam_fnp, bam_bai_fnp, bed_fnp, gdir, primary_genome, dirstub, pub_dir)
    }

    // run PathWeaver only on samples that still need it
    EXTRACT_REGION_ASSEMBLIES(assembly_inputs.todo)

    // all assemblies = freshly computed + already-existing
    assemblies_ch = EXTRACT_REGION_ASSEMBLIES.out.sample_results.mix(assembly_inputs.done)

    // Pass the SINGLE parent dir (rawResults) instead of staging tens of thousands of
    // child dirs -- a huge nxf_stage makes bash segfault, and one dir is bind-mounted
    // into containers natively. Because that dir is the publishDir target (populated
    // asynchronously), gate downstream on a completion guard that waits until every
    // expected assembly is fully written before the directory is read.
    expected_names = assemblies_ch.map { d -> d.name }.collectFile(name: 'expected_assembly_dirs.txt', newLine: true)
    WAIT_FOR_ASSEMBLIES(expected_names, reports_dir)
    assemblies_root = WAIT_FOR_ASSEMBLIES.out.ready.first().map { _ready -> reports_dir }

    // concatenate the results files; publish outside rawResults since it's a `path` input
    // here and below, so writing into it changes its mtime and busts -resume
    CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(assemblies_root, "${results_dir}/reports/concatenatedAssemblies".toString(), file("${results_dir}").baseName, file("${samples_file}"))

    //run population clustering, calling variants on the final haplotypes when variant calling is on
    //run with full region
    PATHWEAVER_POP_CLUSTERING(              file("${results_dir}").baseName, file("${meta_fnp}"), assemblies_root, full_pop_clus_dir.toString(), file("${samples_file}"), file(params.primers_fnp ?: params.empty_file_fnp), genome_twobit_fnp,
                                            variantCallSetup(bed_fnp_ch, genome_fnp, gff_fnp, known_amino_acid_changes_fnp, false))
    //run with inner trimmed region
    PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED(file("${results_dir}").baseName, file("${meta_fnp}"), assemblies_root, trimmed_pop_clus_dir.toString(),
                                            trim_bed_fnp_ch, genome_twobit_fnp, file("${samples_file}"), file(params.primers_fnp ?: params.empty_file_fnp),
                                            variantCallSetup(trim_bed_fnp_ch, genome_fnp, gff_fnp, known_amino_acid_changes_fnp, !params.vc_no_pairwise_comps))
    channel.empty().set { variant_calling_reports_dir }
    channel.empty().set { trimmed_variant_calling_reports_dir }
    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/reports/full/variantCalls")
        trimmed_variant_call_dir = file("${results_dir}/variantCalls")
        meta_fnp_for_variant_calling_ch = channel.fromPath(params.meta_fnp)
        trimmed_meta_fnp_for_variant_calling_ch =  channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).name ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{ dir -> file("${dir}/info/sampleMetaData.tab.txt") }
            trimmed_meta_fnp_for_variant_calling_ch = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.pop_clustering_dir.map{ dir -> file("${dir}/info/sampleMetaData.tab.txt") }
        }
        //the reports across targets from the calls population clustering made
        COMBINE_VARIANT_CALLS (
            bed_fnp_ch,
            PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir,
            "popClustering",
            "variantCalling",
            variant_call_dir.toString(),
            params.vc_extra_args
        )
        COMBINE_VARIANT_CALLS_ON_TRIMMED (
            trim_bed_fnp_ch,
            PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            trimmed_meta_fnp_for_variant_calling_ch,
            PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.pop_clustering_dir,
            "popClustering",
            "variantCalling",
            trimmed_variant_call_dir.toString(),
            params.vc_extra_args
        )
        // rebind emits only when present
        variant_calling_reports_dir         = COMBINE_VARIANT_CALLS.out.reports
        trimmed_variant_calling_reports_dir = COMBINE_VARIANT_CALLS_ON_TRIMMED.out.reports
    }
    emit:
    samples = samples
    assemblies = assemblies_ch
    pop_clustering_res_pop_clustering_dir = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir
    pop_clustering_res_all_selected_clusters_info = PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info
    pop_clustering_res_targets_with_results = PATHWEAVER_POP_CLUSTERING.out.targets_with_results
    variant_calling_reports_dir = variant_calling_reports_dir

    trimmed_pop_clustering_res_pop_clustering_dir = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.pop_clustering_dir
    trimmed_pop_clustering_res_all_selected_clusters_info = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.all_selected_clusters_info
    trimmed_pop_clustering_res_targets_with_results = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.targets_with_results
    trimmed_variant_calling_reports_dir = trimmed_variant_calling_reports_dir
}
