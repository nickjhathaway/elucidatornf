// the arguments population clustering gets to call variants on its final haplotypes the same way VARIANT_CALL_ON_HAP_TABLE does
def popClusteringVariantCallArgs(pairwise_comps) {
    return "--lowVariantCutOff ${params.vc_variant_frequency_cut_off} --occurrenceCutOff ${params.vc_variant_occurrence_cut_off}" +
        " --ignoreSubFields site:LabCross,site:LabControl,site:LabContaminated --exportLabIsolateSeqs" +
        (pairwise_comps ? " --getPairwiseComps" : "") +
        ("" == params.meta_fields_to_calc_pop_diffs ? "" : " --metaFieldsToCalcPopDiffs ${params.meta_fields_to_calc_pop_diffs}")
}

// the variant calling input of population clustering when variant calling is off, an empty file for the genome dir
def noPopClusteringVariantCalls() {
    def empty = file(params.empty_file_fnp)
    return channel.value(tuple(empty, "", empty, empty, empty, ""))
}

// genomes a target in several genomes is placed in first for the PMO, the variant calling genome after any given ones so the target locations match the calls
def pmoPreferredGenomes() {
    def genomes = params.pmo_preferred_genomes ? params.pmo_preferred_genomes.toString().split(",").collect { g -> g.trim() }.findAll { g -> g } : []
    def vc_genome = params.vc_primary_genome ? params.vc_primary_genome.toString().trim() : ""
    if (vc_genome && vc_genome.toLowerCase() != "false" && !genomes.contains(vc_genome)) {
        genomes << vc_genome
    }
    return genomes.join(",")
}
