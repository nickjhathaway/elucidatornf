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
