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

// a comma separated param as a list
def commaListParam(value) {
    return value ? value.toString().split(",").collect { entry -> entry.trim() }.findAll { entry -> entry } : []
}

// the genomes targets are variant called against in order, a target goes to the first it's in: vc_primary_genome then vc_other_genomes
def variantCallGenomes() {
    def genomes = []
    def vc_genome = params.vc_primary_genome ? params.vc_primary_genome.toString().trim() : ""
    if (vc_genome && vc_genome.toLowerCase() != "false") {
        genomes << vc_genome
    }
    commaListParam(params.vc_other_genomes).each { genome ->
        if (!genomes.contains(genome)) {
            genomes << genome
        }
    }
    return genomes
}

// genomes a target in several genomes is placed in first for the PMO, the variant calling genomes first so the target locations match the calls
def pmoPreferredGenomes() {
    def genomes = variantCallGenomes()
    commaListParam(params.pmo_preferred_genomes).each { genome ->
        if (!genomes.contains(genome)) {
            genomes << genome
        }
    }
    return genomes.join(",")
}

// target -> genome from ASSIGN_TARGET_GENOMES' targetGenomes.tsv
def readTargetGenomes(assigned_dir) {
    return file("${assigned_dir}/targetGenomes.tsv").readLines().drop(1).findAll { line -> line.trim() }.collectEntries { line ->
        def cols = line.split("\t")
        [(cols[0]): cols[1]]
    }
}

// the genome a clustered target is variant called against, targets clustered as separate mates (<target>-R1, <target>-R2) use <target>'s
def variantCallGenomeFor(target_genomes, target) {
    if (target_genomes.containsKey(target)) {
        return target_genomes[target]
    }
    if (target.endsWith("-R1") || target.endsWith("-R2")) {
        return target_genomes[target.substring(0, target.length() - 3)] ?: ""
    }
    return ""
}
