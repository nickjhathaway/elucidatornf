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

    workflow.onComplete {
        def outputDir = new File("${results_dir}/run")
        if (!outputDir.exists()) {
            outputDir.mkdirs()
        }
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_params()
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime()
    }
}


  
def record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_params() {
    def output = new File("${params.pw_results_dir}/run/parameters.tsv")
    params.each{ k, v -> 
        output.append("${k}\t${v}\n")
    }
}

/* Record runtime information
 *
 * Records runtime and environment information and writes summary to a tabulated (tsv) file
 *
 */
def record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime() {

    def output = new File("${params.pw_results_dir}/run/runtime.tsv")

    // Append ContainerEngine first as shown in your example
    output.append("PipelineVersion\t${workflow.manifest.version}\n")
    output.append("ContainerEngine\t${workflow.containerEngine}\n")
    output.append("Duration\t${workflow.duration}\n")
    output.append("CommandLine\t${workflow.commandLine}\n")
    output.append("CommitId\t${workflow.commitId}\n")
    output.append("Complete\t${workflow.complete}\n")
    output.append("ConfigFiles\t${workflow.configFiles.join(', ')}\n")
    output.append("Container\t${workflow.container}\n")
    output.append("ErrorMessage\t${workflow.errorMessage}\n")
    output.append("ErrorReport\t${workflow.errorReport}\n")
    output.append("ExitStatus\t${workflow.exitStatus}\n")
    output.append("HomeDir\t${workflow.homeDir}\n")
    output.append("LaunchDir\t${workflow.launchDir}\n")
    output.append("Manifest\t${workflow.manifest}\n")
    output.append("Profile\t${workflow.profile}\n")
    output.append("ProjectDir\t${workflow.projectDir}\n")
    output.append("Repository\t${workflow.repository}\n")
    output.append("Resume\t${workflow.resume}\n")
    output.append("Revision\t${workflow.revision}\n")
    output.append("RunName\t${workflow.runName}\n")
    output.append("ScriptFile\t${workflow.scriptFile}\n")
    output.append("ScriptId\t${workflow.scriptId}\n")
    output.append("ScriptName\t${workflow.scriptName}\n")
    output.append("SessionId\t${workflow.sessionId}\n")
    output.append("Start\t${workflow.start}\n")
    output.append("StubRun\t${workflow.stubRun}\n")
    output.append("Success\t${workflow.success}\n")
    output.append("UserName\t${workflow.userName}\n")
    output.append("WorkDir\t${workflow.workDir}\n")
    output.append("NextflowBuild\t${nextflow.build}\n")
    output.append("NextflowTimestamp\t${nextflow.timestamp}\n")
    output.append("NextflowVersion\t${nextflow.version}\n")
}