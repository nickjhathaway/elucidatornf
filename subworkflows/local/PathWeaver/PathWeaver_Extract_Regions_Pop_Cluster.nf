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


    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(results_dir)
    if (results_dir_obj.exists()){
        results_dir_obj.deleteDir()
    }
    results_dir_obj.mkdirs()

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
    def output = file("${params.pw_results_dir}/run/parameters.tsv")
    output.withWriter { writer -> 
        params.each { k, v -> 
            writer.println("${k}\t${v}")} 
    }
}

/* Record runtime information
 *
 * Records runtime and environment information and writes summary to a tabulated (tsv) file
 *
 */
def record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime() {

    def output = file("${params.pw_results_dir}/run/runtime.tsv")
    output.withWriter { writer -> 
        writer.println("PipelineVersion\t${workflow.manifest.version}")
        writer.println("ContainerEngine\t${workflow.containerEngine}")
        writer.println("Duration\t${workflow.duration}")
        writer.println("CommandLine\t${workflow.commandLine}")
        writer.println("CommitId\t${workflow.commitId}")
        writer.println("Complete\t${workflow.complete}")
        writer.println("ConfigFiles\t${workflow.configFiles.join(', ')}")
        writer.println("Container\t${workflow.container}")
        writer.println("ErrorMessage\t${workflow.errorMessage}")
        writer.println("ErrorReport\t${workflow.errorReport}")
        writer.println("ExitStatus\t${workflow.exitStatus}")
        writer.println("HomeDir\t${workflow.homeDir}")
        writer.println("LaunchDir\t${workflow.launchDir}")
        writer.println("Manifest\t${workflow.manifest}")
        writer.println("Profile\t${workflow.profile}")
        writer.println("ProjectDir\t${workflow.projectDir}")
        writer.println("Repository\t${workflow.repository}")
        writer.println("Resume\t${workflow.resume}")
        writer.println("Revision\t${workflow.revision}")
        writer.println("RunName\t${workflow.runName}")
        writer.println("ScriptFile\t${workflow.scriptFile}")
        writer.println("ScriptId\t${workflow.scriptId}")
        writer.println("ScriptName\t${workflow.scriptName}")
        writer.println("SessionId\t${workflow.sessionId}")
        writer.println("Start\t${workflow.start}")
        writer.println("StubRun\t${workflow.stubRun}")
        writer.println("Success\t${workflow.success}")
        writer.println("UserName\t${workflow.userName}")
        writer.println("WorkDir\t${workflow.workDir}")
        writer.println("NextflowBuild\t${nextflow.build}")
        writer.println("NextflowTimestamp\t${nextflow.timestamp}")
        writer.println("NextflowVersion\t${nextflow.version}")
    }
}