process DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE {
    label 'process_low'

    cpus { ncpus }

    publishDir "${pub_results_dir}", mode: 'copy', overwrite: true, pattern: "targets.bed"
    publishDir "${pub_results_dir}", mode: 'copy', overwrite: true, pattern: "seqs.bed"


    input:
    path seq_table
    path genome_fnp
    val seq_col_name
    val name_col_name
    val target_col_name
    val ncpus
    val pub_results_dir


    output:
    path "seqs.bed", emit: seqs_locs
    path "targets.bed", emit: targets_bed

    script:

    """
    elucidator tableExtractColumns \
            --file allSelectedClustersInfo.tab.txt.gz \
            --delim tab --header \
            --columns ${target_col_name},${name_col_name},${seq_col_name} \
            --getUniqueRows | elucidator createSeqsFromTable \
                                    --nameColumn ${name_col_name} \
                                    --seqColumn ${seq_col_name} \
                                    --file STDIN \
                                    --overWrite \
                                    --out seqs.fasta.gz

    elucidator splitSeqFileWithMeta --fastagz seqs.fasta.gz --metaField ${target_col_name}

    rm -f indvidualBedFnps.txt
    rm -f determineRegionCmds.txt
    rm -f bedFnps.txt
    for x in `elucidator printCol --file allSelectedClustersInfo.tab.txt.gz --delim tab --header --columnName ${target_col_name} --unique --sort`; do
        echo elucidator determineRegionLastz --fasta \${x}.fasta.gz \
                --genome ${genome_fnp} \
                --out  \${x}.bed\
                --name \${x} \
                --keepBestOnly \
                --individualOut \${x}_individual.bed >> determineRegionCmds.txt
        echo \${x}_individual.bed >> indvidualBedFnps.txt
        echo \${x}.bed >> bedFnps.txt
    done;
    elucidator runMultipleCommands --cmdFile determineRegionCmds.txt --numThreads ${ncpus} --raw

    elucidator rBind --skipNonExistFiles --files indvidualBedFnps.txt --delim tab  | elucidator bedCoordSort --bed STDIN --out seqs.bed --overWrite
    elucidator rBind --skipNonExistFiles --files bedFnps.txt --delim tab  | elucidator bedCoordSort --bed STDIN --out targets.bed --overWrite
    """
}


