process GET_GENE_RECORDS_FOR_GENE_IDS {

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "infos"

    label 'process_single'

    input:
    path ids_fnp
    path genome_twobit_fnp
    path gff_fnp
    val pub_dir

    output:
    path "infos", emit: infos
    path "infos/out_allTranscripts.bed", emit: out_allTranscripts_bed
    path "infos/out_allTranscripts.tsv", emit: out_allTranscripts_tsv

    script:

    """
    elucidator gffRecordIDToGeneInfo --gff ${gff_fnp} --2bit ${genome_twobit_fnp}  --overWrite --id $ids_fnp --dout infos --overWriteDir
    elucidator splitColumnContainingMeta --file infos/out_allTranscripts.bed  --column col.6 --delim tab --removeEmptyColumn --addHeader --replacementHeader "#chrom,start,end,name,length,strand" --overWrite --out infos/out_allTranscripts.tsv
    """
}

