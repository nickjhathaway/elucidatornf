process SALMON_QUANT_COMBINE {

    publishDir "${pubdir}", mode: 'copy', overwrite: true

    tag "salmon_quant_combine"

    input:
    tuple path(quant_files), val(pubdir)

    output:
    path("all_quants.tsv.gz")

    script:
    """
    elucidator rBind \
        --contains _salmon_quants_quant.sf.gz \
        --delim tab \
        --header \
        --overWrite \
        --out all_quants.tsv.gz
    """
}

