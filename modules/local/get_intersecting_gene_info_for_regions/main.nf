process GET_INTERSECTING_GENE_INFO_FOR_REGIONS {
    label 'process_single'

    publishDir "${publish_dir}", mode: 'copy', overwrite: true

    input:
    path bed_fnp
    path gff_fnp
    path genome_twobit_fnp
    val publish_dir

    output:
    path "bed_withGeneInfo.bed", emit: bed_withGeneInfo_bed
    path "bed_withGeneInfo.tsv", emit: bed_withGeneInfo_tsv
    path "bed_withGeneAAInfo.bed", emit: bed_withGeneAAInfo_bed

    script:

    """
    cut -f1-6 ${bed_fnp} | elucidator bedGetIntersectingGenesInGff --extraAttributes description,Name --gff ${gff_fnp} --overWrite --bed STDIN --out bed_withGeneInfo.bed
    elucidator splitColumnContainingMeta --column col.6 --delim tab --removeEmptyColumn --addHeader --replacementHeader "#chrom,start,end,name,length,strand" --overWrite  --file bed_withGeneInfo.bed --out bed_withGeneInfo.tsv

    cut -f1-6 ${bed_fnp} | elucidator bedGetOverlappingAminoAcidPositions  --bed STDIN --gff /tank/data/genomes/plasmodium/genomes/pf/info/gff/Pf3D7.gff --genomeTwoBit /tank/data/genomes/plasmodium/genomes/pf/genomes/Pf3D7.2bit --extraAttributes Name,description --out bed_withGeneAAInfo.bed

    """
}

