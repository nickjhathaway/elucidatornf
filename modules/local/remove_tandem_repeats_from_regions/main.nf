process REMOVE_TANDEM_REPEATS_FROM_REGIONS {
    label 'process_single'


    input:
    path bed_fnp
    path genome_twobit_fnp


    output:
    path "out_with_tandems_removed.bed", emit: out_with_tandems_removed

    script:

    """
    elucidator getFastaWithBed --twoBit $genome_twobit_fnp --bed $bed_fnp --out out.fasta
    elucidator runTRF --fasta out.fasta --genomicLocation $bed_fnp --dout trfOutputput --supplement
    cat trfOutputput/genomic_combined.bed | \
            bedtools merge -d 10 | \
            elucidator bed3ToBed6 --bed STDIN | \
            elucidator removeSubRegionsFromBedFile --bed out_allTranscripts.bed --removeBed STDIN  | \
            elucidator filterBedRecordsByLength --bed STDIN  --minLen 30 > out_with_tandems_removed.bed
    """
}

