// main.nf

process ALIGN_AND_SORT {
    // This process will automatically use minimap2.sif and 24GB RAM
    input:
    path reads
    path ref_genome
    path ref_genome_fai

    output:
    path "aligned.bam"

    script:
    """
    minimap2 -a $ref_genome $reads | samtools sort -o aligned.bam -
    """
}

process INDEX_BAM {
    // This process will automatically use minimap2.sif and 8GB RAM from nextflow.config
    input:
    path bam_file

    output:
    // Output both the bam and the index as a pair so Clair3 receives them together
    tuple path(bam_file), path("${bam_file}.bai")

    script:
    """
    samtools index $bam_file
    """
}

process CALL_VARIANTS {
    // This process will automatically use clair3.sif and 24GB RAM
    input:
    // Accept the tuple containing both the BAM and its index
    tuple path(bam_file), path(bai_file)
    path ref_genome
    path ref_genome_fai

    output:
    path "output.vcf"

    script:
    """
    # Clair3 command goes here
    run_clair3.sh --bam_fn=$bam_file --ref_fn=$ref_genome --output=. --threads=${task.cpus} --platform="ont"
    """
}

workflow {
    reads_ch = Channel.fromPath(params.reads)
    ref_ch = Channel.fromPath(params.ref_genome)
    ref_fai_ch = Channel.fromPath(params.ref_genome_fai)

    // 1. Align and sort (with indexed reference for faster processing)
    bam_ch = ALIGN_AND_SORT(reads_ch, ref_ch, ref_fai_ch)
    
    // 2. Index the BAM
    bam_and_index_ch = INDEX_BAM(bam_ch)
    
    // 3. Call variants (passing the paired BAM, index, and reference index channel)
    CALL_VARIANTS(bam_and_index_ch, ref_ch, ref_fai_ch)
}
