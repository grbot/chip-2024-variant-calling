#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC {
    publishDir "${params.output}/", mode: 'copy', overwrite: false
    input: 
    tuple val(id), path(reads)

    output:
    path("fastqc/*")
    
    script:
    """
    mkdir fastqc
    fastqc -t 2 -o fastqc $reads
    """
}

process TRIMGALORE {
    publishDir "${params.output}/", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*_val_*.fq.gz"), emit: trimmed_fq

    script:
    """
    trim_galore --paired --illumina -j 4 $reads --gzip
    """
}

process BWA {
    publishDir "${params.output}/", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*.sam"), emit: align_sam

    script:
    """
    bwa mem -t 4 $params.ref $reads > ${id}.sam
    """
}

process SAMTOBAM {
    publishDir "${params.output}/", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path( "*.sorted.bam"), emit: align_bam
    
    script:
    """
    samtools view -1 $sam -o ${id}.bam 
    samtools sort -m 11G -@ 4 ${id}.bam -o ${id}.sorted.bam
    samtools index ${id}.sorted.bam
    """
}

process MARKDUP {
    publishDir "${params.output}/", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(sorted_bam)

    output:
    tuple val(id), path("*.markdup.bam"), emit: markdup_bam

    script:
    """
    picard MarkDuplicates I=$sorted_bam O=${id}.markdup.bam M=${id}.metrics.txt  
    """
}

process PILEUP {
    publishDir "${params.output}/", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(markdup_bam)

    output:
    tuple val(id), path("${id}.bcf"), emit: call_bcf

    script:
    """
    bcftools mpileup -Ou -f $params.ref $markdup_bam | bcftools call -mv -Ob -o ${id}.bcf
    """
}

process CALLVAR {
    publishDir "${params.output}/", mode: 'copy', overwrite: true

    input:
    tuple val(id), path (call_bcf) 

    output:
    path "*.vcf"

    script:
    """
    bcftools view -i 'QUAL>100 && DP>=50' -v indels $call_bcf -Ov -o ${id}.vcf
    """
}


workflow {
    def fastq = Channel.fromFilePairs(params.fastqs)

    FASTQC(fastq)
    TRIMGALORE(fastq)
    BWA(TRIMGALORE.out.trimmed_fq)
    SAMTOBAM(BWA.out.align_sam)
    MARKDUP(SAMTOBAM.out.align_bam)
    PILEUP(MARKDUP.out.markdup_bam)
    CALLVAR(PILEUP.out.call_bcf)
}
