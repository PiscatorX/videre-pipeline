#! /usr/bin/env nextflow

params.bamdir = "/home/drewx/Documents/videre-pipeline/Bowtie2sam"
params.patt 	="*_trim_?P.bam"
params.output   = "${PWD}/Anvio"
params.contigs  = "/home/drewx/Documents/videre-pipeline/Contigs/MegaHit.fasta"
output = params.output 
contig = file(params.contigs)


reads = params.bamdir +'/'+ params.patt  
Channel.fromPath(reads)
      .ifEmpty{ error "Could not locate pair reads: ${reads}"}
      .set{bam_files}



process reformat_contig{
    
    label 'anvio'
    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: "$output/multiqc_RawReads", mode: 'copy'
    
    input:
        file contig

    output:
	set "anvio_${contig}", "${contig.baseName}.report" into Anvio

"""
    anvi-script-reformat-fasta \
    ${contig} \
    -o anvio_${contig} \
    -l 0 \
    --simplify-names \
    --report-file ${contig.baseName}.report
  
"""
    
}




process init_bam{

    label  'anvio'
    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: output, mode: 'move'
    
    
    input:
	file bam from bam_files
        file contig
    output:
	set "*sorted.bam","*bam-sorted.bam.bai" into anvio
        
"""
    anvi-init-bam \
    ${bam} 

"""
    
}
