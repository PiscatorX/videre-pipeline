#!/bin/env nextflow

params.reads = "$baseDir/data1/*_{1,2}.fq"
output = "$PWD/Videre.O"


Channel.fromFilePairs(params.reads)
       .ifEmpty{ error "Could not locate reads: ${params.reads}"}
       .set{read_pairs}


process fastqc{
    echo true
    publishDir path: output, mode: 'copy'
    input: set pair_id, file(reads) from read_pairs
    


"""
    mkdir -v $pair_id
    ls  -la
    fastqc --extract -f fastq -o $pair_id -t 1 *.fastq

"""

}



// process trimmomatic{
//     echo true
//     publishDir path: output, mode: 'copy'
   
//        input: set pair_id, file(reads) from read_pairs


// """

