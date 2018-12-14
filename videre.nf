#!/bin/env nextflow

params.reads = "SbaseDir/data/"


Channel.
	.fromFilePairs(params.reads)
	.ifEmpty{ error "Could not locate reads: ${params.reads}"}
	.set{read_pairs}


process trimmomatic{
    input: 
}
