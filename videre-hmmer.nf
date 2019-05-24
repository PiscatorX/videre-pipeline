#!/usr/bin/env nextflow

params.peptides = "Salmon/GeneMarkST/*.faa"
params.output	= "hmmer"
params.DB_REF 	= System.getenv('DB_REF')
DB_REF		= params.DB_REF
pfam_hmmerDB    = Channel.value("$DB_REF/Pfam/Pfam-A.hmm")	
output     	= params.output

Channel.fromPath(params.peptides)
    .ifEmpty{ error "Could not locate pair contigs files => ${params.queries_path}" }
    .set{peptide_queries}



process hmmscan{

    echo true
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${output}/hmmer", mode: 'copy'

    input:
	file(pep_file) from peptide_queries
	val pfam_hmmerDB

    output:
        file("*") into hmmscan_out

    script:
	 hmmer_base = pep_file.baseName


"""
 
   hmmscan \
   --tblout ${hmmer_base}.tbl \
   -o ${hmmer_base}.hmmscan \
   -E 100 \
   --cpu ${params.htp_cores} \
   --max \
   --domtblout ${hmmer_base}.domtbl \
   --pfamtblout ${hmmer_base}.pfamtbl \
   ${pfam_hmmerDB} \
   ${pep_file} 
  
"""

}

