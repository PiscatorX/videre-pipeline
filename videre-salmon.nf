#! /usr/bin/env nextflow

//params.readsbase    = "/home/drewx/Documents/subsample"
params.readsbase    = "/home/andhlovu/Novogene/ftpdata.novogene.cn:2300/C101HW18111065/raw_data"
params.pe_patt      = "*_RNA_{1,2}.fq.gz" 
params.DB_REF 	    = System.getenv('DB_REF')
<<<<<<< Updated upstream
params.queries_path = "Contigs"
query_seq           =  file(params.queries_path)
params.output       = "${PWD}/Salmon"
params.cdHit_perc   = 0.98
output              =  params.output
=======
//params.queries_path = "Contigs"
params.queries_path = "/home/andhlovu/Metatranscriptomics_DevOps/megahit_contig/MegaHit/MegaHit.fasta"
query_seq        =  file(params.queries_path)
params.output    = "${PWD}/Salmon"
params.cdHit_perc       = 0.98
output           =  params.output
>>>>>>> Stashed changes
params.salmon_index = true
params.salmon_quant = true
params.gmst         = true


Channel.fromPath(params.queries_path +'/*')
    .into{contig_queries; contig_queries1}

reads = params.readsbase +'/'+ params.pe_patt

Channel.fromFilePairs(reads)
      .ifEmpty{ error "Could not locate pair reads: ${reads}"}
      .set{get_reads}


get_reads.into{reads1; reads2; reads3; readx; reads_cp}



log.info """

Read file pattern 	= ${reads}
Contigs_path            = ${params.queries_path}    
Output			= ${output}
High TP cores    	= ${params.htp_cores}
Midium TP cores    	= ${params.mtp_cores} 
Low TP cores    	= ${params.ltp_cores}
H_mem  			= ${params.h_mem}
Salmon Index            = ${params.salmon_index}
Salmon Quant            = ${params.salmon_quant}
GMST                    = ${params.gmst} 


Reads
=====
"""
readx.each{  if(it instanceof List){println it} }

log.info"""
---------------------------------------------------------
"""




process cd_hit_est{
    
    echo true
    cpus params.htp_cores
    memory "${params.l_mem} GB"
    publishDir path: "${output}/CD-Hit", mode: 'copy'
  
   input:
       file(contigs) from contig_queries
       
    
   output:
       file "${hits_base}.cd_hits" into cd_hits
       file "${hits_base}.cd_hits.clstr" into cdhit_clusters
       val  hits_base into cd_hits_base
       
    script:
	hits_base = "${contigs.baseName}"

"""
    cd-hit-est \
    -i $contigs \
    -c ${params.cdHit_perc} \
    -T ${params.htp_cores} \
    -M 0 \
    -d 0 \
    -r 0 \
    -p 1 \
    -g 1 \
    -o ${hits_base}.cd_hits 
    
"""

}


process contig_initDB{

    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${output}/CD-Hit", mode: 'copy'
    
    input:
	file(cd_hits) from cd_hits
        val hits_base from cd_hits_base
	
    output:	
        file(outfile) into (cd_hits1, cd_hits2, cd_hits3)
        file(tsv_file) into fasta_ref
        file(sqlite_db) into (sqlite_db1, sqlite_db2)

    script:
	sqlite_db  = "${hits_base}X.db"
        outfile    = "${hits_base}X.fasta"
        tsv_file   = "${hits_base}X.db"
   
"""
	
   contig_initDB.py \
   -d ${sqlite_db} \
   -c ${cd_hits} \
   -f\
   -o ${outfile}
   
"""
   
}


   

process salmon_index{
    
    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    storeDir "${params.DB_REF}/Salmon"
    
    input:
	file(cd_hits) from cd_hits1
    
    output:
        file("salmon_index") into salmon_index

    when:
	params.salmon_index == true
	    
"""

    salmon \
    --no-version-check \
    index \
    -t  ${cd_hits}  \
    -i  salmon_index \
    --type quasi \
    -p  ${params.htp_cores} 

           
"""
    
//https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon
//https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-mapping-based-modex
   
}




process salmon_quant{
    
    //echo  true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: output, mode: 'move'
    
    input:
	file(index) from salmon_index
        set pair_id, file(reads) from reads3
      
    output:
        file("salmon_${pair_id}") into salmon_quant

    when:
	params.salmon_quant == true
    
    script:
    	(left, right)=reads
    
"""

    salmon \
    quant \
    --no-version-check \
    -l iu \
    -1 $left \
    -2 $right \
    --index salmon_index \
    --minAssignedFrags 0 \
    --validateMappings \
    --writeUnmappedNames \
    --meta \
    --output salmon_${pair_id} \
    --discardOrphansQuasi \
    -p ${params.htp_cores}


"""
	
}


    
process GeneMarkST{
     
    //echo  true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: "$output/GeneMarkST", mode: 'move'
    input:
	file(cd_hits) from cd_hits2

    output:
	file("${cd_hits}*") into gmst_out
	file("gms.log") into gms_log

    when:
	params.gmst == true
    
"""

      gmst.pl \
      --fnn \
      --faa \
      --format GFF \
      ${cd_hits} \
      --verbose
      
         
"""
//output sent to GhostKoala

}
