dcd#!/usr/bin/env  nextflow


//params.pep_ref 		= "/opt/DB_REF/mmetsp_pep/MMETSP_test.pep.fa"
params.pep_ref 		= "/home/andhlovu/DB_REF/mmetsp_pep/Combined_MMETSP.pep.fa"
//params.nt_ref 		= "/opt/DB_REF/mmetsp_nt/MMETSP_test.nt.fa"
params.nt_ref           = "/home/andhlovu/DB_REF/mmetsp_nt/Combined_MMETSP.nt.fa"
params.output 		= "${PWD}/Diamond"
params.DB_REF 		= System.getenv('DB_REF')
params.taxanodes 	= "/opt/DB_REF/taxonomy/nodes.dmp" 
params.queries_path     =  "/home/drewx/Documents/videre-pipeline/query"
params.diamond_idx 	= true
params.diamond   	= false
params.makeblastdb      = false
params.megablast        = false
diamond_raw       	= file(params.pep_ref)
output           	= params.output
blastdb_raw             = file(params.nt_ref)



Channel.fromPath(params.queries_path +'/*')
       .into{queries ;diamond_queries; blast_queries}




log.info"""

Diamond ref     = ${diamond_raw}
Queries path	= ${params.queries_path}
blastdb_raw     = ${blastdb_raw}
Ouput folder    = ${output}
HTP cores       = ${params.htp_cores} 
MTP cores       = ${params.mtp_cores} 
LTP cores       = ${params.ltp_cores} 
DB Ref          = ${params.DB_REF}
Build dmnd IDX  = ${params.diamond_idx}
RUN Diamond     = ${params.diamond}
Make BlastDB    = ${params.makeblastdb}      

Queries
======
"""

queries.subscribe{ println it }

log.info"""


"""


process diamond_idx{

    echo true
    cpus params.htp_cores
    memory params.h_mem
    storeDir "$params.DB_REF/Diamond"
    
    
    input:
       file diamond_raw

    output:
       file "${dmnd_base}.dmnd" into diamond_idx
       val dmnd_base into diamond_BaseName
	
    when: params.diamond_idx == true

    script:
        dmnd_base = diamond_raw.baseName

"""
    
    diamond \
    makedb \
    --in ${diamond_raw} \
    -d ${dmnd_base} \
    --threads ${params.htp_cores} \
    -v
     
"""
 
}




process makeblastdb{

    echo true
    cpus params.htp_cores
    memory params.h_mem
    storeDir "$params.DB_REF"

    input:
        file blastdb_raw 

    output:
       file('Blast') into blastdb_idx
       val  blastdb into blastdb_BaseName

    when:
        params.makeblastdb == true

    script:
        blastdb = blastdb_raw.baseName
    

"""

    makeblastdb \
    -in ${blastdb_raw} \
    -input_type fasta \
    -dbtype nucl \
    -out Blast/${blastdb} \
    -parse_seqids 
     
    makembindex \
    -input Blast/${blastdb} \
    -iformat blastdb \
    -old_style_index false


"""


}




if  ( params.diamond_idx == false  ){

    dmnd_base = diamond_raw.baseName
    diamond_idx = Channel.fromPath("${params.DB_REF}/Diamond/${dmnd_base}.dmnd")
    diamond_BaseName = Channel.value(dmnd_base)

}



process diamond{

    cpus params.htp_cores
    memory params.h_mem
    publishDir path: output , mode: 'copy'
    
    input:
       file query_seqs from diamond_queries
       file diamond_db from  diamond_idx
       val  dmnd_base from diamond_BaseName
       
    output:
	file(diamond_tag) into diamond_Out
    
    when:
      params.diamond == true

    script:
       diamond_tag  =  query_seqs.getName()+"_dmnd"

       	
"""

   mkdir -v ${ref_tag} 

    diamond \
    blastx \
    -d ${dmnd_base}  \
    --un ${ref_tag}/diamond.unaligned \
    --al ${ref_tag}/diamond.aligned \
    -q ${query_seqs} \
    -o ${diamond_tag}_diamond.out \
    -f 5  \
    --more-sensitive \
    --id 40 \
    --header \
    --max-hsps 1\
    --evalue 1e-5 \
    --block-size 5 \
    --index-chunks 1 \
    --verbose    

"""

    
}



process MegaBlast{
    
    echo true
    cpus params.htp_cores
    memory params.h_mem
    publishDir path: output , mode: 'copy'
    
    input:
       file query_seqs from blast_queries
       file blastdb_file from  blastdb_idx
       val  blastdb_name from blastdb_BaseName
       
    output:
	file('Megablast') into blastOut
    
    when:
        params.megablast == true

    script:
        megablast_tag  =  query_seqs.getName()+"_blast"

       	

"""
 
   blastdbcmd \
   -db Blast/${blastdb_name}  \
   -info

   mkdir -v Megablast
   blastn \
   -query ${query_seqs} \
   -task megablast \
   -db ${params.DB_REF}/Blast/${blastdb_name} \
   -num_threads  ${params.htp_cores} \
   -outfmt 5 \
   -evalue 1e-5 \
   -out Megablast/${megablast_tag}_blast.out \
   -parse_deflines \
   -use_index true
     
"""

    
}
