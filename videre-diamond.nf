#!/usr/bin/env  nextflow

params.pep_ref   = "/opt/DB_REF/mmetsp_pep/MMETSP_test.pep.fa"
params.output    = "${PWD}/Diamond"
params.DB_REF    = System.getenv('DB_REF')
params.taxonmap  = "/opt/DB_REF/mmetsp_pep/seq_taxid.map2"
params.taxonnodes = "/opt/DB_REF/taxonomy/nodes.dmp" 
params.queries_path = "Videre.Out/MegaHit/"
params.diamond_idx = true
params.diamond   =  true
diamond_raw      =  file(params.pep_ref)
query_seq        =  file(params.queries_path)
output           =  params.output


Channel.fromPath(params.queries_path +'/*')
       .into{dmnd_queries; dmnd_queries1 }



log.info"""

Diamond ref     = ${diamond_raw}
Ouput folder    = ${output}
HTP cores       = ${params.htp_cores} 
MTP cores       = ${params.mtp_cores} 
LTP cores       = ${params.ltp_cores} 
DB Ref          = ${params.DB_REF}
Build DB        = ${params.diamond_idx}
RUN Diamond     = ${params.diamond}


Queries
======
"""
dmnd_queries1.subscribe{ println it }

log.info"""


"""

process diamond_idx{

    cpus params.htp_cores
    memory params.h_mem
    echo true
    //storeDir "$params.DB_REF/Diamond"
    input:
        file diamond_raw

    output:
	file "${dmnd_base}.dmnd" into diamond_idx
        val dmnd_base into diamond_BaseName
	
    when:
        params.diamond_idx == true

    script:
        dmnd_base = diamond_raw.baseName


"""
    
    diamond \
    makedb \
    --in ${diamond_raw} \
    -d ${dmnd_base} \
    --threads ${params.htp_cores} \
    --taxonmap ${params.taxonmap} \
    --taxonnodes ${params.taxonnodes} \
    -v
     
"""
 
}



if  ( params.diamond_idx == false  ){

    dmnd_base = diamond_raw.baseName
    diamond_idx = Channel.fromPath("${params.DB_REF}/Diamond/${dmnd_base}.dmnd")
    diamond_BaseName = Channel.value(dmnd_base)

}



process diamond{

    echo true
    cpus params.htp_cores
    memory params.h_mem
    publishDir path: output , mode: 'move'
    
    input:
       file query_seq from dmnd_queries
       file diamond_db from  diamond_idx
       val  dmnd_base from diamond_BaseName
       
    output:
	file(ref_tag) into diamond_matches
    
    when:
      params.diamond == true

    script:
     ref_tag  =  query_seq.getName()+"_dmnd"

       	
"""

   mkdir -v ${ref_tag} 

    diamond \
    blastx \
    -d ${dmnd_base}  \
    --un ${ref_tag}/diamond.unaligned \
    --al ${ref_tag}/diamond.aligned \
    -q ${query_seq} \
    -o ${ref_tag}/matches.dmnd \
    --more-sensitive \
    --evalue 1e-5 \
    --top 90 \
    --id 40 \
    --block-size 5 \
    --index-chunks 1 \
    -v    

"""
    
}


