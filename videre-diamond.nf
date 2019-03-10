#!/usr/bin/env  nextflow

params.pep_ref   = "diamond/MMETSP1449.trinity_out_2.2.0.Trinity.pep.fasta"
//params.output    = "$DB_REF/Diamond"
params.output    = "${PWD}/Diamond"
params.DB_REF    = System.getenv('DB_REF')
params.query_seq = "Videre.Out/MegaHit/MegaHit.fasta.1"


diamond_raw      =  file(params.pep_ref)
query_seq        =  file(params.query_seq)
output           =  params.output   


log.info"""
Diamond ref     = ${diamond_raw}
Ouput folder    = ${output}
HTP cores       = ${params.htp_cores} 
MTP cores       = ${params.mtp_cores} 
LTP cores       = ${params.ltp_cores} 
DB Ref          = ${params.DB_REF}

"""

process diamond_idx{

    echo true
    storeDir "$params.DB_REF/Diamond"
    input:
        file diamond_raw

    output:
	file "${dmnd_base}.dmnd" into diamond_idx
        val dmnd_base into diamond_BaseName
    
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




process diamond{

    publishDir path: output, mode: 'move'
    
    input:
       file query_seq
       file diamond_db from  diamond_idx
       val  dmnd_base from diamond_BaseName
    
    output:
	file('matches.dmnd') into diamond_matches
	file('diamond.unaligned') into diamond_un
	file('diamond.aligned') into diamond_al    

    
"""

    diamond \
    blastx \
    -d ${dmnd_base}  \
    --un diamond.unaligned \
    --al diamond.aligned \
    -q ${query_seq} \
    -o matches.dmnd \
    --more-sensitive \
    --evalue 1e-5 \
    --top 90 \
    --id 40 \
    -v    

"""
    
}



diamond_idx.subscribe{ println it }
