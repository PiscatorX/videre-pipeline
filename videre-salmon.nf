#! /usr/bin/env nextflow


//params.readsbase    = "/home/drewx/Documents/subsample"
params.readsbase    = "/home/andhlovu/Novogene/ftpdata.novogene.cn:2300/C101HW18111065/raw_data"
params.pe_patt      = "*_trim_{1,2}P.fastq"
params.DB_REF 	    = System.getenv('DB_REF')
query_seq           =  file(params.queries_path)
params.output       = "${PWD}/Salmon"
params.cdHit_perc   = 0.98
output              =  params.output
//params.queries_path = "Contigs"
params.queries_path = "/home/andhlovu/Metatranscriptomics_DevOps/megahit_contig/MegaHit/MegaHit.fasta"
query_seq           =  file(params.queries_path)
output              = params.output
DB_REF		    = params.DB_REF
params.bowtie_idx   = true
params.bowtie       = true
params.salmon_index = false
params.salmon_quant = false
params.gmst         = true


Channel.fromPath(params.queries_path +'/*')
    .ifEmpty{ error "Could not locate pair contigs files => ${params.queries_path}" }
    .set{contig_queries}


reads = params.readsbase +'/'+ params.pe_patt

Channel.fromFilePairs(reads)
       .ifEmpty{ error "Could not locate pair reads: ${reads}"}
    .into{reads1; reads2; reads3; reads4; readx}




log.info """

Read file pattern 	= ${reads}
Contigs_path            = ${params.queries_path}    
Output			= ${output}
High TP cores    	= ${params.htp_cores}
Midium TP cores    	= ${params.mtp_cores} 
Low TP cores    	= ${params.ltp_cores}
H_mem  			= ${params.h_mem}
Bowti IDX               = ${params.bowtie_idx}
Bowtie                  = ${params.bowtie}
Salmon Index            = ${params.salmon_index}
Salmon Quant            = ${params.salmon_quant}
GMST                    = ${params.gmst} 


Reads
=====
"""

readx.subscribe{  if(it instanceof List){ println it} }


log.info"""
---------------------------------------------------------
"""


process cd_hit_est{
   
    //echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
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
        file(outfile) into (cd_hits_bowtie, cd_hits_salmon, cd_hits_gmst)
        file(tsv_file) into fasta_ref
        file(sqlite_db) into sqlite_db1

    script:
	sqlite_db  = "${hits_base}X.db"
        outfile    = "${hits_base}X.fasta"
        tsv_file   = "${hits_base}X.tsv"
   
"""
	
   contig_initDB.py \
   -d ${sqlite_db} \
   -c ${cd_hits} \
   -f \
   -o ${outfile}
      
"""
   
}




process bowtie_idx{

    //echo true
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    //storeDir "${DB_REF}/Bowtie"
    publishDir "${DB_REF}/Bowtie"
    
    input:
        file contig_fasta from cd_hits_bowtie
    
    output:  
        set bowtie2_base, file("${bowtie2_base}*") into bowtie_idx
       
        
    when:
	params.bowtie_idx == true
    
    script:
	bowtie2_base =  "bowtie2_${contig_fasta}".replaceFirst(/.fasta/, "")
      
    
"""
    
    bowtie2-build  \
    --large-index \
    --threads ${params.htp_cores} \
    ${contig_fasta} \
    ${bowtie2_base}
    bowtie2-inspect \
    --large-in \
    --summary \
    ${bowtie2_base}  >  bowtie2_${contig_fasta}.idx_stats

    
"""
    
    
}



process bowtie2bam{

    echo true
    tag "${sample}"
    cpus params.htp_cores
    //storeDir "${output}/Bowtie2sam"
    publishDir "${DB_REF}/Bowtie2sam", mode: "copyNoFollow"
    memory "${params.h_mem} GB"

    input:
        set sample, file(reads) from reads3
        set bowtie2_base, file(bowties_idx_files) from bowtie_idx.collect()
         
    output:
	file("${sample}*") into bowtie_sam
        file("${sample}.un")   into bowtie_unaligned
    
    when:
	params.bowtie == true

    script:
        (fwd_reads, rev_reads) = reads 
	fwd_name = fwd_reads.baseName
        rev_name = rev_reads.baseName
    
"""

   bowtie2 \
    --threads ${params.htp_cores} \
    -x ${bowtie2_base} \
    -1 ${fwd_reads} \
    -2 ${rev_reads} \
    --no-unal\
    --time \
    --un ${sample}.un \
    -S ${sample}.sam


    samtools \
    view \
    ${sample}.sam \
    -F 4 \
    -b \
    -o ${sample}.bam 
   
   
 
"""
//     o mySampleFiltered_1P.fq.gz - for paired forward reads
//     o mySampleFiltered_1U.fq.gz - for unpaired forward reads
//     o mySampleFiltered_2P.fq.gz - for paired reverse reads
//     o mySampleFiltered_2U.fq.gz - for unpaired forward reads
//--un <path>        write unpaired reads that didn't align to <path>
// --no-unal          suppress SAM records for unaligned reads
//http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

}




process salmon_index{
    
    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    //storeDir "${params.DB_REF}/Salmon"
    publishDir "${DB_REF}/Salmon"

    input:
	file(cd_hits) from cd_hits_salmon
    
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
        set pair_id, file(reads) from reads2
      
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
	file(cd_hits) from cd_hits_gmst
        file(sqlite_database) from sqlite_db1
	
    output:
	file("${cd_hits}*") into gmst_out
	file(genetable) into genetable
	file("gms.log") into gms_log

    when:
	params.gmst == true

    script:
        base =  cd_hits.baseName
    	genetable = "${base}.gene_tsv"
	gff = "${base}.gff"
	
"""

      gmst.pl \
      --fnn \
      --faa \
      --format GFF \
      ${cd_hits} \
      --verbose

      gff2genetable.py \
      -d ${sqlite_database} \
      -o ${genetable} \
       ${cd_hits}.gff
       
         
"""
//output sent to GhostKoala


}
