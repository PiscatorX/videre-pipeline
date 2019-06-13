#! /usr/bin/env nextflow

/*
 * Copyright (c) 2019, Andrew Ndhlovu.
 *
 * author Andrew Ndhlovu <drewxdvst@outlook.com> 
 *  
 */

params.readsbase   = "/home/andhlovu/MT-assembly-megahit/videre.Out/CD-Hit"
//params.readsbase    = "/home/drewx/Documents/subsample"
params.pe_patt      = "*_RNA_cdhit_{1,2}P.fastq"
params.DB_REF 	    =  "${DB_REF}"
params.output       = "${PWD}/Salmon"
params.cdHit_perc   = 0.98
output              =  params.output
//params.queries_path = "Contigs"
params.queries_path = "/home/andhlovu/MT-contigs"
query_seq           =  file(params.queries_path)
output              = params.output
DB_REF		    = params.DB_REF
params.bowtie_idx   = true
params.bowtie       = true
params.salmon_index = true
params.salmon_quant = true
params.gmst         = true


Channel.fromPath(params.queries_path +'/*')
    .ifEmpty{ error "Could not locate pair contigs files => ${params.queries_path}" }
    .set{contig_queries}


// Channel.fromPath(params.queries_path +'/*')
//     .ifEmpty{ error "Could not locate pair contigs files => ${params.queries_path}" }
//     .into{cd_hits_salmon; cd_hits_bowtie}
//Channel.value("MegaHit").into{contig_basename1contig_basename2,contig_basename3)

reads = params.readsbase +'/'+ params.pe_patt

Channel.fromFilePairs(reads)
       .ifEmpty{ error "Could not locate pair reads: ${reads}"}
       .into{reads1; reads2; reads3; reads4; readsx}




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

//readx.subscribe{  if(it instanceof List){ println it} }

readsx.each{  if(it instanceof List){println it} }

log.info"""
---------------------------------------------------------
"""


process cd_hit_est{
   
    //echo true
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${output}/CD-Hit", mode: 'copy'
  
   input:
       file(contigs) from contig_queries
       
    
   output:
       file(contig_fname) into (cd_hits_bowtie, cd_hits_salmon, cd_hits_gmst)
       file("${contig_basename}.cd_hits.clstr") into cdhit_clusters
       file("${contig_basename}.cd_hits") into cd_hits 
       file(sqlite_database) into sqlite_db
       val contig_basename into (contig_basename1,contig_basename2,contig_basename3)

    script: 
	contig_basename = "${contigs.baseName}_SHB"
        sqlite_database = "${contig_basename}.db"   
	contig_fname = "${contig_basename}.fasta"
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
    -o ${contig_basename}.cd_hits 

    contig_initDB.py \
    -d ${sqlite_database} \
    -c ${contig_basename}.cd_hits \
    -f \
    -o ${contig_fname}

"""

}




process bowtie_idx{

    //echo true
    cpus params.htp_cores
    memory "${params.m_mem} GB"
    publishDir "${DB_REF}/Bowtie"
    
    input:
        val  contig_basename1
        file contig_fasta from cd_hits_bowtie
 	
    output:  
        val bowtie2_base
        file("${bowtie2_base}*") into bowtie_idx
       
    when:
	params.bowtie_idx
    
    script:
        bowtie2_base =  "bowtie2_${contig_basename1}"
      
    
"""
    
    bowtie2-build  \
    --large-index \
    --threads ${params.htp_cores} \
    ${contig_fasta} \
    ${bowtie2_base}  > ${bowtie2_base}.log


    bowtie2-inspect \
    --large-in \
    --summary \
    ${bowtie2_base}  >  bowtie2_${contig_fasta}.idx_stats

    
"""    
    
}



if (!params.bowtie_idx){

   bowtie_idx = Channel.fromPath("${DB_REF}/Bowtie/*")
   
   contig_basename2.map{ fname ->
   			 def base = "bowtie2_${fname}"
   			 return base
   			 }.set{bowtie2_base}

}



 
process bowtie2bam{

    //echo true
    tag "${sample}"
    cpus params.htp_cores
    publishDir "${output}/Bowtie2sam", mode: "copy"
    memory "${params.m_mem} GB"

    input:
        each data from reads3
        file bowtie_idx_files from bowtie_idx.collect()
	val bowtie2_base

    output:
	file("${sample}*") into bowtie_sam
        file("${sample}.un")   into bowtie_unaligned
    
    when:
	params.bowtie

    script:
        (sample, reads) = data
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
//http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

}


process salmon_index{
    //echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    //storeDir "${params.DB_REF}/Salmon"
    publishDir "${DB_REF}", mode: "copy"

    input:
	file(cd_hits) from cd_hits_salmon
    
    output:
        file("salmon_index") into salmon_index

    when:
	params.salmon_index
        
	    
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


if(! params.salmon_index){

     salmon_index = Channel.fromPath("${DB_REF}/salmon_index")

}


process salmon_quant{
    
    //echo  true
    tag "${pair_id}"
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: "${output}/Quant", mode: 'move'
    
    input:
        each data from reads2
	file(index) from salmon_index
        
      
    output:
        file(pair_id) into salmon_quant

    when:
	params.salmon_quant
    
    script:
        (pair_id,reads) = data
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
    --output ${pair_id} \
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
	file sqlite_db

    output:
	file(genetable) into genetable
	file("gms.log") into gms_log
	file("*.faa") into predicted_aa
	file("*.fnn") into predicted_nn
	file("*.gff") into pedicted_gff
	file(sqlite_db) into sqlite_database

    when:
	params.gmst

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
      
      fasta_gff_dedup.py \
      ${cd_hits}.faa \
      ${cd_hits}.gff \
      -f ${cd_hits}_aa.faa \
      -g ${cd_hits}_aa.gff

      fasta_gff_dedup.py \
      ${cd_hits}.fnn \
      ${cd_hits}.gff \
      -f ${cd_hits}_nt.faa \
      -g ${cd_hits}_nt.gff

      gff2genetable.py \
      -d ${sqlite_db} \
      -o ${genetable} \
       ${cd_hits}.gff
     
      
"""

//output sent to GhostKoala


}


