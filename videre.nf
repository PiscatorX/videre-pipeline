#!/usr/bin/env nextflow

/*
 * Copyright (c) 2019, Andrew Ndhlovu.
 *
 * author Andrew Ndhlovu <drewxdvst@outlook.com> 
 *  
 */

params.readtype		= "pe"
params.readsbase 	= "/home/andhlovu/Novogene/ftpdata.novogene.cn:2300/C101HW18111065/raw_data"
//params.readsbase       = "/home/drewx/Documents/subsample"
//params.readsbase        = "/home/andhlovu/subsample"
//params.sortmerna_db   = "${DB_REF}/SILVA/sel_SILVA.fasta"
params.sortmerna_db     = "${DB_REF}/SILVA/SILVA_132_SSURef_Nr99_tax_silva.fasta"
params.sortmerna_idx    = "${DB_REF}/SortMeRNA/SILVA.idx"
//params.sortmerna_db   = "${DB_REF}/SILVA_132_SSURef_Nr99_tax_silva.fasta"
//params.pe_patt 		= "*_{1,2}.fq.gz"
params.pe_patt 		= "*_{1,2}.fq.gz"
params.output  		= "$PWD/videre.Out"
sortmerna_db            = Channel.value(params.sortmerna_db)
sortmerna_idx           = Channel.value(params.sortmerna_idx)
output                  = params.output
DB_REF                  = System.getenv('DB_REF')
params.fastqc  		= false
params.trimmomatic      = true
params.sortmerna_index  = false
params.sortmerna        = true
params.megahit 		= false
params.metaspades 	= false
params.trinity          = false

reads                   = params.readsbase +'/'+ params.pe_patt
fileExt_glob            = "*" + reads.tokenize(".")[-1]
gz_ext                  = ( "."+reads.tokenize(".")[-1] == ".gz") ? ".gz" : ""


if (! sortmerna_db.endsWith('fasta')){

  error "\nSortMeRNA DB filename must end with `fasta` "    

}

Channel.fromFilePairs(reads)
       .ifEmpty{ error "Could not locate pair reads: ${reads}"}
       .into{reads1;
	    reads2;
	    reads3;
	    reads4;
	    readsx}


log.info """

=========================================================
  Videre: A nextlfow pipeline for metranscriptome data
=========================================================
 
Read type  		= ${params.readtype}
Read file pattern 	= ${reads}    
Output			= ${output}
Read ext. glob          = ${fileExt_glob}
High TP cores    	= ${params.htp_cores}
Midium TP cores    	= ${params.mtp_cores} 
Low TP cores    	= ${params.ltp_cores}
High memory  		= ${params.h_mem}
Medium memory           = ${params.m_mem}  	 
Low memory              = ${params.l_mem}  	   

FastQC  		= ${params.fastqc}           
Trimmomatic 		= ${params.trimmomatic}
MegaHit 		= ${params.megahit}          
Metaspades 	        = ${params.metaspades}       
SortmeRNA_idx           = ${params.sortmerna_index}         
SortmeRNA               = ${params.sortmerna}

Assemblers
Megahit			= ${params.megahit}
Metaspades 		= ${params.metaspades}
Trinity                 = ${params.trinity}       

Databases
SortmerRMA_DB           = ${params.sortmerna_db}
---------------------------------------------------------

Reads
=====
"""
readsx.each{  if(it instanceof List){println it} }

log.info"""
---------------------------------------------------------
"""


process fastqc_RawReads{

    //echo true
    cpus params.ltp_cores
    memory "${params.l_mem} GB"	
    input:
	set pair_id, file(reads) from reads1

    output:
       file("RawReadsQC/$pair_id") into fastqc_results

    when:
        params.fastqc


"""
   mkdir -p RawReadsQC/${pair_id} 
   fastqc \
   --extract \
   -f fastq \
   -o RawReadsQC/$pair_id\
   -t 2 ${fileExt_glob}
    
"""
}



process multiqc_RawReads{

    //echo true
    publishDir path: "$output/multiqc_RawReads", mode: 'move'
    cpus params.ltp_cores
    memory "${params.l_mem} GB"

    input:
	file("RawReadsQC/*") from fastqc_results.collect()

    output:
	set file("multiqc*"), file("RawReadsQC/*fastqc*")  into multiqc_report1
    
    when:
        params.fastqc


"""
    
    mv  RawReadsQC/*/*_fastqc*  RawReadsQC
    
    fastqc_combine.pl\
    -v\
    --out\
    RawReadsQC\
    --skip\
    --files 'RawReadsQC/*_fastqc'
     
    multiqc \
    RawReadsQC \
    -v 
    
"""

}



process trimmomatic{
    
    //echo true
    memory "${params.m_mem} GB"
    cpus  params.mtp_cores
    publishDir path: "$output/trimmomatic", mode: 'copy'
   
    
    input:
	set pair_id, file(reads) from reads2

    
    output:
        set pair_id, file("*_1P.fastq"), file("*_2P.fastq") into TrimmedReads
        set file("trim_${pair_id}.log"), file("${pair_id}.log")  into  trim_log
        file("${pair_id}_trim_{1,2}U.fastq") into unpairedReads
	

    when:
        params.trimmomatic

    script:	
    	(fwd, rev)=reads
    	
	    
"""

    $trimmomatic PE\
    $fwd\
    $rev\
    -baseout ${pair_id}_trim.fastq\
    -threads $params.htp_cores\
    -phred33\
    -trimlog  ${pair_id}.log\
    LEADING:10\
    TRAILING:10\
    SLIDINGWINDOW:25:10\
    MINLEN:50  2> trim_${pair_id}.log 
   
"""	
// If the name “mySampleFiltered.fq.gz” is provided, the following 4 file
// names will be used:
//     o mySampleFiltered_1P.fq.gz - for paired forward reads
//     o mySampleFiltered_1U.fq.gz - for unpaired forward reads
//     o mySampleFiltered_2P.fq.gz - for paired reverse reads
//     o mySampleFiltered_2U.fq.gz - for unpaired forward reads
}




if (! params.trimmomatic){

    reads3.map{ [it[0], it[1][0], it[1][1] ] }
        .set{no_trimm}
}


trimmed_reads = ( params.trimmomatic ? TrimmedReads : no_trimm )
    

process build_sortmerRNA_IDX{

    //echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
 
    input:
	val sortmerna_db
        val sortmerna_idx
    output:
        val sortmerna_idx_out
    
    when:
	params.sortmerna_index

    script:
	sortmerna_idx_out  = sortmerna_idx
    
             
"""
    
    indexdb_rna \
    --ref ${sortmerna_db},${sortmerna_idx} \
    -m ${params.h_mem} \
    -v 
     
"""

}


sortmerna_IDX = ( params.sortmerna_index ? sortmerna_idx_out : sortmerna_idx )


process sortmerRNA{

    //echo true
    cpus params.ltp_cores
    publishDir path: "${output}/sortmerna", mode: 'copy'
    memory "${params.m_mem} GB"

    input:
        val sortmerna_IDX
        val sortmerna_db
        set pair_id, file(forward_reads), file(reverse_reads) from trimmed_reads 
    
    output:
	set pair_id, file(forward_reads), file(reverse_reads) into filtered_reads 
        file(forward_reads) into mRNA_reads_fwd
        file(reverse_reads) into mRNA_reads_rev
	file("${pair_id}_lhist*") into lhist
        file("${pair_id}_sortmerna_aligned*") into SortMeRNA_Aligned      
	
       
    when:
	params.sortmerna

    script:
	paired=pair_id + ".fastq"
    
     	
"""

    reformat.sh \
    in=${forward_reads} \
    in2=${reverse_reads} \
    out=${paired}  

    rm -v ${forward_reads} ${reverse_reads}

    sortmerna \
    --ref  ${sortmerna_db},${sortmerna_IDX} \
    --reads ${paired} \
    --paired_out \
    --aligned ${pair_id}_sortmerna_aligned \
    --other mRNA \
    -a ${params.htp_cores} \
    --fastx \
    -m ${params.h_mem} \
    --log -v  

    rm -v ${paired}

    sed -i '/^\$/d' mRNA.fastq

    reformat.sh \
    deleteinput=t \
    in=mRNA.fastq \
    lhist=${pair_id}_lhist \
    out1=${forward_reads} \
    out2=${reverse_reads}

  
"""

}



process fastqc_filtered{
    
    //echo true
    cpus params.ltp_cores
    memory "${params.m_mem} GB"
    publishDir path: "$output/multiqc_filtered", mode: 'copy'
    input:
	set pair_id, file(fwd), file(rev) from filtered_reads

    output:
     	file("filtered/${pair_id}") into fastqc_results2
    when:
        params.fastqc

	
"""

   mkdir -pv  filtered/${pair_id}
   
   fastqc --extract\
   -f fastq \
   -o filtered/$pair_id\
   -t 2 *f*q
   
"""
	    
}


process multiqc_filtered{
    
    //echo true
    cpus params.ltp_cores
    memory "${params.m_mem} GB"
    publishDir path: "$output/multiqc_filtered", mode: 'move'
    
    input:
        file("filtered/*") from fastqc_results2.collect()
        file('*') from trim_log.collect()
	
    output:
    	file("multiqc*")  into multiqc_report2

    when:
        params.fastqc
	
"""

   mv -v  filtered/*/*_fastqc*  filtered/
   
   fastqc_combine.pl \
   -v\
   --out  filtered \
   --skip\
   --files  'filtered/*_fastqc'

   mv -v  trim_*.log  filtered/ 
   
   multiqc filtered\
   -v 
        
"""    

}


if (! params.sortmerna){

    reads4.into{reads_fwd; reads_rev}
    reads_fwd.map{ [it[1][0]] }
        .set{reads_fwd}
    reads_rev.map{ [it[1][1]] }
        .set{reads_rev}
    
}


(params.sortmerna ? mRNA_reads_fwd : reads_fwd).into{mRNA_fwd1; mRNA_fwd2; mRNA_fwd3}
(params.sortmerna ? mRNA_reads_rev : reads_rev).into{mRNA_rev1; mRNA_rev2; mRNA_rev3}


process megahit{
    
    //echo true
    cpus params.htp_cores 
    memory "${params.h_mem} GB"
    publishDir path: output, mode: 'move'
    
    
    input:
	file(all_fwd) from mRNA_fwd1.collect()
        file(all_rev) from mRNA_rev1.collect()
	
    when:
        params.megahit

    output:
        set file("MegaHit"), file('time_megahit') into MegahitOut
        
	
    script:
        fwd=all_fwd.join(",")
        rev=all_rev.join(",")
	   


""" 
  
    /usr/bin/time -v  -o time_megahit megahit \
    -1 $fwd \
    -2 $rev \
    -t ${params.htp_cores} \
    --tmp-dir /tmp \
    -o MegaHit\
    --out-prefix MegaHit \
    --verbose
   
    
"""
//Will fail if k=95 is not reached	   


}




process metaSpades{

    //echo true
    cpus params.htp_cores
    memory "${params.h_mem} GB"
    publishDir path: output, mode: 'move'

    input:
	file(all_fwd) from mRNA_fwd2.collect()
        file(all_rev) from mRNA_rev2.collect()
    
    when:
        params.metaspades
	

    output:
	set file("Metaspades"), file("time_metaspades") into MetaspadesOut
        
	
    script:
	 fwd=all_fwd.join(" ")
         rev=all_rev.join(" ")
    

"""  
    cat ${fwd} > fwd.fastq${gz_ext}
    cat ${rev} > rev.fastq${gz_ext}
    /usr/bin/time -v  -o time_metaspades  metaspades.py \
    -1 fwd.fastq${gz_ext} \
    -2 rev.fastq${gz_ext} \
    --only-assembler \
    -t ${params.htp_cores} \
    -m ${params.h_mem} \
    -o Metaspades
     
"""

}



process Trinity{

    //echo true
    cpus params.htp_cores
    memory "${params.h_mem} GB"
    publishDir path: output, mode: 'move'
    
    
    input:
	file(all_fwd) from mRNA_fwd3.collect()
        file(all_rev) from mRNA_rev3.collect()

    output:
        set file("Trinity"), file('time_Trinity') into Trinity
       
    when:
	params.trinity

    script:
        fwd=all_fwd.join(",")
    rev=all_rev.join(",")
    
	   
"""        
     /usr/bin/time -v  -o time_Trinity  Trinity\
     --seqType fq\
     --left  $fwd\
     --right $rev\
     --max_memory ${params.h_mem}G \
     --CPU ${params.htp_cores}\
     --output Trinity\
     --verbose

"""

}
