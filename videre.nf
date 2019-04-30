#!/usr/bin/env nextflow

params.readtype		= "pe"
//params.readsbase 	= "/home/andhlovu/Novogene/ftpdata.novogene.cn:2300/C101HW18111065/raw_data"
//params.readsbase 	= "/home/drewx/Documents/subsample-tiny"
//params.readsbase      = "/home/andhlovu/subsample"
params.readsbase	= "/home/drewx/Documents/sea-biome/reads.gz"
//params.sortmerna_db   = "${DB_REF}/SILVA/sel_SILVA.fasta"
params.sortmerna_db     = "${DB_REF}/SILVA/SILVA_132_SSURef_Nr99_tax_silva.fasta"
params.sortmerna_idx    = "${DB_REF}/SortMeRNA/SILVA.idx"
//params.sortmerna_db   = "${DB_REF}/SILVA_132_SSURef_Nr99_tax_silva.fasta"
params.pe_patt 		= "*_{1,2}.fq.gz"
//params.pe_patt 		= "*_RNA_{1,2}.fq.gz" 
params.output  		= "$PWD/videre.Out"
sortmerna_db            = Channel.value(params.sortmerna_db)
sortmerna_idx           = Channel.value(params.sortmerna_idx)
output                  = params.output
DB_REF                  = System.getenv('DB_REF')
params.fastqc  		= true
params.trimmomatic      = true
params.sortmerna_index  = false
params.sortmerna        = false
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
	    reads5;
	    reads6;
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
SortmeRNA_idx           = ${params.sortmerna_idx}         

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
        params.fastqc == true


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
        params.fastqc == true


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
    
    echo true
    memory "${params.m_mem} GB"
    cpus  params.mtp_cores
    publishDir path: "$output/trimmomatic", mode: 'copy'
   
    
    input:
	set pair_id, file(reads) from reads2

    
    output:
        set val(pair_id), file("*_1P.fastq"), file("*_2P.fastq") into (TrimmedReads1, TrimmedReads2, TrimmedReads3, TrimmedReads4)
	file("*_1P.fastq") into (fwd_reads1, fwd_reads2, fwd_reads3, fwd_reads4 )
        file('*_2P.fastq') into (rev_reads1, rev_reads2, rev_reads3, rev_reads4)	
        set file("trim_${pair_id}.log"), file("${pair_id}.log")  into  trim_log
        file("${pair_id}_trim_{1,2}U.fastq") into unpairedReads
	file('time_trimmomatic') into time

    when:
        params.trimmomatic == true

    script:	
    	(fwd, rev)=reads
    	
	    
"""

    /usr/bin/time -v -o time_trimmomatic $trimmomatic PE\
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


process build_sortmerRNA_IDX{

    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
 
    input:
	val sortmerna_db
        val sortmerna_idx
    output:
        val sortmerna_idx_out
    
    when:
	params.sortmerna_idx == true

    script:
	sortmerna_idx_out  = sortmerna_idx
    
             
"""
    
    indexdb_rna \
    --ref ${sortmerna_db},${sortmerna_idx} \
    -m ${params.h_mem} \
    -v 
     
"""

}


if (! params.sortmerna_idx) {

    sortmerna_idx  = sortmerna_db.toString().replaceFirst(/fasta/, "idx") 
}


if (! params.trimmomatic){

    reads5.map{ [it[0], it[1][0], it[1][1] ] }
        .set{TrimmedReads1}
    
 }



process sortmerRNA{

    //echo true
    cpus params.htp_cores
    //publishDir path: "${output}/SortMeRNA", mode: 'copy'
    publishDir path:"${DB_REF}/SortMeRNA"
    memory "${params.m_mem} GB"

    input:
        val sortmerna_idx
        val sortmerna_db
        set pair_id, file(forward_reads), file(reverse_reads) from TrimmedReads1 
    
    output:
	set pair_id, file(forward_reads), file(reverse_reads) into filtered_reads 
        file(forward_reads) into (mRNA_fwd1, mRNA_fwd2, mRNA_fwd3)
	file(reverse_reads) into (mRNA_rev1, mRNA_rev2, mRNA_rev3)
        //(mRNA_reads1, mRNA_reads2, mRNA_reads3)
        //file("sortmerna_aligned.fastq*") into SortMeRNA_Aligned      
	// set pair_id, file(forward_reads), file(reverse_reads) into (mRNA_tags1, mRNA_tags2, mRNA_tags3)
        // file("mRNA.fastq") optional true into  (mRNA_reads1, mRNA_reads2, mRNA_reads3)
     	// file("sortmerna_aligned.fastq*") into SortMeRNA_Aligned       
    when:
	params.sortmerna == true

     

"""

  merge-paired-reads.sh \
  ${forward_reads} \
  ${reverse_reads} \
  ${pair_id}.fastq

  sortmerna \
  --ref  ${sortmerna_db},${sortmerna_idx} \
  --reads ${pair_id}.fastq\
  
  --paired_out \
  --aligned sortmerna_aligned \
  --other mRNA \
  -a ${params.htp_cores} \
  --fastx \
  -m ${params.h_mem} \
  --log \
  -v   
  
  rm -v \
  ${forward_reads} \
  ${reverse_reads} \
  ${pair_id}.fastq
 
  reformat.sh \
  in=mRNA.fastq \
  out1=${forward_reads} \
  out2=${reverse_reads}

  #unmerge-paired-reads.sh \
  mRNA.fastq \
  ${forward_reads} \
  ${reverse_reads}
      
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
        params.fastqc == true

	
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
        params.fastqc == true
	
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


if (params.trimmomatic == false) {

    //workaround to test assembly
    //this skips trimmomatic trimming
    set pair_id, file(forward_reads), file(reverse_reads) into (mRNA_reads1, mRNA_reads2, mRNA_reads3)
    reads_cp.into{cp_reads1; cp_reads2}        
    cp_reads1.map{ it[1][0] }
        .into{fwd_reads1; fwd_reads2; fwd_reads3; fwd_reads4}
    cp_reads2.map{ it[1][1] }
        .into{rev_reads1; rev_reads2; rev_reads3; rev_reads4}          
}



process megahit{
    
    echo true
    cpus params.htp_cores 
    memory "${params.h_mem} GB"
    publishDir path: output, mode: 'move'
    
    
    input:
	file(all_fwd) from mRNA_fwd1.collect()
        file(all_rev) from mRNA_rev1.collect()
	
    when:
        params.megahit == true

    output:
        set file("MegaHit"), file('time_megahit') into MegahitOut
        file('MegaHit/MegaHit.fasta') into (megahit_contigs1, megahit_contigs2, megahit_contigs3, megahit_contigs4)
	
    script:
        fwd=all_fwd.join(",")
        rev=all_rev.join(",")
	   


""" 
  
    megahit \
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




// process metaSpades{

//     echo true
//     cpus params.htp_cores
//     memory "${params.h_mem} GB"
//     publishDir path: output, mode: 'move'
//     //storeDir output


//     input:
// 	file(all_fwd) from mRNA_fwd2.collect()
//         file(all_rev) from mRNA_rev2.collect()
    
//     when:
//         params.metaspades == true
	

//     output:
// 	set file("Metaspades"), file("time_metaspades") into MetaspadesOut
        
	
//     script:
// 	 fwd=all_fwd.join(" ")
//          rev=all_rev.join(" ")
    

// """  
//     cat ${fwd} > fwd.fastq${gz_ext}
//     cat ${rev} > rev.fastq${gz_ext}
//     /usr/bin/time -v  -o time_metaspades  metaspades.py \
//     -1 fwd.fastq${gz_ext} \
//     -2 rev.fastq${gz_ext} \
//     --only-assembler \
//     -t ${params.htp_cores} \
//     -m ${params.h_mem} \
//     -o Metaspades
     
// """

// }



// process Trinity{

//     echo true
//     cpus params.htp_cores
//     memory "${params.h_mem} GB"
//     publishDir path: output, mode: 'move'
//     //storeDir output
    
//     input:
// 	file(all_fwd) from mRNA_fwd3.collect()
//         file(all_rev) from mRNA_rev3.collect()

//     output:
//         set file("Trinity"), file('time_Trinity') into Trinity
//         file("Trinity/Trinity.fasta") into trinity_contigs

//     when:
// 	params.trinity ==  true

//     script:
//         fwd=all_fwd.join(",")
//         rev=all_rev.join(",")
	   
// """        
//      /usr/bin/time -v  -o time_Trinity  Trinity\
//      --seqType fq\
//      --left  $fwd\
//      --right $rev\
//      --max_memory ${params.h_mem}G \
//      --CPU ${params.htp_cores}\
//      --output Trinity\
//      --verbose
    

// """

// }

