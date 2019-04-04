#!/usr/bin/env nextflow

params.readtype		= "pe"
params.readsbase 	= "/home/drewx/Documents/subsample"
//params.readsbase 	= "/home/andhlovu/Novogene/ftpdata.novogene.cn:2300/C101HW18111065/raw_data"
//params.readsbase 	= "/home/andhlovu/data"
params.se_patt 		= "*_RNA_1.fq.gz"
params.pe_patt 		= "*_RNA_{1,2}.fq" 
params.output  		= "$PWD/Videre.Out"
DB_REF                  = System.getenv('DB_REF')
params.readqc  		= false
params.megahit 		= true
params.metaspades 	= false
params.trinity          = false
params.quast 		= false
params.trimmomatic      = true
params.sortmerna_idx    = false
params.sortmerna        = false



//params.sortmerna_db     = "${DB_REF}/sel_SILVA.fasta"
params.sortmerna_db     = "${DB_REF}/SILVA_132_SSURef_Nr99_tax_silva.fasta"



if ( params.readtype.toLowerCase() == "se") {

    reads = params.readsbase +'/'+ params.se_patt  
    Channel.fromPath(reads)
	  .ifEmpty{ error "Could not locate pair reads: ${reads}"}
	  .map{it ->[it.baseName, [it]]}
	  .set{get_reads}

}else{

    reads = params.readsbase +'/'+ params.pe_patt
    Channel.fromFilePairs(reads)
	   .ifEmpty{ error "Could not locate pair reads: ${reads}"}
	   .set{get_reads}
}


if (! params.sortmerna_db ){
    
    error("SortMeRNA database not set")
}

fileExt_glob = "*" + reads.tokenize(".")[-1]

gz_ext = ( "."+reads.tokenize(".")[-1] == ".gz") ? ".gz" : ""





log.info """${gz_ext}"""



output = params.output
get_reads.into{reads1; reads2; reads3; readx; reads_cp}

sortmerna_db  = file(params.sortmerna_db)

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
Read QC                 = ${params.readqc}
H_mem  			= ${params.h_mem}
trimmomatic 		= ${params.trimm}

Assemblers
Megahit			= ${params.megahit}
Metaspades 		= ${params.metaspades}

Databases
SortmerRMA_DB           = ${params.sortmerna_db}
---------------------------------------------------------

Reads
=====
"""
readx.each{  if(it instanceof List){println it} }


log.info"""
---------------------------------------------------------
"""




if (params.readqc) {

process fastqc_RawReads{

    cpus 2 
    memory 1G
    publishDir path: output, mode: 'copy'
	
    input:
	set pair_id, file(reads) from reads1

    output:
       file("RawReadsQC/$pair_id") into fastqc_results



"""
   mkdir -p RawReadsQC/$pair_id 
   fastqc\
   --extract\
   -f fastq\
   -o RawReadsQC/$pair_id\
   -t 2 ${fileExt_glob}
    
"""

}




process multiqc_RawReads{

    //echo true
    publishDir path: "$output/multiqc_RawReads", mode: 'copy'
    cpus params.ltp_cores
    memory 2G

    input:
	file("RawReadsQC/**") from fastqc_results.collect()

    output:
	//use set for multiple file reduce channel
	set file("multiqc_report.html"), file("multiqc_data"), file("RawReadsQC/*fastqc*")  into MQC_report1
       

 """  
   mv  RawReadsQC/*/*_fastqc*  RawReadsQC

   fastqc_combine.pl\
   -v\
   --out\
   RawReadsQC\
   --skip\
   --files 'RawReadsQC/*_fastqc'
 
    multiqc\
    RawReadsQC\
    -v 
 """

}
}



if ( params.readtype.toLowerCase() == "se") {
   
process trimmomatic_SE {
    
    //echo true
    publishDir path: "$output/TrimmReads", mode: 'copy'
    cpus params.mtp_cores


    input:
	set pair_id, file(reads) from reads2

    output:
	set file("trim_${pair_id}.log"), file("trim.log")  into  trim_log
	file("*_Trimmed.fastq") into (fwd_reads1, fwd_reads2, fwd_reads3)
	set  val(pair_id), file("*_Trimmed.fastq") into (TrimmedReads1, TrimmedReads2, TrimmedReads3, TrimmedReads4)
	    
"""

    $trimmomatic SE\
    $reads\
    ${pair_id}_Trimmed.fastq\
    -threads $params.htp_cores\
    -phred33\
    -trimlog trim.log\
    LEADING:10\
    TRAILING:10\
    LEADING:10\
    TRAILING:10\
    SLIDINGWINDOW:25:10\
    MINLEN:50  2> trim_${pair_id}.log 

"""


}


}else{




process trimmomatic{
    
    //echo true
    memory params.m_mem
    cpus  params.mtp_cores
    publishDir path: "$output/Trimmomatic", mode: 'copy'
   
    
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
}



if (params.readqc) {
    
process fastqc_TrimmomaticReads{
    
    //echo true
    publishDir path: output, mode: 'copy'
    cpus 2
    memory 1G
    
    input:
	set pair_id, file(fwd), file(rev) from TrimmedReads1

    output:
     	file("TrimmomaticReadsQC/${pair_id}") into fastqc_results2

	
"""
   mkdir -pv TrimmomaticReadsQC/${pair_id}
   fastqc --extract\
   -f fastq\
   -o TrimmomaticReadsQC/$pair_id\
   -t 2 *fastq
   
"""

	    
}




process multiqc_TrimmomaticReads{
    
    //echo true
    cpus params.ltp_cores
    publishDir path: "$output/multiqc_TrimmomaticReads", mode: 'copy'
    
    input:
        file("TrimmomaticReadsQC/*") from fastqc_results2.collect()
        file('*') from trim_log.collect()
	
    output:
    	set file("multiqc_report.html"), file("multiqc_data"), file("TrimmomaticReadsQC/*fastqc*")  into MQC_report2      

	
"""
   mv -v  TrimmomaticReadsQC/*/*_fastqc*  TrimmomaticReadsQC/
   
   fastqc_combine.pl\
   -v\
   --out  TrimmomaticReadsQC\
   --skip\
   --files 'TrimmomaticReadsQC/*_fastqc'

   mv -v  trim_*.log  TrimmomaticReadsQC/ 
   
   multiqc TrimmomaticReadsQC\
   -v 
        
"""    

	    
}
}






if (params.trimm == false) {

    //workaround to test assembly
    //this skips trimmomatic trimming
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
    storeDir output
    
    
    input:
	file(all_fwd) from fwd_reads1.collect()
        file(all_rev) from rev_reads1.collect()
    
    when:
        params.megahit == true
	

    output:
        set file("MegaHit"), file('time_megahit') into MegahitOut
        file('MegaHit/MegaHit.fasta') into (megahit_contigs1, megahit_contigs2, megahit_contigs3, megahit_contigs4)
	
    script:
        fwd=all_fwd.join(",")
        rev=all_rev.join(",")
	   


""" 
  
    /usr/bin/time -v  -o time_megahit  megahit \
    -1 $fwd \
    -2 $rev \
    -t ${params.htp_cores} \
    --tmp-dir /tmp \
    -o MegaHit\
    --out-prefix MegaHit \
    --verbose
    mv MegaHit/MegaHit.contigs.fa  MegaHit/MegaHit.fasta
    #$TRINITY_HOME/util/TrinityStats.pl  MegaHit/MegaHit.fa  | tee MegaHit/contig_Nx_megahit.stats 
    #megahit_toolkit contig2fastg 95 MegaHit/intermediate_contigs/k95.contigs.fa >  MegaHit/k95.fastg     
    
"""
//Will fail if k=95 is not reached	   


}




process metaSpades{

    echo true
    cpus params.htp_cores
    memory "${params.h_mem} GB"
    publishDir path: output, mode: 'copy'
    //storeDir output


    input:
	file(all_fwd) from fwd_reads2.collect()
        file(all_rev) from rev_reads2.collect()
    
    when:
        params.metaspades == true
	

    output:
	set file("Metaspades"), file("time_metaspades") into MetaspadesOut
        file("Metaspades/metaspades.fasta") into (metaspades_contigs1, metaspades_contigs2)
	
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
    mv Metaspades/contigs.fasta  Metaspades/metaspades.fasta
    #$TRINITY_HOME/util/TrinityStats.pl  Metaspades/metaspades.fasta | tee Metaspades/contig_Nx_metaspades.stats 

"""

}




process Trinity{

    echo true
    cpus params.htp_cores
    //memory params.h_mem
    errorStrategy 'ignore'
    //publishDir path: output, mode: 'copy'
    storeDir output
    
    input:
       	file(all_fwd) from fwd_reads3.collect()
        file(all_rev) from rev_reads3.collect()
	

    output:
        set file("Trinity"), file('time_Trinity') into Trinity
        file("Trinity/Trinity.fasta") into trinity_contigs

    when:
	params.trinity ==  true

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
     $TRINITY_HOME/util/TrinityStats.pl  MegaHit/Trinity.fasta | tee Trinity/contig_Nx.stats 


"""


}


// // metaspades_contigs_1 = (params.metaspades) ?  metaspades_contigs1 : Channel.empty()
// // megahit_contigs_1 = (params.megahit) ? megahit_contigs1 :  Channel.empty()



process quast{

    //echo true
    cpus  params.htp_cores
    //errorStrategy 'ignore'
    publishDir path: output, mode: 'copy'
    
    input:
	file(metaspades_contig) from metaspades_contigs1
	file(megahit_contig) from megahit_contigs1
        file(trinity_contig) from trinity_contigs 	
    
    when:
        params.quast == true

    output:
	file("Quast") into QuastOut
        file("MetaQuast") into MetaQuastOut
    
    script:
        contig1  = metaspades_contigs1.val.getName()
        contig2  = megahit_contigs1.val.getName()
        contig3  = trinity_contigs.val.getName()
        all_contigs  = [contig1, contig2, contig3].join(" ")

    
""" 

   /usr/bin/time -v  -o time_quast quast\
   ${all_contigs} \
   -t ${params.htp_cores} \
   -o Quast  


   /usr/bin/time -v  -o time_metaquast metaquast.py \
   ${all_contigs} \
   -t ${params.htp_cores} \
   -o MetaQuast  


"""

}


  
process build_sortmerRNA_IDX{

    echo true
    cpus params.mtp_cores
    storeDir "${DB_REF}/SortMeRNA"
    memory params.m_mem
 
    input:
	file sortmerna_db
    
    output:
	megahit_contigs3
  	file "${sortmerna_idx}*" into SortMeRNA_idx
        
    when:
      params.sortmerna_idx  == true	
            
    script:
        sortmerna_idx = file(params.sortmerna_db).getName().replaceFirst(/fasta/, "idx")

    
"""
    
    indexdb \
    --ref ${sortmerna_db},${sortmerna_idx} \
    -m ${params.h_mem} \
    -v 
    
"""

}


process sortmerRNA{

    // echo true
    // cpus params.mtp_cores
    //publishDir path: "${output}/SortMeRNA", mode: 'copy'
    storeDir "${DB_REF}/SortMeRNA"
    
    memory params.m_mem
    input:
	val SortMeRNA_idx
        //file("MegaHit.contigs.fa") from  megahit_contigs_3
        file contig_fasta from megahit_contigs3
    
    output:
     	file("sortmerna_aligned.fasta") into SortMeRNA_Aligned
        file("sortmerna_contigs.fasta") into SortMeRNA_Contigs
    
    when:
      params.sortmerna == true
            
    script:
	sortmerna_idx = "${DB_REF}/SortMeRNA/" + file(params.sortmerna_db).getName().replaceFirst(/fasta/, "idx")
        
    
"""
    
    sortmerna \
    --ref  ${sortmerna_db},${sortmerna_idx} \
    --reads ${contig_fasta} \
    --aligned sortmerna_aligned \
    --other sortmerna_contigs \
    -a ${params.htp_cores} \
    --fastx \
    -m ${params.m_mem}
"""

}


