#!/usr/bin/env nextflow

   //params.reads = "$baseDir/DataX/*_{1,2}.f*q"
params.reads = "$baseDir/data/reads.{left,right}.fq"
output = "$PWD/Videre.Out"

       
Channel.fromFilePairs(params.reads)
       .ifEmpty{ error "Could not locate or pair reads: ${params.reads}"}
       .set{paired_reads}
paired_reads.into{read_pair1; read_pair2}
	     







process fastqc_RawReads{
    //echo true
    publishDir path: output, mode: 'copy'
	
    input:
	set pair_id, file(reads) from read_pair1

    output:
	file("RawReadsQC/$pair_id") into fastqc_results
    
"""
   
   mkdir -pv RawReadsQC/$pair_id
   fastqc --extract -f fastq  -o RawReadsQC/$pair_id -t 1  *f*q
   
"""

}




process multiqc_RawReads{
    //echo true
    publishDir path: "$output/RawReadsQC/multiqc", mode: 'copy'
    input:
        file("RawReadsQC/*") from fastqc_results.collect()

    output:
    	file("multiqc_report.html") into MQC_report1
    	file("multiqc_data") into MQC_data1
    	file("RawReadsQC/*fastqc*")   into FSQ_results1   

	
"""

   mv  RawReadsQC/*/*_fastqc*  RawReadsQC
   fastqc_combine.pl -v --out  RawReadsQC --skip --files 'RawReadsQC/*_fastqc'
   multiqc RawReadsQC
        
"""    

}




process trimmomatic{

    //echo true
    publishDir path: "$output/TrimmoReads", mode: 'copy'
   
    input:
	set pair_id, file(reads) from read_pair2

    output:
    file("${pair_id}_trim_{1,2}U.fastq") into unpairedReads
    set  val(pair_id), file('*_1P.fastq'), file('*_2P.fastq') into TrimmedReads1,
	                                                           TrimmedReads2,
	                                                           TrimmedReads3
	
   
    script:
	(fwd,rev)=reads

	
"""
     
    trimmomatic PE  $fwd $rev -baseout ${pair_id}_trim.fastq -phred33 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50 2>  trim_out.log
   
     
"""
// If the name “mySampleFiltered.fq.gz” is provided, the following 4 file
// names will be used:
//     o mySampleFiltered_1P.fq.gz - for paired forward reads
//     o mySampleFiltered_1U.fq.gz - for unpaired forward reads
//     o mySampleFiltered_2P.fq.gz - for paired reverse reads
//     o mySampleFiltered_2U.fq.gz - for unpaired

	
}





process fastqc_TrimmReads{
    //echo true
    publishDir path: output, mode: 'copy'
	
    input:
	set pair_id, file(fwd), file(rev) from TrimmedReads1

    output:
     	file("TrimmoReadsQC/$pair_id") into fastqc_results2

	
"""
   
   mkdir -pv TrimmoReadsQC/$pair_id
   fastqc --extract -f fastq  -o TrimmoReadsQC/$pair_id -t 1  *fastq
   
"""

}


process multiqc_TrimmReads{
    //echo true
    publishDir path: "$output/TrimmoReadsQC/multiqc", mode: 'copy'
    input:
        file("TrimmoReadsQC/*") from fastqc_results2.collect()

    output:
    	file("multiqc_report.html") into MQC_report2
    	file("multiqc_data") into MQC_data2
    	file("TrimmoReadsQC/*fastqc*")   into FSQ_results2   

	
"""

   mv  TrimmoReadsQC/*/*_fastqc*  TrimmoReadsQC
   fastqc_combine.pl -v --out  TrimmoReadsQC --skip --files 'TrimmoReadsQC/*_fastqc'
   multiqc TrimmoReadsQC
        
"""    

}




process megahit{
    
    echo true
    publishDir path: "$output/MegaHit", mode: 'copy'
    input:
	set  pair_id, file(fwd), file(rev) from TrimmedReads2

    output:
       file("$pair_id") into megahitResults

	
"""

    megahit -1 $fwd  -2 $rev  -m 0.75  -t 1  -o $pair_id  
    

"""

}




process trinity{

    echo true
    publishDir path: "$output/Trinity", mode: 'copy'
	
    input:
        set pair_id, file(fwd), file(rev) from TrimmedReads3
    output:
         file("trinity_${pair_id}") into trinityOut


"""

   Trinity --seqType fq --max_memory 4G --no_normalize_reads \
          --left $fwd  --right $rev --output trinity_${pair_id} --CPU 4
   $TRINITY_HOME/util/TrinityStats.pl  trinity_${pair_id}/Trinity.fasta | tee trinity_${pair_id}/Trinity.stats 
   $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl  trinity_${pair_id}/transcripts.TMM.EXPR.matrix \
trinity_${pair_id}/Trinity.fasta | tee  trinity_${pair_id}/Trinity.ExN50.stats
 
"""

}


