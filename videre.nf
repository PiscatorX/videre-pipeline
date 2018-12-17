#!/usr/bin/env nextflow

params.reads = "$baseDir/data1/*_{1,2}.f*q"
output = "$PWD/Videre.Out"

       
Channel.fromFilePairs(params.reads)
       .ifEmpty{ error "Could not locate reads: ${params.reads}"}
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
    publishDir path: "$output/trimmoQC", mode: 'copy'
   
    input:
	set pair_id, file(reads) from read_pair2

    output:
	file("${pair_id}_trim*") into trimmedreads
   
    script:
	(fwd,rev)=reads

	
"""
     
    trimmomatic PE  $fwd $rev -baseout ${pair_id}_trimmer.fastq -phred33 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50 2>  trim_out.log
   
     
"""
// If the name “mySampleFiltered.fq.gz” is provided, the following 4 file
// names will be used:
//     o mySampleFiltered_1P.fq.gz - for paired forward reads
//     o mySampleFiltered_1U.fq.gz - for unpaired forward reads
//     o mySampleFiltered_2P.fq.gz - for paired reverse reads
//     o mySampleFiltered_2U.fq.gz - for unpaired

	
}



// process fastqc_RawReads{
//     echo true
//     publishDir path: output, mode: 'copy'
	
//     input:
// 	set pair_id, file(reads) from read_pair1

//     output:
// 	file("raw_reads/$pair_id") into fastqc_results
    
// """
   
//    mkdir -pv raw_reads/$pair_id
//    fastqc --extract -f fastq  -o raw_reads/$pair_id  -t 1  *q
//    fastqc_combine.pl -v -o raw_reads/$pair_id --skip --files "raw_reads/$pair_id/*_fastqc"
//    multiqc raw_reads/$pair_id
// """
// }



// process {
 
//     echo true
//     publishDir path: "$output/trimmoQC", mode: 'copy'
	
//     input:
// 	file() 
	

// }
