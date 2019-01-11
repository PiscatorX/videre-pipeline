#!/usr/bin/env nextflow

//params.reads = "$baseDir/DataX/*_{1,2}.f*q"
params.reads = "$baseDir/data/reads?.{left,right}.fq"
params.nr_faa = Channel.fromPath("/home/drewx/Documents/videre-pipeline/data/nr.faa")	
params.cd_hit_threads=4
params.threads = 4
output = "$PWD/Videre.Out"
cd_hit_clusters = Channel.from(0.80, 0.98)
       
Channel.fromFilePairs(params.reads)
       .ifEmpty{ error "Could not locate or pair reads: ${params.reads}"}
       .set{paired_reads}
       paired_reads.into{read_pair1; read_pair2; read_pair3}


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



process trimmomatic {

    //echo true
    publishDir path: "$output/TrimmoReads", mode: 'copy'
   
    input:
	set pair_id, file(reads) from read_pair2

    output:
        file("*_1P.fastq") into (fwd_reads1, fwd_reads2)          
	file('*_2P.fastq') into (rev_reads1, rev_reads2)
        file("${pair_id}_trim_{1,2}U.fastq") into unpairedReads
	set  val(pair_id), file("*_1P.fastq"), file("*_2P.fastq") into (TrimmedReads1, TrimmedReads2, TrimmedReads3)
	    
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

    //echo true
    publishDir path: output, mode: 'copy'
    input:
        val all_fwd from fwd_reads1.flatMap{ it.getName() }.collect()
	val all_rev from rev_reads1.flatMap{ it.getName() }.collect()
	file(reads) from TrimmedReads3.collect()
	

    output:
        file("MegaHit") into MegahitOut
	file('MegaHit/MegaHit.contigs.fa') into Contigs_fasta
	    
    script:
	f =  all_fwd.join(',')
	r =  all_rev.join(',')

	    
	    
"""       
     megahit -1 $f  -2  $r -m  0.75  -t 4  -o MegaHit  --out-prefix MegaHit --verbose
     $TRINITY_HOME/util/TrinityStats.pl  MegaHit/final.contigs.fa | tee MegaHit/contig_Nx.stats 
     megahit_toolkit contig2fastg 95 MegaHit/intermediate_contigs/k95.contigs.fa >  MegaHit/k95.fastg     

"""
	    
}




// process trinity{

//     echo true
//     publishDir path: "$output/Trinity", mode: 'copy'
	
//     input:
//         set pair_id, file(fwd), file(rev) from TrimmedReads3
//     output:
//          file("trinity_${pair_id}") into trinityOut


// """

//    Trinity --seqType fq --max_memory 4G --no_normalize_reads \
//           --left $fwd  --right $rev --output trinity_${pair_id} --CPU 4
//    $TRINITY_HOME/util/TrinityStats.pl  trinity_${pair_id}/Trinity.fasta | tee trinity_${pair_id}/Trinity.stats 
//    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl  trinity_${pair_id}/transcripts.TMM.EXPR.matrix \
// trinity_${pair_id}/Trinity.fasta | tee trinity_${pair_id}/Trinity.ExN50.stats
 
// """

// }




process cd_hit{
    
    //echo true
   maxForks 1
   publishDir path: "$output/CDHIT", mode: 'copy'
  
   input:
       each cluster_perc from  cd_hit_clusters
       file(contigs) from Contigs_fasta
   
   output:
       file "Cd_Hit_${cluster_perc}.cd_hits" into (cd_hits1, cd_hits2, cd_hits3)
       file "Cd_Hit_${cluster_perc}.cd_hits.clstr" into cdhit_clusters
       
       
"""
    cd-hit-est -i $contigs -c $cluster_perc -T $params.cd_hit_threads -d 0 -r 0 -p 1 -g 1  -o \
    Cd_Hit_${cluster_perc}.cd_hits 
    
"""

}



process salmon_index{
    
    //echo  true
    publishDir path: output, mode: 'copy'
    input:
	file(cd_hits) from cd_hits1.collect()
       
    output:
        file("Salmon") into salmon_index
    
"""

    salmon index -t  Cd_Hit_0.98.cd_hits -i  Salmon/transcripts_index --type quasi -p 4 -k 31     
    

"""
	
}




process salmon_quant{
    
    echo  true
    publishDir path: output, mode: 'copy'
    input:
	file(index) from salmon_index
        set pair_id, file(reads) from read_pair3
      
    output:
        file("Salmon") into Salmon_quant

    
    script:
    	(left, right)=reads

    
"""

    salmon quant -l IU --index Salmon/transcripts_index -1 $left  -2 $right --validateMappings --meta  --output Salmon -p 4

"""
	
}




    
process gmst{
     
    echo  true
    publishDir path: "$output/Gmst", mode: 'copy'
    input:
	file(cd_hits) from cd_hits2.collect()

    output:
	file("Cd_Hit_0.98.cd_hits.*") into gmst_out
	file("gms.log") into gms_log
		    
"""
      gmst.pl --fnn -faa   Cd_Hit_0.98.cd_hits  --verbose
         

"""
//output sent to GhostKoala

}




process diamond{    
    echo true
    maxForks 1
    publishDir path: "$output/Diamond", mode: 'copy'
    input:
       file(nr_faa) from params.nr_faa
       file(reads_fna) from cd_hits3.collect()

    output:
       file('matches.dmnd') into diamond_matches
       file('nr.dmnd')      into diamond_ref
       file('diamond.unaligned') into diamond_un
       file('diamond.aligned') into diamond_al
    

"""

    diamond makedb --in $nr_faa -d nr --threads $params.threads -v
    diamond blastx -d nr \
--un diamond.unaligned \
--al diamond.aligned \
-q  Cd_Hit_0.98.cd_hits \
-o matches.dmnd \
--more-sensitive \
--evalue 1e-5 \
--top 90 \
--id 40 \
-v    

"""
    
}


