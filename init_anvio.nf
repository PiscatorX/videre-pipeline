#! /usr/bin/env nextflow

params.sample_name      = "St_Helana_Bay"
params.DB_REF 	        = System.getenv('DB_REF')
DB_REF		        = params.DB_REF
params.output   	= "${PWD}/Anvio"
params.bamdir		= "${DB_REF}/Bowtie2sam"
params.bam_patt 	= "*.bam"
params.contig  		= "/home/drewx/Documents/videre-pipeline/Salmon/CD-Hit/MegaHitX.fasta"
params.gene_table 	= "/home/drewx/Documents/videre-pipeline/Salmon/GeneMarkST/MegaHitX.gene_tsv"
params.min_contig       = 100
output 			= params.output 
contig 			= file(params.contig)
gene_table 		= file(params.gene_table)

reads = params.bamdir +'/'+ params.bam_patt  
Channel.fromPath(reads)
      .ifEmpty{ error "Could not locate bam files: ${reads}"}
      .into{bam_files1; bam_files2}



process gen_contigDB{
    label 'anvio'
    echo true
    cpus params.ltp_cores
    memory "${params.m_mem} GB"
    //errorStrategy 'ignore'
    publishDir path: output, mode: 'copyNoFollow'
    
    input:
        file contig
        file gene_table
	
    output:
	file "${contig.baseName}_anvio.db"  into (contigDB1, contigDB2)  

    
"""

  
    anvi-gen-contigs-database \
    -f ${contig} \
    -o ${contig.baseName}_anvio.db \
    --external-gene-calls ${gene_table} \
    -n 'Anvio database: ${contig}'

    #anvi-run-hmms \
    #-c ${contig.baseName}_anvio.db

"""

 
}



process init_bam{

    label  'anvio'
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: output, mode: 'copyNoFollow'
    
    
    input:
	file bam from bam_files1
        file contig
    
    output:
	set file("${bam_dir}/*.bam"), file("${bam_dir}/*.bai")  into anvio_bam

   script:
      bam_dir =  bam.baseName

"""
    mkdir -v ${bam_dir}
    anvi-init-bam \
    ${bam} \
    -o ${bam_dir}/${bam} 

"""
    
}



process anvi_profile{

    label 'anvio'
    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: output, mode: 'copy'
    
    input:
	set file(bam), file(index) from  anvio_bam
	file contig_db from  contigDB1 

    output:
	file("${sample}") into anvi_profiles
	//file("${sample}/PROFILE.db") into anvi_profiles
        set file("${sample}/AUXILIARY-DATA.db"),  file("${sample}/RUNLOG.txt") into profile_data
        
    
    script:
	bam_filename = bam.getName()
        sample = bam.baseName
	
	
"""
 
  anvi-profile \
  -i ${bam} \
  -c ${contig_db} \
  --min-contig-length ${params.min_contig} \
  --output-dir $sample \
  --sample-name $sample
    
"""

}

  


process anvi_merge{

    label 'anvio'
    echo true
    cpus params.mtp_cores
    memory "${params.m_mem} GB"
    publishDir path: output, mode: 'move'

    input:
	file(profile) from anvi_profiles.collect()
        file contig_db from contigDB2 

    script:
	

    
"""
    anvi-merge */PROFILE.db  \
    -o SAMPLES-MERGED \
    -c ${contig_db}    
    --sample-name

"""

}

