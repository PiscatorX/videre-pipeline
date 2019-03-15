























process cd_hit_est{
    
   //echo true
    cpus  params.htp_cores
    publishDir path: "${output}/CD-Hit", mode: 'copy'
  
   input:
       //file(contigs) from metaspades_contigs2
       file(contigs) from megahit_contigs2
       
   output:
       file "Cd_Hit_${params.cdHit_perc}.cd_hits" into (cd_hits1, cd_hits2, cd_hits3)
       file "Cd_Hit_${params.cdHit_perc}.cd_hits.clstr" into cdhit_clusters


"""
    cd-hit-est \
    -i $contigs \
    -c ${params.cdHit_perc} \
    -T ${params.htp_cores} \
    -d 0 \
    -r 0 \
    -p 1 \
    -g 1 \
    -o Cd_Hit_${params.cdHit_perc}.cd_hits 
     ls -l
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
        set pair_id, file(reads) from reads3
      
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



