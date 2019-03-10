


process diamond{
    
    echo true
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


"""
    
}



    // diamond blastx -d nr \
    // --un diamond.unaligned \
    // --al diamond.aligned \
    // -q  Cd_Hit_0.98.cd_hits \
    // -o matches.dmnd \
    // --more-sensitive \
    // --evalue 1e-5 \
    // --top 90 \
    // --id 40 \
    // -v    
