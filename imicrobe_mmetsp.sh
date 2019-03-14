#!/usr/bin/env bash

#this script recursively downloads the imicrobe MMETSP data using sample and data identifiers
#http://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/104/sample-attr.tab
#cut -f 1,2  sample-attr.tab >  sample_id.map
#https://de.cyverse.org/anon-files//iplant/home/shared/imicrobe/projects/104/samples/1662/MMETSP0982.pep.fa



data_type=$1
data_file=$2

if [ -z data_type ]

then
    echo ".cds.fa"
    echo ".cds.fa.uproc.kegg"
    echo ".cds.fa.uproc.kegg.annotated.tab"
    echo ".cds.fa.uproc.pfam28"
    echo ".cds.fa.uproc.pfam28.annotated.tab"
    echo ".cds.json"
    echo ".fastq.tar"
    echo ".nt.fa"
    echo ".pep.fa"
    echo "_clean.fasta"
exit 1

fi

if [ -z data_file ]

then
    echo "No datafile provided!"
    echo "provide .tsv datafile with 2 columns: Sample_id mmetsp_id"
fi



while read sample_id mmetsp_id
do

  wget https://de.cyverse.org/anon-files//iplant/home/shared/imicrobe/projects/104/samples/${sample_id}/${mmetsp_id}.pep.fa
  
done <  $2

# MMETSP0982-Undescribed-sp-CCMP2436.1.fastq
# MMETSP0982-Undescribed-sp-CCMP2436.2.fastq
# MMETSP0982.cds.fa
# MMETSP0982.cds.fa.uproc.kegg
# MMETSP0982.cds.fa.uproc.kegg.annotated.tab
# MMETSP0982.cds.fa.uproc.pfam28
# MMETSP0982.cds.fa.uproc.pfam28.annotated.tab
# MMETSP0982.cds.json
# MMETSP0982.fastq.tar
# MMETSP0982.nt.fa
# MMETSP0982.pep.fa
# MMETSP0982_clean.fasta

