#!/usr/nin/bash

#this script recursively downloads the imicrobe MMETSP data using sample and data identifiers
#http://datacommons.cyverse.org/browse/iplant/home/shared/imicrobe/projects/104/sample-attr.tab

$1=data_type

while read sample_id,mmetsp_ref
do

  echo sample_id,mmetsp_ref

done <  sample_id.map

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
