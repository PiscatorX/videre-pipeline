#!/usr/bin/env bash

cut -d " "  -f 1,3  $1   | grep  TAXON_ID | tr  -d ">"  |  sed -e 's|/TAXON_ID=||'  >  $2
wc  -l  $2
cut -d " "  -f 1,4  $1   | grep  TAXON_ID | tr  -d ">"  |  sed -e 's|/TAXON_ID=||'  >> $2
wc -l   $2
