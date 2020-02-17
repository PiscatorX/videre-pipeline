# Videre 

Videre is a [nexflow](https://www.nextflow.io/) computational pipeline for metatranscriptome data assembly and analysis

## **Getting Started**

The videre pipeline take raw RNA fastq reads and assembles them into a metatranscrptome. To get started clone/download this repository to your local machine and install the prerequisites in the list below. 

## Prerequisites

### Nextflow

Nextflow requires Java to to run, for more details are avialable from [https://www.nextflow.io/](https://www.nextflow.io/)

###  Miniconda

Due to the multiple number of programs required for this pipeline [miniconda](https://docs.conda.io/en/latest/miniconda.html)  is highly recommended for package and environment management. 

###  Installing

All tools and dependencies used in this pipeline are listed in the  [videre.yml](https://github.com/PiscatorX/videre-pipeline/blob/master/videre.yml) file. The environment may be created using ``conda env create -f videre.yml``

### Running the pipeline

There are main pipeline scripts are 
* videre.nf is the main nextflow pipeline script for sequence read assembly to generate the metascriptome.   
* videre-salmon.nf is the nextflow script predicting the protein coding regions from the metatranscriptome.
* videre-diamond.nf is for annotating the protein coding regions and assigning taxonomy using the [MMETSP](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001889)
* videre-hmmer.nf is for assigning Pfam (protein domain families) using  hmmscan


## Authors

* **Andrew Ndhlovu** 