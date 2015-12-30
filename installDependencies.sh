#!/usr/bin/env bash
#Ubuntu Linux Dependencies
sudo apt-get update && apt-get upgrade -y
sudo apt-get install -y gcc-4.6-base && apt-get install -y cpp-4.6 && apt-get install -y libgomp1 && apt-get install -y libquadmath0 && apt-get install -y libc6-dev && apt-get install -y build-essential && apt-get install -y zlib1g-dev && apt-get install -y libc6-dev
sudo apt-get install -y libcurl4-openssl-dev && apt-get install -y libxml2-dev && apt-get install -y curl && apt-get install -y cmake
sudo apt-get install -y zlibc && apt-get install -y zlib1g-dev && apt-get install -y libhdf5-dev
sudo apt-get update
sudo apt-get update && apt-get install -y git
sudo mkdir /KallistoSource && cd /KallistoSource && git clone https://github.com/pachterlab/kallisto.git && cd ./kallisto && mkdir ./build && cd ./build && cmake .. && make && make install
sudo echo "deb http://cran.stat.ucla.edu/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update && apt-get install -y r-base && apt-get install -y r-base-dev && apt-get update 
#Artemis Dependency
sudo R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BiocInstaller")'
sudo R -e 'library(BiocInstaller);biocLite(c("rtracklayer","XML","biomaRt","RCurl"),ask=FALSE)'
sudo R -e 'library(BiocInstaller);biocLite(c("parallel","jsonlite","qusage","GenomeInfoDb","limma","Biobase","SummarizedExperiment","SPIA","EnrichmentBrowser","clusterProfiler","rhdf5","matrixStats","GenomicRanges","GenomicFeatures","Matrix","KEGGREST","beeswarm","tools","graphite","roxygen2","knitr"),ask=FALSE )'
sudo R -e 'library(BiocInstaller);biocLite(c("Homo.sapiens","Mus.musculus","RUVSeq","erccdashboard"),ask=FALSE)'
#TxDbLite dependency
sudo  R -e 'library(BiocInstaller);biocLite(c("DBI","RSQLite","ensembldb","rtracklayer","stringdist","Biostrings","OrganismDbi","Rsamtools"),ask=FALSE)'
sudo apt-get install samtools -y




