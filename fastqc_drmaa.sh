#!/bin/bash
set -x

#Source modules for current shell
source $MODULESHOME/init/bash

#Loads specified  fastqc version
module load fastqc/$4

#Make assemblies output director if it doesn't exist
mkdir -p $3/$1

#Runs fastqc with $1=SAMPLE, $2=RAW_DATA_DIR, $3=QC_PATH, $4 VERSION $5 COMPRESSION

comp=$5
#create config file
if [ "$comp" = "GZIP" ]; then
	fastqc -o $3/$1 $2/$1.fastq.gz
else
	fastqc -o $3/$1 $2/$1.fastq
fi

exit 0
