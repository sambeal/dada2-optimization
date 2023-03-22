#!/bin/bash

cd /Users/samanthabeal/Documents/MSc/Bioinformatics/dada2-optimization


################ Input data ################

# list files in folder
# extract part of name before second underscore and find unique
ls input | cut -d_ -f1,2 | sort | uniq

# store list of names in bash variable and check contents
samples=$(ls input | cut -d_ -f1,2 | sort | uniq)
echo $samples

################ Primer removal ################


# count number of sequences across all files in folder
cd ../../output3/cutadapt
gzip *.fastq

gzcat *.fastq.gz | grep -c "^@M00" 

for s in $samples;
do
            echo "${s}_L001_R1_001.fastq.gz" \
            count=$(gunzip -c ${s}_L001_R1_001.fastq.gz | grep -c "^@M00" ) \
       >> SeqNumber_primerremoval.txt;
done

gunzip *.fastq.gz

#output-1 = nothing specified, background default parameters: 
#output-2 = cutadapt -e 0.2,  --discard-untrimmed,  --minimum-length 1:
#output-3 = cutadapt -e 0.1,  --discard-untrimmed,  --minimum-length 1:
