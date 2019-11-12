#!/bin/bash

NUMPROC=30
MEMORY=225
PROJ=${PWD##*/}
USR_KMER="71,81,91,99,121,127"

#suggested file structure
# project
# |____raw
# |   |____genome1
# |   |____genome2
# |   |____genome3
# |   |____...
# |____fastqc
# |____assembly
# |____quast
# |____etc
# |____...
#ddefine resources to use

#Trim galore: 
#trim adapter seqs, overrep seqs, qual < 20 (also dedup w/ erika's suggestion?)
if [ ! -d "trimming" ]; then mkdir trimming; fi
if [ ! -d "assembly" ]; then mkdir assembly; fi
if [ ! -d "quast" ]; then mkdir quast; fi
      
#####TRIMMING#####-----------------------------------------------------------
if test -n "$(find ./raw -maxdepth 1 -name '*R2*' -print -quit)"; then
    parallel -j  $NUMPROC trim_galore --paired --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" ` :::+ ` find  ./raw  -name "*_R2*.f*.gz" `
    find ./raw -name "*R1*.gz"* | sort | uniq > fwds
    find ./raw -name "*R2*.gz"* | sort | uniq > revs
else
   parallel  -j $NUMPROC trim_galore --paired --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" `
fi 
 
#####SPADES ASSEMBLY#####
if test -n "$(find ./raw -maxdepth 1 -name '*R2*' -print -quit)"; then
    echo reverse reads found. proceeding with paired end assembly
    #setspades.R
    echo "if (\"here\" %in% row.names(installed.packages())){" > spades_yaml.R
    echo "library(here)" >> spades_yaml.R
    echo "}  else {" >> spades_yaml.R
    echo "install.packages(\"here\")" >> spades_yaml.R
    echo "library(here)}" >> spades_yaml.R
    echo "workingdir <- paste0(here(),\"/\")" >> spades_yaml.R
    echo "fwds <- read.table(paste0(workingdir,\"/fwds\"))" >> spades_yaml.R
    echo "revs <- read.table(paste0(workingdir,\"/revs\"))" >> spades_yaml.R
    echo "setwd(workingdir)" >> spades_yaml.R
    echo "write(paste0('[ \n',  '   { \n', '     orientation: \"fr\", \n', '     type: \"paired-end\", \n', '     right reads: ['), file = \"libraries.yaml\", append = F)" >> spades_yaml.R
    echo "for(i in 1:nrow(fwds)){write(paste0('       \"',fwds[i,1],'\",'), file = \"libraries.yaml\", append = T)}" >> spades_yaml.R
    echo "write(paste0('       ], \n', '       left reads: ['), file = \"libraries.yaml\", append = T)" >> spades_yaml.R
    echo "for(i in 1:nrow(revs)){write(paste0('       \"',revs[i,1],'\",'), file = \"libraries.yaml\", append = T)}" >> spades_yaml.R                   
    echo "write(paste0('       ] \n', '   } \n', ']'), file = \"libraries.yaml\", append = T)" >> spades_yaml.R

    chmod 755 spades_yaml.R
    Rscript spades_yaml.R

    #paired-end assembly
    spades.py --dataset libraries.yaml -k "${USR_KMER[@]}" --careful -t 30 -m 300 -o ./assembly/$PROJ_usr_kmer
    spades.py --dataset libraries.yaml --careful -t $NUMPROC -m $MEMORY -o ./assembly/$PROJ_default_kmer
else
    echo no reverse reads found. proceeding with single end assebmly
    #single end assembly
    pools=(` find ./trimming -name "*R1*.gz" | sort | uniq | cat `) 
    libnum=( `seq 1 "${#pools[@]}"` )
    libnum=("${libnum[@]/#/--s}")
    unset reads
    for (( i=0; i<${#libnum[*]}; ++i)); do reads+=( ${libnum[$i]} ${pools[$i]} ); done
    #kmers="$(seq -s ',' 21 2 127)"
    #spades.py "${reads[@]}" -k "$kmers" --careful -t $NUMPROC -m $MEMORY -o ./assembly/ipyRAD_all_kmer_output

    spades.py "${reads[@]}" -k "${USR_KMER[@]}" --careful -t $NUMPROC -m $MEMORY -o ./assembly/${PROJ}usr_kmer
    spades.py "${reads[@]}" --careful -t $NUMPROC -m $MEMORY -o ./assembly/${PROJ}_default_kmer
fi

#####QUAST#####-----------------------------------------------------------
#per KMER
for genome in `find ./assembly -maxdepth 1 -name *$PROJ* -type d | cut -d "/" -f3`; do
    KMERS=(` find ./assembly/$genome/ -maxdepth 1 -name "K*" | cut -d "/" -f4 | sort | uniq `)
    cp ./assembly/${genome}/contigs.fasta ./quast/"${genome}_contigs.fasta"
         for KMER in `find ./assembly/$genome/ -maxdepth 1 -name "K*" | cut -d "/" -f4 | sort | uniq`; do
             cp ./assembly/${genome}/${KMER}/final_contigs.fasta ./quast/${genome}_${KMER}.fasta
        done
    by_kmer=("${KMERS[@]/#/./quast/${genome}_}")
    by_kmer_files=("${by_kmer[@]/%/.fasta}")
    quast.py -e -t $NUMPROC -m $MEMORY ./quast/${genome}_contigs.fasta "${by_kmer_files[@]}" -o ./quast/${genome}_results/
done

#all default KMER selection
all_user=(` find ./quast/*usr*contigs.fasta  -type f `)
quast.py -e -t $NUMPROC -m $MEMORY  "${all_user[@]}" -o ./quast/usr_spades/

#spades default settings
#all default KMER selection

all_default=(` find ./quast/*default*contigs.fasta  -type f `)
quast.py -e -t $NUMPROC -m $MEMORY "${all_default[@]}" -o ./quast/default_spades/  

#########STOP######-----------------------------------------------------------
#before proceeding: evaluate your quast results. Choose desired de novo, copy to ./assemly with a name to include *chosen.fasta

