#!/bin/bash

NUMPROC=48
MEMORY=350
PROJ=${PWD##*/}
#USR_KMER="71,81,91,99,121,127"

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
if test -n "$(find ./raw -name '*R2*' -print -quit)"; then
    parallel -j $NUMPROC trim_galore --paired --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" ` :::+ ` find  ./raw  -name "*_R2*.f*.gz" `
    FWDS=`find ./trimming -name "*R1*.gz" | sort | uniq`
    REVS=(` find ./trimming -name "*R2*.gz" | sort | uniq `)
else
   parallel -j $NUMPROC trim_galore --illumina --fastqc -o trimming/ ::: ` find  ./raw  -name "*_R1*.f*.gz" `
   FWDS=`find ./trimming -name "*R1*.gz" | sort | uniq`
fi 
 
#####SPADES ASSEMBLY#####
if test -n "$(find ./raw -name '*R2*' -print -quit)"; then
    echo reverse reads found. proceeding with paired end assembly using $NUMPROC cores and $MEMORY GB memory
    for file1 in $FWDS; do
        file2=${file1/R1_001_val_1/R2_001_val_2}
        tempname="${file1/.\/trimming\//}"
        outname="${tempname/_*_R*.fq.gz/}_output"
        echo $outname
        #spades.py --careful -k 71,81,91,99,121,127 -1 $file1 -2 $file2 -t $NUMPROC -m $MEMORY -o ./assembly/${outname}_usr_kmer 
        spades.py --careful -1 $file1 -2 $file2 -t $NUMPROC -m $MEMORY -o ./assembly/${outname}_default_kmer
    done
else
    echo no reverse reads found. proceeding with single end assebmly using $NUMPROC cores and $MEMORY GB memory
    for file1 in $FWDS; do
        tempname="${file1/.\/trimming\//}"
        outname="${tempname/_*_R*.fq.gz/}_output"
        #spades.py --careful -k 71,81,91,99,121,127 -1 $file1 -t $NUMPROC -m $MEMORY -o ./assembly/${outname}_usr_kmer
        spades.py --careful -1 $file1 -t $NUMPROC -m $MEMORY -o ./assembly/${outname}_default_kmer
    done
fi

#####QUAST#####-----------------------------------------------------------
#per KMER
for genome in `find ./assembly -maxdepth 1 -name "*kmer" -type d | cut -d "/" -f3`; do
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
#all_user=(` find ./quast/*usr*contigs.fasta  -type f `)
#quast.py -e -t $NUMPROC -m $MEMORY  "${all_user[@]}" -o ./quast/usr_spades/

#spades default settings
#all default KMER selection

all_default=(` find ./quast/*default*contigs.fasta  -type f `)
quast.py -e -t $NUMPROC -m $MEMORY "${all_default[@]}" -o ./quast/default_spades/  

#########STOP######-----------------------------------------------------------
#before proceeding: evaluate your quast results. Choose desired de novo, copy to ./assemly with a name to include *chosen.fasta

