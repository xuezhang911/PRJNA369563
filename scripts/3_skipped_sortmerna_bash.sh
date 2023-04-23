#!/bin/bash
set -eu
if [ $# -ne 2 ]; then
 echo "this scripts expects 2 arguments"
 exit 1
fi
if [ ! -f $1 ]; then
 echo "the forward file does not exist"
 exit 1
fi
if [ ! -f ${2} ]; then
 echo "the forward file does not exist"
 exit 1
fi
set -eu
a1=$(basename ${1/_1.fq.gz/})
if [ -d /crex/proj/snic2021-23-14/Xue/PRJNA369563/sortmerna/kvdb ]; then
 rm -rf /crex/proj/snic2021-23-14/Xue/PRJNA369563/sortmerna/kvdb/*
fi
module load bioinfo-tools SortMeRNA/4.3.3
sortmerna \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5.8s-database-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5s-database-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-16s-id95.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-23s-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-16s-id90.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-23s-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-18s-id95.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-28s-id98.fasta \
--reads $1 \
--reads $2 \
--workdir /crex/proj/snic2021-23-14/Xue/PRJNA369563/sortmerna --fastx \
--other /crex/proj/snic2021-23-14/Xue/PRJNA369563/sortmerna/$a1 \
--paired_in --threads 8 --out2
mv /crex/proj/snic2021-23-14/Xue/PRJNA369563/sortmerna/out/aligned.log /crex/proj/snic2021-23-14/Xue/PRJNA369563/sortmerna/out/${a1}-aligned.log


