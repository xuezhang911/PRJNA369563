#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 39:00:00
#SBATCH -J trimreads
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

cd /crex/proj/snic2021-23-14/Xue/PRJNA369563/fastq
module load bioinfo-tools TrimGalore/0.6.1 FastQC/0.11.9 MultiQC/1.12
for i in $(find . -name "*_1.fastq.gz");
do
a1=$(basename ${i/_1.fastq.gz/})
printf "this_is\t$a1\n"
trim_galore -q 35 --phred33 --stringency 4 --illumina —paired $i ${a1}_2.fastq.gz --gzip -o ../clean --basename $a1
done
fastqc ../clean/*.fq.gz
multiqc ../clean/


