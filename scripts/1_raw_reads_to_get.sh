#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 96:00:00
#SBATCH -J raw_reads
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

cd /crex/proj/snic2021-23-14/Xue/PRJNA369563
module load bioinfo-tools sratools/3.0.0
for i in $(seq 593 673);
do fastq-dump  --split-files --gzip GSM2474$i
done