#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 27:00:00
#SBATCH -J index
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools rsem/1.3.3
rsem-prepare-reference --gtf /crex/proj/snic2021-23-14/Xue/genome/h19.gtf \
/crex/proj/snic2021-23-14/Xue/genome/GRCh37.p13.genome.fa \
/crex/proj/snic2021-23-14/Xue/genome/ref_rsem