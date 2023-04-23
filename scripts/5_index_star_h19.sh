#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t  27:00:00
#SBATCH -J index
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools star/2.7.9a
STAR --runMode genomeGenerate \
--runThreadN 6 \
--genomeDir /crex/proj/snic2021-23-14/Xue/Genome/human \
--genomeFastaFiles /crex/proj/snic2021-23-14/Xue/genome/GRCh37.p13.genome.fa \
--sjdbGTFfile /crex/proj/snic2021-23-14/Xue/genome/gencode.v19.annotation.gtf \
--sjdbOverhang 99 \
--sjdbGTFtagExonParentTranscript Parent