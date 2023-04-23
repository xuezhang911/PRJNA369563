#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 27:00:00
#SBATCH -J index
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools star/2.7.9a
STAR --runMode genomeGenerate \
--runThreadN 6 \
--genomeDir /crex/proj/snic2021-23-14/Xue/spike_in \
--genomeFastaFiles /crex/proj/snic2021-23-14/Xue/spike_in/ERCC92.fa \
--sjdbGTFfile /crex/proj/snic2021-23-14/Xue/spike_in/ERCC92.gtf

