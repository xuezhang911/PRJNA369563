#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 72:00:00
#SBATCH -J RNAseq_alignment
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools star/2.7.9a
cd /crex/proj/snic2021-23-14/Xue/PRJNA369563/clean
for i in $(find .  -name "*_1.fq.gz");
do
echo $i
a1=$(basename ${i/_1.fq.gz/})
printf "this_is\t$a1\n"
STAR --genomeDir /crex/proj/snic2021-23-14/Xue/Genome/human \
--runThreadN 6 \
--readFilesIn $i ./${a1}_2.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /crex/proj/snic2021-23-14/Xue/PRJNA369563/star_out2/$a1 \
--limitBAMsortRAM 16000000000 \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--quantMode TranscriptomeSAM GeneCounts
done