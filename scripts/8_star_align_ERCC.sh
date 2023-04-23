#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 39:00:00
#SBATCH -J spike_in_alignment
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools star/2.7.9a
cd /crex/proj/snic2021-23-14/Xue/PRJNA369563/star_out2
for i in $(find .  -name "*Unmapped.out.mate1");
do
echo $i ${i/Unmapped.out.mate1/Unmapped.out.mate2}
a1=$(basename ${i/Unmapped.out.mate1//})
printf "this_is\t$a1\n"
STAR --genomeDir /crex/proj/snic2021-23-14/Xue/spike_in \
--runThreadN 12 \
--readFilesIn $i ${i/Unmapped.out.mate1/Unmapped.out.mate2} \
--outFileNamePrefix /crex/proj/snic2021-23-14/Xue/PRJNA369563/star_ERCC/$a1 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--quantMode TranscriptomeSAM GeneCounts
done