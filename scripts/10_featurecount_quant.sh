#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 27:00:00
#SBATCH -J quantification
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools subread/2.0.3
for i in $(find /crex/proj/snic2021-23-14/Xue/PRJNA369563/star_out2/ -name "*Aligned.toTranscriptome.out.bam");
do
echo $i
a1=$(basename ${i/Aligned.toTranscriptome.out.bam/})
printf "this_is\t$a1\n"
featureCounts \
-a /crex/proj/snic2021-23-14/Xue/genome/gencode.v19.annotation.gtf \
-o /crex/proj/snic2021-23-14/Xue/PRJNA369563/featurecounts/${a1}_counts_output.txt \
-T 10 \
-p /crex/proj/snic2021-23-14/Xue/PRJNA369563/star_out2/${a1}Aligned.toTranscriptome.out.bam
done