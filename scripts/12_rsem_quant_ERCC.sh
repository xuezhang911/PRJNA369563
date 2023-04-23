#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 27:00:00
#SBATCH -J quantification
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools rsem/1.3.3
for i in $(find /crex/proj/snic2021-23-14/Xue/PRJNA369563/star_ERCC/ -name "*Aligned.toTranscriptome.out.bam");
do
echo $i
a1=$(basename ${i/Aligned.toTranscriptome.out.bam/})
printf "this_is\t$a1\n"
rsem-calculate-expression \
--paired-end \
--no-bam-output \
--alignments -p 16 -q /crex/proj/snic2021-23-14/Xue/PRJNA369563/star_ERCC/${a1}Aligned.toTranscriptome.out.bam \
/crex/proj/snic2021-23-14/Xue/spike_in/ref_rsem_ERCC \
/crex/proj/snic2021-23-14/Xue/PRJNA369563/rsem_ERCC_2/${a1}_rsem
done
