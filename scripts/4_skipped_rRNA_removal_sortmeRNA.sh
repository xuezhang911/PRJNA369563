#!/bin/bash -l
#SBATCH -A naiss2023-22-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 36:00:00
#SBATCH -J sortme
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

for i in $(find /crex/proj/snic2021-23-14/Xue/PRJNA369563/clean  -name "*_1.fq.gz");
do
echo $i
bash /crex/proj/snic2021-23-14/Xue/scripts/sortme_bash2.sh  $i ${i/_1.fq.gz/_2.fq.gz}
done