
#!/bin/bash -l
module load bioinfo-tools rsem/1.3.3
rsem-prepare-reference --gtf /crex/proj/snic2021-23-14/Xue/spike_in/ERCC92.gtf \
/crex/proj/snic2021-23-14/Xue/spike_in/ERCC92.fa \
/crex/proj/snic2021-23-14/Xue/spike_in/ref_rsem_ERCC