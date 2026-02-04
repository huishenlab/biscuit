mkdir -p ref

wget -qO - https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/chr1.fa.gz | \
gunzip -c | \
bgzip -c - > ref/chr1.fa.gz

samtools faidx ref/chr1.fa.gz
