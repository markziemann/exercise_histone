#!/bin/bash
# bedtools
# pigz
# parallel


mkdir -p ref \
  && cd ref \
  && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz \
  && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz \
  && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz \
  && kallisto index -i gencode.v46.transcripts.fa.gz.idx gencode.v46.transcripts.fa.gz \
  && bwa index GRCh38.primary_assembly.genome.fa.gz

zcat gencode.v46.primary_assembly.annotation.gtf.gz \
| grep 'exon_number 1;' \
| awk '$7=="+"' \
| cut -d ';' -f1 \
| cut -f1,4,9 \
| sed 's/gene_id "//' \
| tr -d '"' \
| sort -u > tss.tsv

zcat gencode.v46.primary_assembly.annotation.gtf.gz \
| grep 'exon_number 1;' \
| awk '$7=="-"' \
| cut -d ';' -f1 \
| cut -f1,5,9 \
| sed 's/gene_id "//' \
| tr -d '"' \
| sort -u >> tss.tsv

# expand TSS windows with slop
gunzip -k GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa && rm  GRCh38.primary_assembly.genome.fa
cut -f-2 GRCh38.primary_assembly.genome.fa.fai > GRCh38.primary_assembly.genome.fa.g

awk '{OFS="\t"}{print $1,$2,$2,$3}' tss.tsv \
| bedtools slop -b 2000  -g GRCh38.primary_assembly.genome.fa.g \
| bedtools sort > tss.bed

# merge overlapping TSS promoters
for GENE in $(cut -f4 tss.bed | sort -u ) ; do
  grep -w $GENE tss.bed \
  | bedtools merge \
  | sed "s/$/\t${GENE}/"
done > tss_merged.bed


echo "GeneID Chr Start End Strand" | tr ' ' '\t' > tss_merged.saf
awk '{OFS="\t"}{print $4,$1,$2,$3,"."}' tss_merged.bed >> tss_merged.saf
