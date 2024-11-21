#!/bin/bash

cd raw_data/h3k9ac/

IDX=../../ref/GRCh38.primary_assembly.genome.fa.gz

TSS=../../ref/tss_merged.saf

for FQZ in *fq.gz ; do

  BAM=$(echo $FQZ | sed 's/.fq.gz/.bam/')

  bwa mem -t 16 $IDX $FQZ | samtools sort -@8 -o $BAM -

done

for BAM in *bam ; do samtools index $BAM & done ; wait

featureCounts -f -O -T 16 -d 30 -Q 10 -F SAF -a $TSS -o h3k9_counts.tsv *bam

cd ../h3k36ac/

for FQZ in *fq.gz ; do

  BAM=$(echo $FQZ | sed 's/.fq.gz/.bam/')

  bwa mem -t 16 $IDX $FQZ | samtools sort -@8 -o $BAM -

done

for BAM in *bam ; do samtools index $BAM & done ; wait

featureCounts -f -O -T 16 -d 30 -Q 10 -F SAF -a $TSS -o h3k36_counts.tsv *bam

for BAM in *bam ; do
  POS=$BAM.pos
  bamToBed -i $BAM \
  | awk ' {OFS="\t"} { if ($6=="+") { print $1,$2+99,$2+100 } else { print $1,$3-100,$3-99 } } ' \
  | awk '$2>0 && $3>0' \
  | bedtools sort \
  |  bedtools closest -D b -a - -b ../../ref/tss1bp.bed > $POS
done
