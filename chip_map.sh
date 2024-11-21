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

# TSS read density plot

# get read position relative to closest TSS - the goal is to obtain a list of canonical TSSs
for BAM in *bam ; do
  POS=$BAM.pos
  bamToBed -i $BAM \
  | awk ' {OFS="\t"} { if ($6=="+") { print $1,$2+99,$2+100 } else { print $1,$3-100,$3-99 } } ' \
  | awk '$2>0 && $3>0' \
  | bedtools sort \
  | bedtools closest -D b -a - -b ../../ref/tss1bp.bed > $POS
done

# Find out which TSS has most reads annotated in each sample
for POS in *pos ; do
  grep chr $POS \
  | cut -f 7 \
  | sort \
  | uniq -c \
  | sort -k1nr \
  | tr '.' '\t' \
  | awk '!arr[$2]++{print $2"."$3}' > $POS.tss
done

# Curate a list of canonical TSSs (most reads) for each gene
awk '{OFS="\t"}{print $0,$0}' *tss  \
| sed 's#\.#\t#' \
| sort \
| uniq -c \
| sort -k1nr \
| awk '!arr[$2]++ {print $4}' > maintss.txt

# Get the BED file for next step
grep -wFf maintss.txt ../../ref/tss1bp.bed > maintss.bed


# Now finally get the read position relative to the canonical TSS
for BAM in *bam ; do
  POS2=$BAM.pos2
  bamToBed -i $BAM \
  | awk ' {OFS="\t"} { if ($6=="+") { print $1,$2+99,$2+100 } else { print $1,$3-100,$3-99 } } ' \
  | awk '$2>0 && $3>0' \
  | bedtools sort \
  | bedtools closest -D b -a - -b ../../ref/tss1bp.bed > $POS2
done
