#!/bin/bash

cd raw_data/h3k9ac/

IDX=../../ref/GRCh38.primary_assembly.genome.fa.gz

for FQZ in *fq.gz ; do

  BAM=$(echo $FQZ | sed 's/.fq.gz/.bam/')

  bwa mem -t 16 $IDX $FQZ | samtools sort -@8 -o $BAM -

done

for BAM in *bam ; do samtools index $BAM & done ; wait

cd ../h3k36ac/

for FQZ in *fq.gz ; do

  BAM=$(echo $FQZ | sed 's/.fq.gz/.bam/')

  bwa mem -t 16 $IDX $FQZ | samtools sort -@8 -o $BAM -

done

for BAM in *bam ; do samtools index $BAM & done ; wait
