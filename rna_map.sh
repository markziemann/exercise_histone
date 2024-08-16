#!/bin/bash

cd ../raw_data/rna/

IDX=../../ref/gencode.v46.transcripts.fa.gz.idx

for FQZ in *fq.gz ; do

  OUTDIR=$(echo $FQZ | sed 's/.fq.gz//')

  kallisto quant --single -l 150 -s 50 -t 8 -i $IDX -o $OUTDIR $FQZ

done

for TSV in */*abundance.tsv ; do
  NAME=$(echo $TSV | cut -d '/' -f1)
  cut -f1,4 $TSV | sed 1d | sed "s/^/${NAME}\t/"
done | gzip > 3col.tsv.gz

cd ../..
