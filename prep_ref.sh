#!/bin/bash

mkdir -p ref \
  && cd ref \
  && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz \
  && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz \
  && wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz \
  && kallisto index -i gencode.v46.transcripts.fa.gz.idx gencode.v46.transcripts.fa.gz \
  && bwa index GRCh38.primary_assembly.genome.fa.gz



