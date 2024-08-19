FROM bioconductor/bioconductor_docker:RELEASE_3_19

# Update apt-get

RUN apt-get update \
        && apt upgrade -y \
        && apt install -y nano git libncurses-dev xorg openbox \
        && apt install -y samtools bwa kallisto bedtools parallel pigz \
        && pip install magic-wormhole           \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install CRAN packages

RUN Rscript -e 'install.packages(c("zoo","tidyverse","reshape2","gplots","MASS","eulerr","kableExtra","vioplot","pkgload","beeswarm"))'

# Install bioconductor packages

RUN Rscript -e 'BiocManager::install(c("DESeq2","fgsea","limma","topconfects","mitch"))'

# Install featurecounts
RUN mkdir /app \
  && cd /app \
  && wget -O download.tar.gz https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz/download \
  && tar xf download.tar.gz \
  && cp -r subread-2.0.6-Linux-x86_64/bin/* /usr/local/bin \
  && cd / 

# Install gtftools
RUN cd /app \
  && wget https://www.genemine.org/codes/GTFtools_0.9.0.zip \
  && unzip GTFtools_0.9.0.zip \
  && echo 'alias gtftools="python3 /app/GTFtools_0.9.0/gtftools.py"' >> /root/.bashrc \
  && . /root/.bashrc \
  && cd /

# get a clone of the code repo
RUN git clone https://github.com/markziemann/exercise_histone.git

# Set the container working directory

ENV DIRPATH /exercise_histone

#COPY raw_data $DIRPATH

WORKDIR $DIRPATH
