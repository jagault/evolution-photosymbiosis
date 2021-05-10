FROM rocker/r-ver:3.5.3

RUN apt-get update && apt-get install -y \
    libmagick++-dev \
    libmpfr-dev \
    libcurl4-openssl-dev

RUN install2.r \
    data.table \
    phangorn \
    here \
    phytools \
    corHMM

WORKDIR /home/analysis

ENTRYPOINT ["Rscript"]