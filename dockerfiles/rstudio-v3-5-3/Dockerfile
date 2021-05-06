FROM rocker/rstudio:3.5.3

RUN apt-get update && apt-get install -y \
    libmagick++-dev \
    libmpfr-dev

RUN install2.r \
    data.table \
    phangorn \
    parallel \
    here \
    phytools \
    corHMM