FROM rocker/rstudio:4.0.3

RUN apt-get update && apt-get install -y \
    libxml2 \
    libglpk-dev
    

RUN install2.r \
    data.table \
    phytools \
    phangorn \
    corHMM \
    parallel \
    here