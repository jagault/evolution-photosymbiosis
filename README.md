# Evolution of Photosymbiosis

This repo contains data and scripts necessary to replicate the results from our publication.  

Build the image via the following:

`docker build -t jagault/evolution-photosymbiosis .`

Or pull the image directly from dockerhub via the following:

`docker pull jagault/evolution-photosymbiosis`

Use the following to run the image in a container. 

`docker run -p 8787:8787 -v /home/jordan/git_repos/evolution-photosymbiosis:/home/rstudio -e PASSWORD=test jagault/evolution-photosymbiosis`