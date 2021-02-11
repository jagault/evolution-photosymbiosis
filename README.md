# Evolution of Photosymbiosis

This repo contains data and scripts necessary to replicate the results from our publication.  

### Running analysis in Docker container  

The Dockerfile contains instructions to build an image with Rstudio with R version 4.0.3. All packages and libraries necessary for the analysis are included in the image. The image can be downloaded directly from Dockerhub or built from the Dockerfile. 

Pull the image directly from Dockerhub via the following:

`docker pull jagault/evolution-photosymbiosis`

Or build the image via the following:

`docker build -t jagault/evolution-photosymbiosis .`  

Once the image is built, it can be run via the following command:

`docker run -d -p 8787:8787 -v /path/to/repo:/home/rstudio -e PASSWORD=yourpasswordhere jagault/evolution-photosymbiosis`  

The use of `-d` makes the container run in the background. `-p 8787:8787` opens a port through which rstudio can be accessed through an internet browser. On a linux machine, rstudio can be accessed via `http://localhost:8787`.  

If you are using a Windows or Mac machine, you can access rstudio via `http://yourip:8787`. Your ip address should be displayed when opening up the Docker Quickstart terminal.  

Rstudio no longer supports the use of the default password when none is specified. The user must now specify a password following the `-e` tag. 

The use of `-v /path/to/repo:/home/rstudio` will set up a bind mount between the project repo and rstudio within the container. All files within the repo will be available within Rstudio and saved changes will persist outside of the container.  

Because the corHMM analyses take a very long time to run, it is recommended that they are not run through the Rstudio console. Rather it is better to source the scripts to R through the terminal. This way the analyses can be monitored via a log file and the R session within Rstudio can be refreshed freely. The analyses can also be run in the background through the terminal via `tmux` or `screen`. Neither are included with the container but can be installed within it if desired. R and all files within the repo are available through the terminal window included with Rstudio.