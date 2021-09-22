# Evolution of Photosymbiosis

This repo contains data and scripts necessary to replicate the results from Gault et al. (2021) and Gault (2019). 

### Running analysis in Docker container  

We provide three different docker images in the `dockerfiles` directory that can be used to replicate the analysis in different environments. `rstudio-v3-5-3` contains a dockerfile that allows the user to run the analysis in Rstudio. `rver-v3-5-3` contains a dockerfile that allows the user to run the analysis in command line R. `rscript-v3-5-3` contains a dockerfile that allows the user to batch run the analysis on a cluster by using an entrypoint to R. More details on the three images can be found below. All docker images can be built directly from the dockerfiles or can be pulled from Dockerhub.  


#### rstudio-v3-5-3

The Dockerfile contains instructions to build an image with Rstudio with R version 3.5.3. All packages and libraries necessary for the analysis are included in the image. The image can be downloaded directly from Dockerhub or built from the Dockerfile. 

Pull the image directly from Dockerhub via the following:

`docker pull jagault/evolution-photosymbiosis:rstudiov3.5.3`

Or, from within the directory containing the Dockerfile, build the image via the following:

`docker build . -t jagault/evolution-photosymbiosis:rstudiov3.5.3`  

Once the image is built, it can be run within the main repo directory `/evolution-photosymbiosis` via the following command:

`docker run -d -p 8787:8787 -v /path/to/repo:/home/rstudio -e PASSWORD=yourpasswordhere jagault/evolution-photosymbiosis:rstudiov3.5.3`  

The use of `-d` makes the container run in the background. `-p 8787:8787` opens a port through which rstudio can be accessed through an internet browser. On a linux machine, rstudio can be accessed via `http://localhost:8787`.  

If you are using a Windows or Mac machine, you can access rstudio via `http://yourip:8787`. Your ip address should be displayed when opening up the Docker Quickstart terminal.  

Rstudio no longer supports the use of the default password when none is specified. The user must now specify a password following the `-e` tag. 

The use of `-v /path/to/repo:/home/rstudio` will set up a bind mount between the project repo and rstudio within the container. All files within the repo will be available within Rstudio and saved changes will persist outside of the container.  

Because the corHMM analyses take a very long time to run, it is recommended that they are not run through the Rstudio console. Rather it is better to source the scripts to R through the terminal. This way the analyses can be monitored via a log file and the R session within Rstudio can be refreshed freely. The analyses can also be run in the background through the terminal via `tmux` or `screen`. Neither are included with the container but can be installed within it if desired. R and all files within the repo are available through the terminal window included with Rstudio. This is the method we used for running all analyses in Gault 2019. 

#### rver-v3-5-3

The Dockerfile contains instructions to build an image with command line base R version 3.5.3. All packages and libraries necessary for the analysis are included in the image. The image can be downloaded directly from Dockerhub or built from the Dockerfile. 

Pull the image directly from Dockerhub via the following:

`docker pull jagault/evolution-photosymbiosis:rver3.5.3`

Or, from within the directory containing the Dockerfile, build the image via the following:

`docker build . -t jagault/evolution-photosymbiosis:rver3.5.3`  

Once the image is built, it can be run within the main repo directory `/evolution-photosymbiosis` via the following command:

`docker run -ti jagault/evolution-photosymbiosis:rver3.5.3`  

#### rscriptv3-5-3

This image is built with an entrypoint to Rscript. Running this container automatically calls Rscript which will run a specified R script with optional arguments. The R script must be in the same directory from which you are running the container or a path must be provided. 

The image can be pulled from Docker Hub or built with the following commands. 

`docker pull jagault/evolution-photosymbiosis:rscriptv3.5.3`

`docker build . -t jagault/evolution-photosymbiosis:rscriptv3.5.3` 

Below is an example command for running `stree-1rate.R`. This script takes the posterior distribution of supertrees and associated traits and runs the corHMM analysis with a single rate category. The ouput is a `.rds` file that contains the results. This script takes two additional arguments. The first argument designates which tree to select (in this case `1`) and the second designates which replicate which will be used to name the output file (in this case `2`). The output file in this case will be named `stree-1rate-t1-r2.rds`. 

`docker run -ti -v /directory/containing/script/and/data:/home/analysis jagault/evolution-photosymbiosis:rscriptv3.5.3 stree-1rate.R 1 2`

It is easy to loop through trees with any desired number of replicates by replacing the arguments `1` and `2`. The analysis can thus be submitted as batch jobs to a cluster. This is the method we used to run all analyses for Gault et al. 2021. 

### Structure of repository 

This repo was originally for the analysis for my masters thesis, Gault 2019. For the publication Gault et al. 2021 we expanded and reran the analyses. All scripts for these analyses can be found in analysis_rerun. More details on this directory can be found below. 

`/R` contains helper functions used through the scripts in all anlayses.  

`/data` contains the phylogenetic trees and trait data used for all analyses. 

`/analysis` contains the scripts to perform the analyses with associated results for Gault 2019.  
  * `/stree_corHMM` contains script to run the corHMM analysis across the supertree posterior distribution.
  * `/stree_asr` contains the script `stree_asr.R` which performs ancestral state reconstruction based on the rates from the corHMM analysis and the script `stree_nodeframe.R` which summarizes the results of the ancestral state reconstruction across the supertree posterior distribution.  
  * `/mtree_corHMM` contains script to run the corHMM analysis across the molecular tree posterior distribution.
  * `/mtree_asr` contains the script `mtree_asr.R` which performs ancestral state reconstruction based on the rates from the corHMM analysis and the script `mtree_nodeframe.R` which summarizes the results of the ancestral state reconstruction across the molecular posterior distribution.  
  * `/mtree_stree_facultative_removed` contains scripts to rerun the above analyses across the supertree and molecular tree posterior distribution with facultative species removed.  
  
`/figures` contains the scripts necessary to reproduce the raw figures from Gault 2019. 

`/analysis_rerun` contains all scripts and results from Gault et al. 2021. 
   * `/stree_corHMM` contains scripts that were used to run the analysis using the method described above for the `rscriptv3.5.3` docker image. Also contains scripts to process the output data. Because the replicates for each tree were used to thoroughly sample the likelihood surface, we only wanted to keep the result with the highest likelihood for each tree. This is done with the process scripts. 
   * `/stree_asr` contains scripts to perform ancestral state reconstructions across the individual trees for each model. Also contains scripts to summarize the results on a single consensus tree for selected nodes. 
   * /mtree_corHMM and /mtree_asr contain the same scripts as above but for running the analysis across the molecular trees. 
   * /results contains scripts for creating figures, summarizing rates, and counting number of gains/losses for the supertrees and molecular trees. 

