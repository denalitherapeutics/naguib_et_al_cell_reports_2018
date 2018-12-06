options(repos = c(CRAN = "http://cran.rstudio.com"),
        Ncpus = parallel::detectCores())
library(switchr)
source("https://bioconductor.org/biocLite.R")
manifest <- loadManifest('/tmp/manifest.txt')  # within docker container
switchTo("supt4h1", seed = manifest)
