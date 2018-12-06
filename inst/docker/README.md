This directory contains the necessary instructions to create a docker image
with versioned R packages that can be used to reproduce the analyses published
in Naguib et al, Cell Reports, 2019.

# Retrieving the docker image from dockerhub

The final, precompiled docker image is available via
[dockerhub](https://hub.docker.com/)
and can be retrieved with the following shell command:

```
docker pull denalitherapeutics/supt4h1
```

# Rebuilding the docker image from scratch

Alternatively, the docker image can be rebuilt using the files in this
directory. For the original analysis the
[switchr](https://cran.r-project.org/web/packages/switchr/index.html)
R package was used to document the exact version and installation source
for each additional R package.

During the `docker build` process (see below), the `switchr` manifest is
used to recreate the same R library within the docker container.

To start the docker build process, execute the following command in the home
directory of this R package:

```bash
docker build -t supt4h1 inst/docker
```

# Starting a docker container

The docker image is built on the `rocker/rstudio:3.5` image, which contains
the Rstudio IDE.

To start a container running Rstudio and share the contents of the current
working directory in the /home/rstudio/SUPT4H1 directory in the container,
execute

```bash
docker run --rm \
    -v $PWD:/home/rstudio/SUPT4H1 \
    -d -p 8787:8787 \
    supt4h1
```

browse to `http://localhost:8787` and log in as user `rstudio` (password:
rstudio).

# Activating the switchr library

To access the `switchr` library created from third-party packages, execute
the following code within you Rstudio session (within the container):

```r
library(switchr)
switchTo('supt4h1')
```
