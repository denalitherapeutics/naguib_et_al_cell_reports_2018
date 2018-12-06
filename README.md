# SUPT4H1 companion R package for Naguib et al, Cell Reports, 2018

This repository contains the `SUPT4H1` R package that can be used to reproduce
analyses and figures from Naguib *et al*, which has been accepted for
publication in
[Cell Reports](https://www.cell.com/cell-reports/home)
.

## Installation

To install this R package and the necessary dependencies, please use the
`BiocManager::install` command:

```r
install.packages(c("BiocManager", "remotes"))
library(BiocManager)
BiocManager::install("denalitherapeutics/naguib_et_al_cell_reports_2018",
                     dependencies = TRUE)
```

## R markdown documents

This repository's `inst/rmd` directory includes R markdown files that
reproduce the analysis results and figures of the paper. After the `SUPT4H1` R
package has been installed in your library, you can locate the location of the
`.Rmd` files on your system with the following command:

```r
dir(system.file("rmd", package = "SUPT4H1"))
```

and open e.g. the `figure_1.Rmd` file for editing from the R console:

```r
file.edit(system.file("rmd", "figure_1.Rmd", package = "SUPT4H1"))
```

Each of the R markdown documents can be `knit` into an html report, either
separately for each Figure or a single report with all results
(`all_figures.Rmd`).

## Reproducing the original R environment

The following sections outline how to reproduce the *exact* R environment used
for the original analysis.

For this analysis, a number of open-source packages have been used, e.g. the
`ggplot2` or `limma` packages that are available from
[CRAN](https://cran.r-project.org/)
and
[Bioconductor](https://www.bioconductor.org/)
respectively.

To ensure reproducibility, the version number and the provenance of all R
packages has been recorded in the `inst/docker/manifest.txt` file.

### Docker image

We provide a docker image that contains R 3.5.0, all third party R packages 
that were used and the 
[Rstudio](https://www.rstudio.com/products/rstudio/)
IDE. 

#### Using the precompiled docker image

The final, precompiled docker image is available via
[dockerhub](https://hub.docker.com/)
and can be retrieved with the following shell command:

```bash
docker pull denalitherapeutics/naguib_et_al_cell_reports_2018
```

To explore the code in this repository within the docker container, you can
clone the repository onto your local machine and share the directory with
a docker container, e.g. with the following shell commands:

```bash
git clone git@github.com:denalitherapeutics/naguib_et_al_cell_reports_2018.git SUPT4H1
cd SUPT4H1
docker run --rm -v $PWD:/home/rstudio/SUPT4H1 -d -p 8787:8787 denalitherapeutics/naguib_et_al_cell_reports_2018
```

You can access an Rstudio session that runs within the docker container by
pointing your web browser at http://localhost:8787.
(Username: `rstudio`, password: `rstudio`).

The required R packages are available in a separate `supt4h1` compute
environment that was created with the 
[switchr](https://cran.r-project.org/web/packages/switchr/index.html)
package. It can be activated at the beginning of any interactive R
session with the following command:

```r
library(switchr)
switchTo("supt4h1")
```

#### Recreating the `supt4h1` compute environment

The docker image was created with the help of the 
[switchr](https://github.com/gmbecker/switchr) R package because it allows
users to install specific versions of R packages. 

The `inst/docker` directory contains the `manifest.txt` file, which
tracks the provenances and versions of all R packages used for the original 
analysis.

If you would like to create the `supt4h1` environment on your local system,
please install the `switchr` package and point it at the `manifest.txt`
file.

```r
install.packages(c("BiocManager", "remotes", "switchr"))
library(switchr)
manifest_file <- system.file("docker", "manifest.txt", package = "SUPT4H1")
manifest <- loadManifest(manifest_file)
```

The `switchr` manifest contains information about all required R packages and
the specific version used in the analysis.

```
> manifest

A seeding manifest (SessionManifest object)

Describes a cohort of 89 package versions. 
89 packages are listed in the underlying package manifest

Package versions:
    name            version    
1   "AnnotationDbi" "1.42.1"   
2   "assertthat"    "0.2.0"    
3   "backports"     "1.1.2"    
4   "base64enc"     "0.1-3"    
... "..."           "..."      
85  "docopt"        "0.4.5"    
86  "littler"       "0.3.3"    
87  "RCurl"         "1.95-4.10"
88  "RJSONIO"       "1.3-0"    
89  "switchr"       "0.12.7"   
```

To create and populate the `supt4h1` environment (or switch to a previously
created one), please use the `switchTo` command:

```r
switchTo("supt4h1", seed = manifest)
```

We used the same code during the `docker build` process (see below), to
create the `~/.switchr/supt4h1` R library within the docker container.

If you would like to rebuild the `supt4h1` docker image from scratch, execute
the following command in the `inst/docker` directory of the `SUPT4H1` package:

```bash
docker build -t supt4h1 inst/docker
```
