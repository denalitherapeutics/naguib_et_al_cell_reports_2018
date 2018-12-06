FROM rocker/rstudio:3.5.0
MAINTAINER Thomas Sandmann, sandmann@dnli.com
RUN apt-get update && apt-get install -y \
    libbz2-dev \
    libicu-dev \
    libpcre3-dev \
    liblzma-dev \
    zlib1g-dev
USER rstudio
COPY ["manifest.txt", "install.R", "install_supt4h1.R", "/tmp/"]
RUN Rscript -e "install.packages('switchr')"
RUN Rscript "/tmp/install.R"
RUN Rscript "/tmp/install_supt4h1.R"
USER root
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
