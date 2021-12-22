# This first stage is from openanalytics
# https://github.com/openanalytics/r-base/blob/b90354b20a1ddd392958bacfa500237ea197907d/Dockerfile
# I have update the R version
FROM ubuntu:20.04

LABEL MAINTAINER Thiago Britto-Borges "thiago.brittoborges@uni-heidelberg.de"

# Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
RUN useradd docker \
    && mkdir /home/docker \
    && chown docker:docker /home/docker \
    && addgroup docker staff

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    apt-utils \
    ed \
    less \
    locales \
    vim-tiny \
    wget \
    ca-certificates \
    apt-transport-https \
    gsfonts \
    gnupg2 \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" > /etc/apt/sources.list.d/cran.list

# note the proxy for gpg
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

ENV R_BASE_VERSION 4.1.1

# Now install R and littler, and create a link for littler in /usr/local/bin
# Also set a default CRAN repo, and make sure littler knows about it too
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    littler \
    r-cran-littler \
    r-base=${R_BASE_VERSION}* \
    r-base-core=${R_BASE_VERSION}* \
    r-base-dev=${R_BASE_VERSION}* \
    r-recommended=${R_BASE_VERSION}* \
    && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
    && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
    && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*



RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libfftw3-dev \
    && rm -rf /var/lib/apt/lists/*

# install renv 
ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# set options
COPY Rprofile.site /usr/lib/R/etc/

# setup app
COPY . /root/magnetique/
WORKDIR /root/magnetique/

# install R environment
RUN R -e 'renv::restore()'
RUN R -e 'renv::install("markdown")'
RUN R -e 'renv::install("future")'


EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('root/magnetique/app.R')"]
