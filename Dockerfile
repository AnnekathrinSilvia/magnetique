FROM openanalytics/r-base

MAINTAINER Thiago Britto-Borges "thiago.brittoborges@uni-heidelberg.de"

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*

# install r enviroment 
ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY renv.lock ./
RUN R -e 'renv::restore()'

# set options
COPY Rprofile.site /usr/lib/R/etc/

# setup app
COPY . /root/magnetique/

# setup data
WORKDIR /root/magnetique/
RUN download_data.sh

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('root/magnetique/app.R')"]