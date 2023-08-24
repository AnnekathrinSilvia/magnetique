FROM shinyverse_criu_base:latest

WORKDIR /app

# setup app
COPY app.R /app
COPY utils.R /app
COPY database.R /app
COPY www/ /app/www/
COPY tours/ /app/tours/
COPY renv.lock /app/renv.lock

COPY deploy/rstart.R rstart.R
COPY deploy/restore.sh /restore.sh
COPY deploy/launch_shiny_app.sh launch_shiny_app.sh

USER root
RUN mkdir -p /criu_dumps \
    && chmod 755 /restore.sh launch_shiny_app.sh  \
    && mkdir -p /data 

EXPOSE 3838

CMD /bin/bash /restore.sh
